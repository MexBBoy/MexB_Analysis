#!/usr/bin/env python3
"""Stage 7 - tunnel analysis for MexB.

Bottleneck radius by threshold connectivity on a clearance grid, widest-path
tracing, constriction lining residues, and a pseudo-atom trace with the local
radius in the B-factor column.

Method
------
A regular grid is built over the whole trimer. For each voxel the *clearance*
is the distance to the nearest atomic van der Waals surface, i.e.
min_a(|v - x_a| - r_a). A tunnel of radius R exists between the seed and bulk
solvent iff the seed voxel and the grid boundary lie in one connected
component of {clearance >= R}. The bottleneck is the largest such R, found by
bisection; the path is then the shortest route through that component.

The calculation is run on the trimer. An isolated protomer has an artificially
open protomer-protomer interface and the path finder escapes through it.
"""
from __future__ import annotations

import argparse
import heapq
import os
import sys

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import (  # noqa: E402
    CXDIR, DBP, DETERGENTS, PBP, SUBDOMAINS, TABLES, Structure, centroid,
    coords, fmt, load_structures, vdw, write_csv,
)

STEP = 0.6          # grid spacing, angstrom
MARGIN = 6.0        # padding around the trimer bounding box
R_LO, R_HI = 0.4, 6.0
BISECT_TOL = 0.005


# ------------------------------------------------------------------- grid

class ClearanceGrid:
    def __init__(self, atoms, step=STEP, margin=MARGIN, verbose=True):
        self.step = step
        X = coords(atoms)
        self.radii = np.array([vdw(a.element) for a in atoms])
        self.origin = X.min(axis=0) - margin
        hi = X.max(axis=0) + margin
        self.shape = tuple(int(np.ceil((hi[i] - self.origin[i]) / step)) + 1
                           for i in range(3))
        n = int(np.prod(self.shape))
        if verbose:
            print(f"    grid {self.shape} = {n/1e6:.1f} M voxels "
                  f"@ {step} A")
        gx = self.origin[0] + step * np.arange(self.shape[0])
        gy = self.origin[1] + step * np.arange(self.shape[1])
        gz = self.origin[2] + step * np.arange(self.shape[2])
        tree = cKDTree(X)
        # k nearest centres; the surface-nearest atom is always among them
        # because the vdW radii here span only 1.20-1.80 A.
        k = 12
        clear = np.empty(n, dtype=np.float32)
        # chunk over i-slices to keep peak memory bounded
        per_slice = self.shape[1] * self.shape[2]
        chunk = max(1, int(2_000_000 // per_slice))
        Y, Z = np.meshgrid(gy, gz, indexing="ij")
        plane = np.stack([Y.ravel(), Z.ravel()], axis=-1)
        for i0 in range(0, self.shape[0], chunk):
            i1 = min(self.shape[0], i0 + chunk)
            xs = gx[i0:i1]
            pts = np.empty(((i1 - i0) * per_slice, 3), dtype=np.float64)
            for j, xv in enumerate(xs):
                sl = slice(j * per_slice, (j + 1) * per_slice)
                pts[sl, 0] = xv
                pts[sl, 1:] = plane
            d, idx = tree.query(pts, k=k, workers=-1)
            clear[i0 * per_slice:i1 * per_slice] = (
                d - self.radii[idx]).min(axis=1)
            del pts, d, idx
        self.clearance = clear.reshape(self.shape)
        self.atoms = atoms
        self.X = X

    def index_of(self, point):
        i = np.round((np.asarray(point) - self.origin) / self.step)
        i = np.clip(i, 0, np.array(self.shape) - 1).astype(int)
        return tuple(i)

    def point_of(self, idx):
        return self.origin + self.step * np.asarray(idx, dtype=float)

    def free_seed(self, point, min_clear=1.2, search=8.0):
        """Nearest voxel to `point` with clearance >= min_clear."""
        c = self.index_of(point)
        rad = int(np.ceil(search / self.step))
        sl = tuple(slice(max(0, c[i] - rad),
                         min(self.shape[i], c[i] + rad + 1)) for i in range(3))
        sub = self.clearance[sl]
        ok = np.argwhere(sub >= min_clear)
        if len(ok) == 0:
            return None, None
        off = np.array([sl[i].start for i in range(3)])
        cand = ok + off
        d = np.linalg.norm((cand - np.array(c)) * self.step, axis=1)
        best = cand[int(np.argmin(d))]
        return tuple(best), float(self.clearance[tuple(best)])


def boundary_mask(shape):
    m = np.zeros(shape, dtype=bool)
    m[0, :, :] = m[-1, :, :] = True
    m[:, 0, :] = m[:, -1, :] = True
    m[:, :, 0] = m[:, :, -1] = True
    return m


STRUCT3 = ndimage.generate_binary_structure(3, 1)   # 6-connectivity

BULK_CLEARANCE = 3.0    # a voxel this open, connected to the box edge, is bulk
FLOOR = 0.4             # ignore voxels tighter than this


def bulk_region(clearance):
    """Voxels that are unambiguously bulk solvent: open (>= BULK_CLEARANCE)
    and connected to the edge of the box."""
    mask = clearance >= BULK_CLEARANCE
    lab, n = ndimage.label(mask, structure=STRUCT3)
    edge = set(np.unique(lab[boundary_mask(clearance.shape)])) - {0}
    if not edge:
        return np.zeros_like(mask)
    return np.isin(lab, list(edge))


def widest_path(clearance, seed, bulk, floor=FLOOR):
    """Exact max-min (widest-path) search from `seed` to the bulk region.

    Returns (bottleneck, path_indices). Dijkstra over the max-min semiring:
    the key of a voxel is the minimum clearance along the best path reaching
    it, and voxels are finalised in decreasing order of that key, so the
    first bulk voxel popped carries the true bottleneck.
    """
    shape = clearance.shape
    ny, nz = shape[1], shape[2]

    def flat(i, j, k):
        return (i * ny + j) * nz + k

    def unflat(f):
        k = f % nz
        j = (f // nz) % ny
        i = f // (ny * nz)
        return i, j, k

    flatclear = clearance.ravel()
    flatbulk = bulk.ravel()
    sf = flat(*seed)
    if flatclear[sf] < floor:
        return None, None
    best = {sf: float(flatclear[sf])}
    prev = {}
    done = set()
    pq = [(-float(flatclear[sf]), sf)]
    strides = (ny * nz, -ny * nz, nz, -nz, 1, -1)
    while pq:
        negb, u = heapq.heappop(pq)
        if u in done:
            continue
        done.add(u)
        b = -negb
        if flatbulk[u]:
            path = [u]
            while path[-1] in prev:
                path.append(prev[path[-1]])
            path.reverse()
            return b, [unflat(f) for f in path]
        ui, uj, uk = unflat(u)
        for d, st in enumerate(strides):
            v = u + st
            # guard the array edges so we never wrap between rows/planes
            if d < 2:
                vi = ui + (1 if d == 0 else -1)
                if vi < 0 or vi >= shape[0]:
                    continue
            elif d < 4:
                vj = uj + (1 if d == 2 else -1)
                if vj < 0 or vj >= ny:
                    continue
            else:
                vk = uk + (1 if d == 4 else -1)
                if vk < 0 or vk >= nz:
                    continue
            if v in done:
                continue
            cv = float(flatclear[v])
            if cv < floor:
                continue
            nb = b if b < cv else cv
            if nb > best.get(v, -1e9):
                best[v] = nb
                prev[v] = u
                heapq.heappush(pq, (-nb, v))
    return None, None


def smooth_path(pts, n=3):
    p = np.array(pts, dtype=float)
    if len(p) < 5:
        return p
    for _ in range(n):
        q = p.copy()
        q[1:-1] = (p[:-2] + p[1:-1] + p[2:]) / 3.0
        p = q
    return p



# ------------------------------------------------- continuous refinement

class ExactClearance:
    """Continuous clearance field min_a(|p - x_a| - r_a), batched."""

    def __init__(self, atoms):
        self.X = coords(atoms)
        self.r = np.array([vdw(a.element) for a in atoms])
        self.tree = cKDTree(self.X)

    def __call__(self, pts, k=16):
        pts = np.atleast_2d(np.asarray(pts, dtype=float))
        d, idx = self.tree.query(pts, k=k, workers=-1)
        return (d - self.r[idx]).min(axis=1)


def _basis(t):
    t = t / (np.linalg.norm(t) + 1e-12)
    a = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(a, t)) > 0.9:
        a = np.array([0.0, 1.0, 0.0])
    u = np.cross(t, a)
    u /= np.linalg.norm(u) + 1e-12
    v = np.cross(t, u)
    return u, v


def refine_path(pts, clear_fn, iters=4, span=0.8, n=7):
    """Centre each path point on the local medial axis, moving only
    perpendicular to the local tangent, with Laplacian smoothing between
    sweeps so the path stays short instead of wandering into side cavities.

    Voxel-centre sampling under-reports narrow gaps and does so unevenly,
    which can move the apparent constriction onto the wrong residue;
    refinement removes that bias.
    """
    P = np.array(pts, dtype=float)
    if len(P) < 5:
        return P
    anchor = P.copy()
    for it in range(iters):
        s_ = span * (0.6 ** it)
        offs = np.linspace(-s_, s_, n)
        grid2 = np.array([(a, b) for a in offs for b in offs])
        idxs = list(range(1, len(P) - 1))
        cands = np.empty((len(idxs), len(grid2), 3))
        for j, i in enumerate(idxs):
            u, v = _basis(P[i + 1] - P[i - 1])
            cands[j] = P[i] + np.outer(grid2[:, 0], u) + \
                np.outer(grid2[:, 1], v)
        vals = clear_fn(cands.reshape(-1, 3)).reshape(len(idxs), len(grid2))
        # keep each point within 2.0 A of where the grid path put it
        drift = np.linalg.norm(cands - anchor[idxs][:, None, :], axis=2)
        vals = np.where(drift > 2.0, -1e9, vals)
        best = np.argmax(vals, axis=1)
        for j, i in enumerate(idxs):
            P[i] = cands[j, best[j]]
        P = smooth_path(P, n=1)
    return P


def densify(P, spacing=0.15):
    out = [P[0]]
    for a, b in zip(P[:-1], P[1:]):
        d = np.linalg.norm(b - a)
        m = max(1, int(np.ceil(d / spacing)))
        for t in range(1, m + 1):
            out.append(a + (b - a) * (t / m))
    return np.array(out)


# --------------------------------------------------------------- annotation

def lining_residues(atoms, point, radius, extra=1.5):
    """Residues with an atom whose vdW surface is within `extra` of a sphere
    of `radius` at `point`. Returns [(chain,resseq,resname,clearance)].

    `atoms` must be the same set used as obstructions, so that a ligand
    stripped from the calculation is not reported as lining the tunnel it
    was removed from.
    """
    out = {}
    for a in atoms:
        if a.is_hydrogen:
            continue
        d = float(np.linalg.norm(np.array(a.xyz) - point))
        clr = d - vdw(a.element)
        if clr <= radius + extra:
            key = (a.chain, a.resseq, a.resname)
            if key not in out or clr < out[key]:
                out[key] = clr
    return sorted(((c, r, n, v) for (c, r, n), v in out.items()),
                  key=lambda x: x[3])


def path_lining(atoms, pts, radii, extra=1.0):
    """Residues lining anywhere along the path (same atom set as used as
    obstructions)."""
    tree = cKDTree(pts)
    out = {}
    for a in atoms:
        if a.is_hydrogen:
            continue
        d, i = tree.query(np.array(a.xyz))
        clr = d - vdw(a.element)
        if clr <= float(radii[i]) + extra:
            key = (a.chain, a.resseq, a.resname)
            if key not in out or clr < out[key]:
                out[key] = clr
    return sorted(((c, r, n, v) for (c, r, n), v in out.items()),
                  key=lambda x: x[3])


def trimer_axis(s):
    """Pseudo-3-fold axis: normal to the plane of the three chain centroids,
    oriented periplasm-positive (towards the docking domains)."""
    cens = []
    for ch in s.chains:
        ca = s.ca(ch)
        cens.append(np.array([v for v in ca.values()]).mean(axis=0))
    cens = np.array(cens)
    c = cens.mean(axis=0)
    u, sv, vt = np.linalg.svd(cens - c)
    axis = vt[2]
    # orient towards the docking domain (periplasm)
    dock = []
    for ch in s.chains:
        ca = s.ca(ch)
        dock.append(centroid(ca, list(range(181, 278))))
    dock = np.array(dock).mean(axis=0)
    if np.dot(dock - c, axis) < 0:
        axis = -axis
    return c, axis


def channel_call(s, path_pts, lining, axis_c, axis, seed_h):
    """Tentative CH1/CH2/CH3 assignment from exit height and lining."""
    exit_pt = path_pts[-1]
    h_exit = float(np.dot(exit_pt - axis_c, axis))
    # composition of the outer half of the lining
    groups = {"PC": 0, "PN": 0, "TM": 0, "other": 0}
    for (c, r, n, v) in lining:
        if r in SUBDOMAINS["PC1"] or r in SUBDOMAINS["PC2"]:
            groups["PC"] += 1
        elif r in SUBDOMAINS["PN1"] or r in SUBDOMAINS["PN2"]:
            groups["PN"] += 1
        elif 337 <= r <= 565 or 876 <= r <= 1030 or r <= 35:
            groups["TM"] += 1
        else:
            groups["other"] += 1
    tot = max(1, sum(groups.values()))
    fPC, fPN, fTM = (groups["PC"] / tot, groups["PN"] / tot,
                     groups["TM"] / tot)
    conf = "tentative"
    if fTM > 0.45 and h_exit < seed_h:
        call = "CH2 (membrane / central-cavity facing)"
    elif fPC > fPN and fPC > 0.2:
        call = "CH1 (PC1/PC2 periplasmic cleft)"
    elif fPN > fPC and fPN > 0.2:
        call = "CH3 (PN1/PN2 groove)"
    else:
        call = "unassigned"
        conf = "uncertain - no dominant lining group"
    return call, conf, h_exit, fPC, fPN, fTM


def write_trace(path, pts, radii, chain="T", resname="TUN"):
    with open(path, "w") as fh:
        fh.write("REMARK  tunnel trace; B-factor column = local radius (A)\n")
        for i, (p, r) in enumerate(zip(pts, radii), start=1):
            fh.write(f"HETATM{i:5d}  O   {resname} {chain}{i:4d}    "
                     f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00{r:6.2f}"
                     f"           O\n")
        fh.write("END\n")


# -------------------------------------------------------------------- main

def analyse(s, mode, step, max_tunnels=3, verbose=True):
    """mode: 'protein' (ligands stripped) or 'withlig' (ligands obstruct)."""
    prot = [a for a in s.protein_atoms if not a.is_hydrogen]
    lig = [a for a in s.het_atoms if not a.is_hydrogen]
    atoms = prot + (lig if mode == "withlig" else [])
    print(f"\n  {s.name} [{mode}] {len(atoms)} obstructing heavy atoms")
    grid = ClearanceGrid(atoms, step=step)
    clear_fn = ExactClearance(atoms)
    axis_c, axis = trimer_axis(s)
    bulk = bulk_region(grid.clearance)
    print(f"    bulk region: {bulk.sum()/1e6:.1f} M voxels "
          f"(clearance >= {BULK_CLEARANCE} A, connected to the box edge)")

    # Seeds: one per chain. Ligand centroid where that chain has a ligand,
    # otherwise the substrate-site position transferred from the pocket
    # residues of the same chain (DBP/PBP midpoint).
    seeds = {}
    ligs_by_chain = {}
    for (lch, lres, lname, ats) in s.ligands():
        heavy = [a for a in ats if not a.is_hydrogen]
        ligs_by_chain.setdefault(lch, []).append(
            (lname, lres, coords(heavy).mean(axis=0)))
    for ch in s.chains:
        if ch in ligs_by_chain:
            for (lname, lres, cen) in ligs_by_chain[ch]:
                seeds[(ch, f"{lname}{lres}")] = (cen, "ligand centroid")
        else:
            ca = s.ca(ch)
            cen = 0.5 * (centroid(ca, DBP) + centroid(ca, PBP))
            seeds[(ch, "site")] = (cen, "transferred DBP/PBP midpoint")

    rows = []
    for (ch, tag), (cen, origin) in sorted(seeds.items()):
        sidx, sclr = grid.free_seed(cen)
        if sidx is None:
            print(f"    {ch}/{tag}: no free voxel near the seed - skipped")
            rows.append([s.name, mode, ch, tag, 0, origin, "", "", "", "",
                         "", "", "", "", "", "", "no free seed voxel"])
            continue
        clearance = grid.clearance.copy()
        seed_h = float(np.dot(grid.point_of(sidx) - axis_c, axis))
        for rank in range(1, max_tunnels + 1):
            R_grid, path = widest_path(clearance, sidx, bulk)
            if R_grid is None:
                if rank == 1:
                    print(f"    {ch}/{tag}: no path to bulk at "
                          f">= {FLOOR} A")
                    rows.append([s.name, mode, ch, tag, 0, origin,
                                 fmt(sclr), "", "", "", "", "", "", "", "",
                                 "",
                                 f"no path to bulk at >= {FLOOR} A"])
                break
            raw = np.array([grid.point_of(q) for q in path])
            ref = refine_path(smooth_path(raw, n=1), clear_fn)
            pts = densify(ref, spacing=0.15)
            radii = clear_fn(pts)
            length = float(np.linalg.norm(np.diff(pts, axis=0),
                                          axis=1).sum())
            j = int(np.argmin(radii))
            narrow, R = pts[j], float(radii[j])
            lin = lining_residues(atoms, narrow, R)
            own = [x for x in lin if x[0] == ch]
            top = (own or lin)[:6]
            call, conf, h_exit, fPC, fPN, fTM = channel_call(
                s, pts, lin, axis_c, axis, seed_h)
            cons = ";".join(f"{n}{r}{c}:{v:.2f}" for c, r, n, v in top)
            pathlin = path_lining(atoms, pts, radii)
            chains_seen = sorted({c for (c, r, n, v) in pathlin})
            print(f"    {ch}/{tag} #{rank}: bottleneck {R:.2f} A "
                  f"(grid {R_grid:.2f}), path {length:.1f} A, exit "
                  f"({raw[-1][0]:.0f},{raw[-1][1]:.0f},{raw[-1][2]:.0f})")
            print(f"       constriction: {cons}")
            print(f"       -> {call} [{conf}]; path lined by "
                  f"chain(s) {'/'.join(chains_seen)}")
            stp = max(1, len(pts) // 400)
            write_trace(os.path.join(
                CXDIR, f"{s.name}_{mode}_{ch}_{tag}_t{rank}_tunnel.pdb"),
                pts[::stp], radii[::stp])
            rows.append([s.name, mode, ch, tag, rank, origin, fmt(sclr),
                         fmt(R), fmt(R_grid), fmt(length, 1),
                         f"{narrow[0]:.2f},{narrow[1]:.2f},{narrow[2]:.2f}",
                         cons, call, conf, fmt(h_exit, 1),
                         "/".join(chains_seen), ""])
            # Block the constriction (not the exit): routes that share a
            # gate are the same tunnel with a different tail, so closing the
            # gate is what forces the search onto a genuinely different
            # passage.
            gi = grid.index_of(narrow)
            rad = int(np.ceil(3.5 / grid.step))
            sl = tuple(slice(max(0, gi[i] - rad),
                             min(clearance.shape[i], gi[i] + rad + 1))
                       for i in range(3))
            clearance[sl] = -1.0
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step", type=float, default=STEP)
    ap.add_argument("--structure", default=None)
    a = ap.parse_args()
    structs = load_structures()
    if a.structure:
        structs = [s for s in structs if a.structure in s.name]
    rows = []
    print("=== Stage 7: tunnels (trimer) ===")
    for s in structs:
        has_det = any(n in DETERGENTS for (_, _, n, _) in s.ligands())
        modes = ["protein"]
        if s.ligands():
            modes.append("withlig")
        for mode in modes:
            rows += analyse(s, mode, a.step)
        if has_det:
            print(f"  ({s.name} carries detergent; the protein/withlig "
                  f"difference quantifies occlusion of the exit route)")
    write_csv(os.path.join(TABLES, "tunnels.csv"),
              ["structure", "mode", "chain", "seed", "tunnel_rank",
               "seed_origin", "seed_clearance_A", "bottleneck_radius_A",
               "bottleneck_radius_grid_only_A", "geodesic_path_length_A",
               "narrowest_point_xyz", "constriction_lining_clearance_A",
               "channel_call", "assignment_confidence",
               "exit_height_along_axis_A", "chains_lining_path",
               "note"], rows)
    print(f"\nwrote results/tables/tunnels.csv ({len(rows)} rows)")


if __name__ == "__main__":
    main()
