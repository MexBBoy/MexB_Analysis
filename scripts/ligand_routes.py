#!/usr/bin/env python3
"""The route each ligand takes: the line of least deformation, per molecule.

The widest-path traces answer "where is the roomiest way into this protomer".
That is a property of the protein: ask it about a 20-atom antibiotic or a
69-atom detergent in the same structure and it returns the same line. It is
the right question for "how open is this protomer", and the wrong one for
"how did this molecule get there".

A hard size gate does not fix it. Every ligand here is 4.7 to 11.0 A across
while the channel it sits in necks to 2.1-2.9 A, so a rigid gate blocks all
of them, the probe clamps to the widest available radius and every route
collapses back onto the widest path. That is not an artefact: these are
snapshots with the ligand already bound, and a rigid protein has no opening
wide enough to have admitted it. The channel must breathe either way.

So the question worth asking is not "where does it fit" but "where does it
have to force the protein open the least". Two changes make the route a
property of the molecule:

  - Local girth. A detergent does not thread as a rigid rod; what has to pass
    a given constriction is its widest locally rigid piece - a sugar ring, an
    alkyl segment - not the cross section of the whole molecule. We slide a
    window along the ligand's long axis and take the widest local cross
    section, measured to van der Waals surfaces.

  - A deformation cost instead of a gate. Stepping through a voxel costs
    max(0, r_lig - clearance): how far the channel must open there, for this
    ligand. The route minimises the integral of that cost, with a small term
    for length so that among equally unstrained routes the direct one wins.

The cost field is ligand specific, so the routes genuinely differ. A slim
substrate pays nothing through a 3 A region and goes straight; a bulky one
pays everywhere and bends towards the wide parts even at the cost of extra
distance. The integral itself is the interesting number - the total breathing
each molecule demands of the protein to reach its own site.

Writes results/tables/ligand_routes.csv and a trace per ligand under
results/chimerax/lig_<pdb>_<chain>_<resname>_tunnel.pdb
"""
from __future__ import annotations

import heapq
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tunnels as T
import per_structure_tunnels as pst
from cleft_to_ligand import MOUTH_R, free_point
from published_pockets import LINING, PDBDIR, load_channel, pocket_ligands
from mexb_common import (CXDIR, DBP, STRUCT_DIR, TABLES, Structure, apply_rt,
                         centroid, coords, fmt, kabsch, vdw, write_csv)


PASS_FLOOR = 1.0    # below this the grid is protein, not a squeezable gap
LAMB = 40.0         # weight on deformation against path length


def _spine(X, u, k=7):
    """A curved centre line through the ligand, atoms ordered along `u`."""
    o = np.argsort(u)
    Z = X[o]
    k = max(3, min(k, len(Z)))
    pad = np.vstack([np.repeat(Z[:1], k // 2, 0), Z,
                     np.repeat(Z[-1:], k // 2, 0)])
    ker = np.ones(k) / k
    return np.column_stack([np.convolve(pad[:, d], ker, mode="valid")
                            for d in range(3)])


def _dist_to_polyline(P, C):
    """Distance from each point of P to the polyline C."""
    A, B = C[:-1], C[1:]
    AB = B - A
    L2 = np.einsum("ij,ij->i", AB, AB)
    L2[L2 == 0] = 1e-9
    d = P[:, None, :] - A[None, :, :]
    t = np.clip(np.einsum("ijk,jk->ij", d, AB) / L2, 0.0, 1.0)
    proj = A[None, :, :] + t[:, :, None] * AB[None, :, :]
    return np.linalg.norm(P[:, None, :] - proj, axis=2).min(1)


def girth(heavy):
    """The ligand's tube radius about its own curved centre line.

    A flexible molecule threads a channel along its own length, so what has to
    clear a constriction is its thickness about the line it follows, not its
    extent about a straight axis. Atoms are ordered along the first principal
    axis, smoothed into a centre line, and each atom's distance to that line
    plus its own van der Waals radius gives the tube the molecule sweeps.

    Measuring to a straight axis instead charges a curled detergent for its
    own curvature: the bend carries atoms far from the axis even where the
    chain is locally thin, so the straight figure reports a girth the molecule
    never actually has to present. Both are returned - `whole` is the straight
    measure, kept so the difference is visible.

    Returns (tube radius, straight-axis girth, long-axis length).
    """
    X = coords(heavy)
    rv = np.array([vdw(a.element) for a in heavy])
    if len(X) < 4:
        return float(rv.max()), float(rv.max()), 0.0
    Y = X - X.mean(0)
    _, _, vt = np.linalg.svd(Y, full_matrices=False)
    u = Y @ vt[0]
    whole = float((np.linalg.norm(Y @ vt[1:].T, axis=1) + rv).max())
    tube = float((_dist_to_polyline(X, _spine(X, u)) + rv).max())
    return tube, whole, float(np.ptp(u))


def sphere_radius(heavy):
    X = coords(heavy)
    d = np.linalg.norm(X - X.mean(0), axis=1)
    return float((d + np.array([vdw(a.element) for a in heavy])).max())


def least_strain(grid, bulk, rlig, src, dst):
    """Route from src to dst minimising how far the protein must open.

    Dijkstra on the cost  (1 + LAMB * max(0, rlig - clearance)) ds. The first
    term is plain path length, so among routes that never strain the protein
    the shortest wins; the second is the opening this ligand demands at each
    voxel, which dominates when it is non-zero. Voxels below PASS_FLOOR are
    protein rather than a squeezable gap and are not crossed at any price.

    Bulk solvent is barred except at the target. Outside the protein the
    clearance is large, so every bulk voxel is unstrained and costs only its
    length: left open, the cheapest route leaves by the nearest surface and
    travels around the outside of the protein to the mouth, which is not a
    channel and makes the search explore the whole solvent box. Barring bulk
    keeps the route interior until it surfaces at the cleft, where it ends.

    Returns (voxel path, integrated strain in A^2, peak strain in A).
    """
    clear = grid.clearance
    ok = (clear >= PASS_FLOOR) & (~bulk | dst)
    pen = np.maximum(0.0, rlig - clear)
    off = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1)
           for c in (-1, 0, 1) if (a, b, c) != (0, 0, 0)]
    step = grid.step
    sh = grid.shape
    D = {src: 0.0}
    prev = {}
    q = [(0.0, src)]
    while q:
        d, v = heapq.heappop(q)
        if dst[v]:
            out = [v]
            while v != src:
                v = prev[v]
                out.append(v)
            out = out[::-1]
            strain = 0.0
            peak = 0.0
            for i in range(1, len(out)):
                a, b = out[i - 1], out[i]
                ds = step * float(np.sqrt(sum((a[k] - b[k]) ** 2
                                              for k in range(3))))
                strain += ds * float(pen[b])
                peak = max(peak, float(pen[b]))
            return out, strain, peak
        if d > D.get(v, np.inf):
            continue
        for o in off:
            w = (v[0] + o[0], v[1] + o[1], v[2] + o[2])
            if not (0 <= w[0] < sh[0] and 0 <= w[1] < sh[1]
                    and 0 <= w[2] < sh[2]) or not ok[w]:
                continue
            ds = step * float(np.sqrt(o[0] ** 2 + o[1] ** 2 + o[2] ** 2))
            nd = d + ds * (1.0 + LAMB * float(pen[w]))
            if nd < D.get(w, np.inf):
                D[w] = nd
                prev[w] = v
                heapq.heappush(q, (nd, w))
    return None, np.nan, np.nan


def main():
    RP, rarc, rtot = load_channel()
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    rca = ref.ca("E")
    want = pst.panel_protomers()

    print("=== the route each ligand can take, gated on its own girth ===")
    rows = []
    for (pid, ch), nm in sorted(want.items(), key=lambda x: x[1]):
        t0 = time.time()
        path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            path = os.path.join(PDBDIR, f"{pid}.pdb")
        if not os.path.exists(path):
            continue
        s = Structure(path)
        mca = s.ca(ch)
        common = [r for r in LINING if r in mca and r in rca]
        if len(common) < 25:
            print(f"  {nm}: too few lining CA in common"); continue
        R, t = kabsch(np.array([rca[r] for r in common]),
                      np.array([mca[r] for r in common]))
        mouth = apply_rt(R, t, RP[-1][None, :])[0]

        cen = {c: coords([a for a in s.protein_atoms if a.chain == c
                          and a.name.strip() == "CA"]).mean(0) for c in s.chains}
        near = sorted(cen, key=lambda c: float(
            np.linalg.norm(cen[c] - cen[ch])))[:3]
        atoms = [a for a in s.protein_atoms
                 if not a.is_hydrogen and a.chain in set(near)]
        grid = T.ClearanceGrid(atoms, step=T.STEP, verbose=False)
        clear_fn = T.ExactClearance(atoms)
        bulk = T.bulk_region(grid.clearance)

        gx = [grid.origin[i] + grid.step * np.arange(grid.shape[i])
              for i in range(3)]
        d2 = ((gx[0][:, None, None] - mouth[0]) ** 2
              + (gx[1][None, :, None] - mouth[1]) ** 2
              + (gx[2][None, None, :] - mouth[2]) ** 2)
        dst = bulk & (d2 <= MOUTH_R ** 2)
        if not dst.any():
            print(f"  {nm}: no bulk voxel near the mouth"); continue

        ca = s.ca(ch)
        dbp = centroid(ca, DBP)
        pl = [(rn, h) for (c, rn, h) in pocket_ligands(s) if c == ch]
        if not pl:
            print(f"  {nm}: no pocket ligand"); continue

        for (rn, heavy) in sorted(pl, key=lambda x: float(np.linalg.norm(
                coords(x[1]).mean(0) - dbp))):
            rloc, rwhole, length = girth(heavy)
            rsph = sphere_radius(heavy)
            seed, clr = free_point(clear_fn, coords(heavy).mean(0))
            if seed is None:
                print(f"  {nm} {rn}: no free voxel near the ligand"); continue
            sidx = grid.free_seed(seed)[0]

            # how wide a rigid body could get there at all, for reference
            Rmax, _ = T.widest_path(grid.clearance, sidx, dst)
            vox, strain, peak = least_strain(grid, bulk, rloc, sidx, dst)
            if vox is None:
                print(f"  {nm} {rn}: no route to the cleft mouth"); continue
            pts = T.densify(T.refine_path(T.smooth_path(
                np.array([grid.point_of(q) for q in vox]), n=1), clear_fn),
                spacing=0.15)
            rad = clear_fn(pts)
            arc = np.concatenate([[0.0], np.cumsum(
                np.linalg.norm(np.diff(pts, axis=0), axis=1))])
            total = float(arc[-1])
            straight = float(np.linalg.norm(pts[-1] - pts[0]))
            # the stretch that actually has to open, on the refined trace
            need = np.maximum(0.0, rloc - rad)
            frac = float((need > 0.05).mean())

            out = os.path.join(CXDIR, f"lig_{pid}_{ch}_{rn}_tunnel.pdb")
            T.write_trace(out, pts, rad)
            rows.append([pid, ch, nm, rn, len(heavy), fmt(rloc), fmt(rwhole),
                         fmt(rsph), fmt(length),
                         fmt(Rmax if Rmax is not None else np.nan),
                         fmt(total), fmt(straight),
                         fmt(total / straight if straight > 0 else np.nan),
                         fmt(strain), fmt(peak), fmt(100.0 * frac),
                         fmt(float(rad.min())), os.path.basename(out)])
            print(f"  {nm:16} {rn:>4} {len(heavy):3d} atoms  girth {rloc:.2f} A"
                  f"  route {total:5.1f} A ({total/straight:.2f}x direct)"
                  f"  strain {strain:6.1f} A^2  peak {peak:.2f} A"
                  f"  {100*frac:3.0f}% pinched")
        print(f"    [{nm} {time.time()-t0:.0f}s]")

    write_csv(os.path.join(TABLES, "ligand_routes.csv"),
              ["pdb", "chain", "ligand", "resname", "heavy_atoms",
               "local_girth_A", "whole_girth_A", "sphere_radius_A",
               "long_axis_A", "widest_rigid_A", "route_length_A",
               "straight_line_A", "tortuosity", "strain_integral_A2",
               "peak_opening_A", "percent_pinched", "route_min_radius_A",
               "trace_file"], rows)
    print("\nwrote results/tables/ligand_routes.csv")


if __name__ == "__main__":
    main()
