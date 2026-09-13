#!/usr/bin/env python3
"""The widest route from the periplasmic cleft to each bound ligand.

Traces so far were seeded at the ligand and allowed to leave by whatever
opening was widest. That is the right question for "where can this protomer
let something out", but it is why the rows of P19 start at different depths:
all but ampicillin's route breaks out of the side of the porter domain rather
than at the periplasmic mouth, so it never covers the shallow part of the
path.

This asks the transport question instead - what is the widest route from the
periplasmic cleft to where the substrate actually sits - by fixing both ends:
the seed at the ligand, and the exit at the cleft mouth of the reference
channel, mapped into each structure's own frame by superposition.

Two stages, because a max-min search alone is not enough. Widest-path gives
the bottleneck R*, but among routes that all achieve R* it has no preference
and wanders: the same search once returned a 220 A excursion for a 15 A
separation. Taking the shortest route through the voxels with clearance >= R*
gives the same bottleneck by a direct path.

Writes results/tables/cleft_routes.csv and a trace per protomer under
results/chimerax/cleft_<pdb>_<chain>_tunnel.pdb
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
from published_pockets import LINING, PDBDIR, load_channel, pocket_ligands
from mexb_common import (CXDIR, DBP, PBP, STRUCT_DIR, TABLES, Structure,
                         apply_rt, centroid, coords, fmt, kabsch, write_csv)

MOUTH_R = 9.0       # bulk voxels this close to the mapped mouth are the exit
END_TRIM = 3.0
REACH = 11.0


def free_point(clear_fn, p, reach=REACH, step=0.8, floor=1.4):
    g = np.arange(-reach, reach + 1e-9, step)
    off = np.array([[x, y, z] for x in g for y in g for z in g])
    off = off[np.linalg.norm(off, axis=1) <= reach]
    r = clear_fn(p + off)
    k = int(np.argmax(r - 0.02 * np.linalg.norm(off, axis=1)))
    return (p + off[k], float(r[k])) if r[k] >= floor else (None, 0.0)


def shortest_within(grid, thr, src, dst):
    """Shortest route from src to any dst voxel, through clearance >= thr."""
    ok = grid.clearance >= thr
    off = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1)
           for c in (-1, 0, 1) if (a, b, c) != (0, 0, 0)]
    step = grid.step
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
            return out[::-1]
        if d > D.get(v, np.inf):
            continue
        for o in off:
            w = (v[0] + o[0], v[1] + o[1], v[2] + o[2])
            if not (0 <= w[0] < grid.shape[0] and 0 <= w[1] < grid.shape[1]
                    and 0 <= w[2] < grid.shape[2]) or not ok[w]:
                continue
            nd = d + step * float(np.sqrt(o[0] ** 2 + o[1] ** 2 + o[2] ** 2))
            if nd < D.get(w, np.inf):
                D[w] = nd
                prev[w] = v
                heapq.heappush(q, (nd, w))
    return None


def ligand_positions(rows):
    """Every pocket ligand placed on its protomer's cleft route.

    Cheap second pass over the traces already written: for each ligand atom,
    the nearest point of the route, so a ligand gets a span and a mean along
    the same axis the tunnel is drawn on.
    """
    out = []
    for (pid, ch, nm, *_rest) in rows:
        f = os.path.join(CXDIR, f"cleft_{pid}_{ch}_tunnel.pdb")
        if not os.path.exists(f):
            continue
        P, rad = pst.read_trace(f)
        arc = np.concatenate([[0.0], np.cumsum(
            np.linalg.norm(np.diff(P, axis=0), axis=1))])
        path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            path = os.path.join(PDBDIR, f"{pid}.pdb")
        s = Structure(path)
        for (c, rn, h) in pocket_ligands(s):
            if c != ch:
                continue
            X = coords(h)
            d = np.linalg.norm(X[:, None, :] - P[None, :, :], axis=2)
            j = d.argmin(1)
            off = d[np.arange(len(X)), j]
            # the trace runs ligand -> cleft, so distance from the cleft
            # is the arc measured back from its far end
            frm = arc[-1] - arc[j]
            out.append([pid, ch, nm, rn, len(X), fmt(float(frm.mean())),
                        fmt(float(frm.min())), fmt(float(frm.max())),
                        fmt(float(off.mean())), fmt(float(off.min())),
                        fmt(float(np.interp(arc[j].mean(), arc, rad)))])
    write_csv(os.path.join(TABLES, "cleft_ligands.csv"),
              ["pdb", "chain", "ligand", "resname", "heavy_atoms",
               "along_route_A", "along_route_min_A", "along_route_max_A",
               "mean_offset_A", "closest_offset_A", "radius_there_A"], out)
    print("  wrote results/tables/cleft_ligands.csv")


def main():
    chan = load_channel()
    RP, rarc, rtot = chan
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    rca = ref.ca("E")
    want = pst.panel_protomers()

    print("=== the widest route from the cleft to each ligand ===")
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
        # reference -> this structure, so the reference mouth can be placed
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

        # the exit: bulk voxels around the reference mouth, in this frame
        gx = [grid.origin[i] + grid.step * np.arange(grid.shape[i])
              for i in range(3)]
        d2 = ((gx[0][:, None, None] - mouth[0]) ** 2
              + (gx[1][None, :, None] - mouth[1]) ** 2
              + (gx[2][None, None, :] - mouth[2]) ** 2)
        dst = bulk & (d2 <= MOUTH_R ** 2)
        if not dst.any():
            print(f"  {nm}: no bulk voxel within {MOUTH_R:.0f} A of the mouth")
            continue

        pl = [h for (c, rn, h) in pocket_ligands(s) if c == ch]
        if not pl:
            print(f"  {nm}: no pocket ligand"); continue
        ca = s.ca(ch)
        dbp = centroid(ca, DBP)
        seed0 = min((coords(h).mean(0) for h in pl),
                    key=lambda c: float(np.linalg.norm(c - dbp)))
        seed, clr = free_point(clear_fn, seed0)
        if seed is None:
            print(f"  {nm}: no free voxel near the ligand"); continue
        sidx = grid.free_seed(seed)[0]

        Rg, _ = T.widest_path(grid.clearance, sidx, dst)
        if Rg is None:
            print(f"  {nm}: no route from the ligand to the cleft mouth")
            continue
        vox = shortest_within(grid, Rg - 1e-6, sidx, dst)
        pts = T.densify(T.refine_path(T.smooth_path(
            np.array([grid.point_of(q) for q in vox]), n=1), clear_fn),
            spacing=0.15)
        rad = clear_fn(pts)
        arc = np.concatenate([[0.0], np.cumsum(
            np.linalg.norm(np.diff(pts, axis=0), axis=1))])
        total = float(arc[-1])
        depth = total - arc                  # 0 at the cleft, deep at the end
        mid = (depth >= END_TRIM) & (depth <= total - END_TRIM)
        neck = float(rad[mid].min()) if mid.any() else float(rad.min())

        out = os.path.join(CXDIR, f"cleft_{pid}_{ch}_tunnel.pdb")
        T.write_trace(out, pts, rad)
        rows.append([pid, ch, nm, fmt(total), fmt(Rg), fmt(neck),
                     fmt(float(rad[0])), fmt(clr), os.path.basename(out)])
        print(f"  {nm:16} {pid} {ch}: {total:5.1f} A from the cleft to the "
              f"ligand, bottleneck {neck:.2f} A  ({time.time() - t0:.0f}s)")

    ligand_positions(rows)

    write_csv(os.path.join(TABLES, "cleft_routes.csv"),
              ["pdb", "chain", "ligand", "length_A", "grid_bottleneck_A",
               "route_bottleneck_A", "radius_at_ligand_A",
               "seed_clearance_A", "trace_file"], rows)
    print("\nwrote results/tables/cleft_routes.csv")


if __name__ == "__main__":
    main()
