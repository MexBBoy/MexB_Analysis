#!/usr/bin/env python3
"""Carry each cleft-to-ligand route on to the far end of its own pocket.

A route traced from the cleft stops at the ligand, because that is where the
search was seeded. The pocket carries on past it, so this asks how much
further, and writes the answer as an extension of the same trace.

"The far end of the pocket" is defined by burial rather than by direction: a
flood fill from bulk solvent through water-sized free space gives every voxel
a path distance from the outside, and the target is the most buried voxel of
the cavity the ligand sits in, within REACH of it. That is a measurement, not
a guess at where the channel would go if it continued - and it is worth
knowing that the answer is usually a side recess rather than a continuation,
which is why these extensions are short.

The answer comes out of the fill itself rather than out of another search:
the extension is how much more buried the cavity's deepest point is than the
ligand. A max-min search to a single buried voxel explores the whole grid and
takes tens of minutes per protomer, and what it returns is a route into a
side recess rather than a continuation of the channel.

Writes results/tables/pocket_end.csv
"""
from __future__ import annotations

import os
import sys
import time
import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tunnels as T
import per_structure_tunnels as pst
from published_pockets import PDBDIR
from mexb_common import CXDIR, STRUCT_DIR, TABLES, Structure, coords, fmt, \
    write_csv

FLOOR = 1.4         # water-sized: what counts as free space for the fill
REACH = 18.0        # how far from the ligand to look for the pocket's end


def depth_from_bulk(grid, free, bulk):
    """Path distance, in voxel steps, from bulk solvent to every free voxel.

    Grown one shell at a time with a binary dilation rather than walked with
    a queue: same answer, but the queue version touches ten million voxels
    one at a time in Python and takes tens of minutes per protomer.
    """
    seen = free & bulk
    dist = np.full(grid.shape, -1, np.int32)
    dist[seen] = 0
    cross = ndimage.generate_binary_structure(3, 1)
    for step in range(1, 4000):
        grown = ndimage.binary_dilation(seen, cross) & free & ~seen
        if not grown.any():
            break
        dist[grown] = step
        seen |= grown
    return dist


def main():
    rows_in = pst.rows_of(os.path.join(TABLES, "cleft_routes.csv"))
    if not rows_in:
        print("  run cleft_to_ligand.py first")
        return
    print("=== from the ligand on to the end of its pocket ===")
    out = []
    for r in rows_in:
        t0 = time.time()
        pid, ch, nm = r["pdb"], r["chain"], r["ligand"]
        trace = os.path.join(CXDIR, r["trace_file"])
        if not os.path.exists(trace):
            continue
        P, rad = pst.read_trace(trace)      # runs ligand -> cleft

        path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            path = os.path.join(PDBDIR, f"{pid}.pdb")
        s = Structure(path)
        cen = {c: coords([a for a in s.protein_atoms if a.chain == c
                          and a.name.strip() == "CA"]).mean(0)
               for c in s.chains}
        near = sorted(cen, key=lambda c: float(
            np.linalg.norm(cen[c] - cen[ch])))[:3]
        atoms = [a for a in s.protein_atoms
                 if not a.is_hydrogen and a.chain in set(near)]
        grid = T.ClearanceGrid(atoms, step=T.STEP, verbose=False)
        clear_fn = T.ExactClearance(atoms)
        bulk = T.bulk_region(grid.clearance)
        free = grid.clearance >= FLOOR
        dist = depth_from_bulk(grid, free, bulk)

        src = grid.free_seed(P[0])[0]       # the ligand end of the route
        if src is None or dist[src] < 0:
            print(f"  {nm}: no free voxel at the ligand"); continue
        rad_v = int(np.ceil(REACH / grid.step))
        sl = tuple(slice(max(0, src[i] - rad_v),
                         min(grid.shape[i], src[i] + rad_v + 1))
                   for i in range(3))
        sub = dist[sl].copy()
        gx = [np.arange(sl[i].start, sl[i].stop) for i in range(3)]
        Q = np.stack(np.meshgrid(*gx, indexing="ij"), -1) * grid.step \
            + grid.origin
        sub[np.linalg.norm(Q - P[0], axis=-1) > REACH] = -1

        # How much further the cavity goes, as a length rather than a path.
        # The flood fill already gives every voxel its distance from bulk, so
        # the extension is how much more buried the deepest point of this
        # cavity is than the ligand. Tracing a route to that point is both
        # expensive - a max-min search to one buried voxel explores the whole
        # grid, and did - and beside the point, since the point is a side
        # recess rather than a continuation of the path.
        deepest = int(sub.max())
        ext = max(0.0, (deepest - int(dist[src])) * grid.step)
        k = np.unravel_index(int(np.argmax(sub)), sub.shape)
        tgt = tuple(sl[i].start + k[i] for i in range(3))
        side = float(np.linalg.norm(grid.point_of(tgt) - P[0]))
        out.append([pid, ch, nm, fmt(ext), fmt(float(grid.clearance[tgt])),
                    fmt(float(deepest * grid.step)),
                    fmt(float(int(dist[src]) * grid.step)), fmt(side),
                    r["trace_file"]])
        print(f"  {nm:16} {pid} {ch}: cavity carries on {ext:4.1f} A past "
              f"the ligand, to a point {side:4.1f} A away with "
              f"{grid.clearance[tgt]:.2f} A of room  "
              f"({time.time() - t0:.0f}s)")

    write_csv(os.path.join(TABLES, "pocket_end.csv"),
              ["pdb", "chain", "ligand", "extension_A", "pocket_end_room_A",
               "pocket_end_burial_A", "ligand_burial_A",
               "straight_line_to_end_A", "trace_file"], out)
    print("\nwrote results/tables/pocket_end.csv")


if __name__ == "__main__":
    main()
