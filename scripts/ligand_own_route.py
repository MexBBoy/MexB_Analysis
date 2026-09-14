#!/usr/bin/env python3
"""Chloramphenicol's own route out, and what it says about its site.

P21 traces one route per protomer, in at the CH1 periplasmic cleft, and every
ligand sits on it except chloramphenicol, which stays 13 A off however the
porter pockets are threaded. Its site is lined 64 percent by PN1/PN2 and only
5 percent by PC1/PC2, which looked like it bound at the CH3 groove - a
different entrance that no rerouting through the porter pockets could reach.

This was written to trace that entrance, and it refuted the idea. The route is
seeded at the ligand and allowed to leave by whatever opening is widest, with
the mouth it finds then named by the subdomains lining its last 12 A, the same
test P18 uses - so the exit is measured, not assumed. The measured mouth is
CH1: 86 percent PC1/PC2 and 0 percent PN1/PN2, the same periplasmic cleft
every other row already uses.

So the PN-rich lining describes the chamber chloramphenicol sits in, not the
way it got there. It occupies a side branch off the CH1 route rather than a
separate entrance, which is why a line drawn through the porter pockets misses
it by 13 A while a line seeded at the ligand passes within 1.4 A. The file is
named for the route rather than for the channel it was expected to find.

From the ligand the route carries on the way a substrate travels, through the
distal pocket and out at the funnel on the trimer's three-fold axis, so the
result is comparable end to end with the CH1 rows.

Writes results/tables/ch3_tunnels.csv, results/tables/ch3_tunnel_ligands.csv
and a trace per protomer under results/chimerax/ch3_<pdb>_<chain>.pdb
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tunnels as T
from all_channels import exit_call
from cleft_to_ligand import free_point
from full_tunnels import FUNNEL_R, ABOVE, ball, leg, trimer_axis
from published_pockets import PDBDIR, pocket_ligands
from mexb_common import (CXDIR, DBP, STRUCT_DIR, TABLES, Structure, centroid,
                         coords, fmt, write_csv)

WANT = {("21FP", "B"): "Chloramphenicol"}


def main():
    print("=== chloramphenicol's own route out, exit measured not assumed ===")
    rows, ligrows = [], []
    for (pid, ch), nm in sorted(WANT.items(), key=lambda x: x[1]):
        t0 = time.time()
        path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            path = os.path.join(PDBDIR, f"{pid}.pdb")
        if not os.path.exists(path):
            print(f"  {nm}: no structure"); continue
        s = Structure(path)

        cen = {c: coords([a for a in s.protein_atoms if a.chain == c
                          and a.name.strip() == "CA"]).mean(0) for c in s.chains}
        near = sorted(cen, key=lambda c: float(
            np.linalg.norm(cen[c] - cen[ch])))[:3]
        atoms = [a for a in s.protein_atoms
                 if not a.is_hydrogen and a.chain in set(near)]
        grid = T.ClearanceGrid(atoms, step=T.STEP, verbose=False)
        clear_fn = T.ExactClearance(atoms)
        bulk = T.bulk_region(grid.clearance)
        axcen, ax, allca = trimer_axis(s)
        top = float(((allca - axcen) @ ax).max())

        pl = [(rn, h) for (c, rn, h) in pocket_ligands(s) if c == ch]
        if not pl:
            print(f"  {nm}: no pocket ligand"); continue
        rn, heavy = pl[0]
        seed, _ = free_point(clear_fn, coords(heavy).mean(0))
        if seed is None:
            print(f"  {nm}: no free voxel at the ligand"); continue
        sidx = grid.free_seed(seed)[0]

        # out of the ligand by whatever opening is widest: the entrance is
        # found, not assumed, and then named by what lines it
        A, ra = leg(grid, clear_fn, sidx, bulk)
        if A is None:
            print(f"  {nm}: the ligand reaches no opening at all"); continue
        pa = T.densify(T.refine_path(T.smooth_path(A, n=1), clear_fn),
                       spacing=0.15)
        name, fpc, fpn, ftm, fdk = exit_call(atoms, pa, clear_fn(pa))

        # on from the ligand the way a substrate travels: distal, then funnel
        exi = None
        for up in ABOVE:
            fp = axcen + ax * (top + up)
            for r in FUNNEL_R:
                m = ball(grid, bulk, fp, r)
                if m.any():
                    exi = m
                    break
            if exi is not None:
                break
        ca = s.ca(ch)
        q, _ = free_point(clear_fn, centroid(ca, DBP))
        M = B = None
        if q is not None and exi is not None:
            sd = grid.free_seed(q)[0]
            tgt = np.zeros(grid.shape, bool); tgt[sd] = True
            M, rm = leg(grid, clear_fn, sidx, tgt)
            B, rb = leg(grid, clear_fn, sd, exi)
        if M is None or B is None:
            print(f"  {nm}: entrance found ({name}) but no way on to the funnel")
            continue

        pts = np.vstack([A[::-1], M[1:], B[1:]])
        pts = T.densify(T.refine_path(T.smooth_path(pts, n=1), clear_fn),
                        spacing=0.15)
        rad = clear_fn(pts)
        arc = np.concatenate([[0.0], np.cumsum(
            np.linalg.norm(np.diff(pts, axis=0), axis=1))])
        total = float(arc[-1])
        out = os.path.join(CXDIR, f"ch3_{pid}_{ch}.pdb")
        T.write_trace(out, pts, rad)
        rows.append([pid, ch, nm, fmt(total), fmt(min(ra, rm, rb)),
                     fmt(float(rad.min())), name, fmt(100 * fpc), fmt(100 * fpn),
                     fmt(100 * ftm), fmt(100 * fdk), os.path.basename(out)])
        print(f"  {nm}: entrance is {name}")
        print(f"    lining of that mouth: PC {100*fpc:.0f}%  PN {100*fpn:.0f}%"
              f"  TM {100*ftm:.0f}%  dock {100*fdk:.0f}%")
        print(f"    {total:.1f} A end to end, neck {min(ra, rm, rb):.2f} A"
              f"  ({time.time()-t0:.0f}s)")

        for (c, r2, h) in pocket_ligands(s):
            if c != ch:
                continue
            X = coords(h)
            d = np.linalg.norm(X[:, None, :] - pts[None, :, :], axis=2)
            j = d.argmin(1)
            off = d[np.arange(len(X)), j]
            ligrows.append([pid, ch, nm, r2, len(X), fmt(float(arc[j].mean())),
                            fmt(float(arc[j].min())), fmt(float(arc[j].max())),
                            fmt(float(off.mean())), fmt(float(off.min())),
                            fmt(float(np.interp(arc[j].mean(), arc, rad))),
                            fmt(total)])
            print(f"    {r2} now sits {off.min():.2f} A from this line, "
                  f"{100*arc[j].mean()/total:.0f}% along")

    write_csv(os.path.join(TABLES, "ch3_tunnels.csv"),
              ["pdb", "chain", "ligand", "length_A", "leg_bottleneck_A",
               "trace_min_radius_A", "entrance", "mouth_PC_pct", "mouth_PN_pct",
               "mouth_TM_pct", "mouth_dock_pct", "trace_file"], rows)
    write_csv(os.path.join(TABLES, "ch3_tunnel_ligands.csv"),
              ["pdb", "chain", "ligand", "resname", "heavy_atoms",
               "along_tunnel_A", "along_min_A", "along_max_A", "mean_offset_A",
               "closest_offset_A", "radius_there_A", "tunnel_length_A"],
              ligrows)
    print(f"\nwrote results/tables/ch3_tunnels.csv ({len(rows)})")


if __name__ == "__main__":
    main()
