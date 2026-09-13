#!/usr/bin/env python3
"""Rebuild ligand_routes.csv from the traces ligand_routes.py has written.

The search is the expensive part and each route is written to disk as soon as
it is found, but the table was only assembled at the end of the run, so an
interrupted run threw away every completed search. Nothing in the table needs
the grid: the trace carries its own radius per point, so the strain integral
is just max(0, r_lig - radius) integrated along the trace, and the ligand
supplies the rest. This recovers the table from whatever traces exist.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import per_structure_tunnels as pst
from ligand_routes import girth, sphere_radius
from published_pockets import PDBDIR, pocket_ligands
from mexb_common import CXDIR, STRUCT_DIR, TABLES, Structure, coords, fmt, write_csv


def main():
    rows = []
    for (pid, ch), nm in sorted(pst.panel_protomers().items(), key=lambda x: x[1]):
        path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            path = os.path.join(PDBDIR, f"{pid}.pdb")
        if not os.path.exists(path):
            continue
        s = Structure(path)
        seen = {}
        for (c, rn, heavy) in pocket_ligands(s):
            if c != ch:
                continue
            f = os.path.join(CXDIR, f"lig_{pid}_{ch}_{rn}_tunnel.pdb")
            if not os.path.exists(f):
                continue
            # several copies of one resname share a trace name; take each once
            if rn in seen:
                continue
            seen[rn] = True
            P, rad = pst.read_trace(f)
            arc = np.concatenate([[0.0], np.cumsum(
                np.linalg.norm(np.diff(P, axis=0), axis=1))])
            ds = np.diff(arc)
            rloc, rwhole, length = girth(heavy)
            need = np.maximum(0.0, rloc - rad)
            strain = float((need[1:] * ds).sum())
            peak = float(need.max())
            frac = float((need > 0.05).mean())
            total = float(arc[-1])
            straight = float(np.linalg.norm(P[-1] - P[0]))
            rows.append([pid, ch, nm, rn, len(heavy), fmt(rloc), fmt(rwhole),
                         fmt(sphere_radius(heavy)), fmt(length), "",
                         fmt(total), fmt(straight),
                         fmt(total / straight if straight > 0 else np.nan),
                         fmt(strain), fmt(peak), fmt(100.0 * frac),
                         fmt(float(rad.min())), os.path.basename(f)])
            print(f"  {nm:16} {rn:>4} route {total:5.1f} A  strain {strain:6.1f}"
                  f" A^2  peak {peak:.2f} A  {100*frac:3.0f}% pinched")
    write_csv(os.path.join(TABLES, "ligand_routes.csv"),
              ["pdb", "chain", "ligand", "resname", "heavy_atoms",
               "local_girth_A", "whole_girth_A", "sphere_radius_A",
               "long_axis_A", "widest_rigid_A", "route_length_A",
               "straight_line_A", "tortuosity", "strain_integral_A2",
               "peak_opening_A", "percent_pinched", "route_min_radius_A",
               "trace_file"], rows)
    print(f"\nwrote results/tables/ligand_routes.csv ({len(rows)} routes)")


if __name__ == "__main__":
    main()
