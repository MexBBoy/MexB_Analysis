#!/usr/bin/env python3
"""Where each pocket-lining residue sits along the reference entry channel.

Gives the poster panels a residue scale to hang off, and answers a question
worth answering explicitly: does depth along the channel separate the
proximal pocket from the distal one? It does not. Both sets have residues at
26-35 A and again at 62-63 A, because the channel winds and arc length from
the mouth is not a coordinate that distinguishes two pockets sitting either
side of it. Only the residues whose side chains lie close to the centreline
are meaningful to place on this axis at all, so an offset column is written
and the panels filter on it.

Writes results/tables/channel_residues.csv.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from published_pockets import channel_depth, load_channel
from mexb_common import (DBP, PBP, STRUCT_DIR, SWITCH_LOOP, TABLES, Structure,
                         fmt, write_csv)

BACKBONE = {"N", "CA", "C", "O", "OXT"}


def main():
    chan = load_channel()
    if chan is None:
        print("  entry channel trace missing - run tunnels.py first")
        return
    s = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    sets = {"DBP": set(DBP), "PBP": set(PBP), "switch": set(SWITCH_LOOP)}

    rows = []
    for rid in sorted(set(DBP) | set(PBP) | set(SWITCH_LOOP)):
        at = s.residue_atoms("E", rid)
        if not at:
            continue
        sc = [a for a in at if a.name.strip() not in BACKBONE] or at
        cen = np.mean([a.xyz for a in sc], axis=0)
        depth, off = channel_depth(cen, chan)
        if depth is None:
            continue
        where = ";".join(k for k, v in sets.items() if rid in v)
        rows.append([rid, at[0].resname, where, fmt(depth), fmt(off)])

    write_csv(os.path.join(TABLES, "channel_residues.csv"),
              ["resseq", "resname", "site", "depth_from_entrance_A",
               "offset_from_channel_A"], rows)

    print("=== pocket-lining residues on the entry channel ===\n")
    near = [r for r in rows if float(r[4]) <= 8.0]
    print(f"  {len(near)} of {len(rows)} have a side chain within 8 A of the "
          f"centreline:\n")
    for r in sorted(near, key=lambda r: float(r[3])):
        print(f"    {r[1]}{r[0]:<5} depth {float(r[3]):5.1f} A  "
              f"offset {float(r[4]):4.1f} A  {r[2]}")

    print("\n  --- does depth separate the two pockets? ---")
    for lab in ("DBP", "PBP"):
        v = sorted(float(r[3]) for r in rows if lab in r[2])
        print(f"    {lab}: {len(v)} residues spanning {min(v):.0f}-{max(v):.0f} A"
              f", median {np.median(v):.0f}")
    print("    The ranges overlap almost completely. Depth measures distance "
          "from the\n    periplasmic mouth along a winding path; it is not a "
          "pocket coordinate.\n    Pocket identity has to come from contacts, "
          "as it does in the panels.")
    print("\nwrote results/tables/channel_residues.csv")


if __name__ == "__main__":
    main()
