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
from scipy.spatial import cKDTree

from mexb_common import (DBP, PBP, STRUCT_DIR, SUBDOMAINS, SWITCH_LOOP,
                         TABLES, Structure, coords, fmt, write_csv)

BACKBONE = {"N", "CA", "C", "O", "OXT"}
KEYS = ["PN1", "PN2", "PC1", "PC2"]
NEAR = 9.0     # a channel point is "lined by" atoms within this radius
SMOOTH = 25    # trace points, about 5 A of path


def write_subdomains(s, chan):
    """Which porter subdomain lines the channel, as a function of depth.

    For every point on the trace, the porter-subdomain membership of the
    protein atoms around it, smoothed along the path. The raw per-point
    assignment flips between neighbours; smoothed it resolves into a clean
    gradient from the PC1/PC2 cleft at the mouth to PN2/PC1 at the deep end,
    which is the anatomy the channel is usually described by.
    """
    import collections
    P, arc, total = chan
    at = [a for a in s.protein_atoms if a.chain == "E" and not a.is_hydrogen]
    tree = cKDTree(coords(at))
    sub = {}
    for k, v in SUBDOMAINS.items():
        for r in v:
            sub.setdefault(r, []).append(k)

    F = np.zeros((len(P), len(KEYS)))
    for i, pt in enumerate(P):
        c = collections.Counter()
        for j in tree.query_ball_point(pt, NEAR):
            for k in sub.get(at[j].resseq, []):
                c[k] += 1
        tot = sum(c.values()) or 1
        F[i] = [c[k] / tot for k in KEYS]

    depth = total - arc
    o = np.argsort(depth)
    d, F = depth[o], F[o]
    ker = np.ones(SMOOTH) / SMOOTH
    Sm = np.vstack([np.convolve(F[:, j], ker, mode="same")
                    for j in range(len(KEYS))]).T
    Sm = Sm / np.maximum(Sm.sum(1, keepdims=True), 1e-9)

    rows = [[fmt(d[i])] + [fmt(Sm[i, j], 3) for j in range(len(KEYS))]
            for i in range(len(d))]
    write_csv(os.path.join(TABLES, "channel_subdomains.csv"),
              ["depth_from_entrance_A"] + [f"fraction_{k}" for k in KEYS],
              rows)

    print("\n  --- which subdomain lines the channel, by depth ---")
    for lo in range(0, 65, 10):
        m = (d >= lo) & (d < lo + 10)
        if m.sum() < 3:
            continue
        v = Sm[m].mean(0)
        rank = sorted(zip(KEYS, v), key=lambda x: -x[1])
        print(f"    {lo:2d}-{lo + 10:2d} A: "
              + "  ".join(f"{k} {100 * x:3.0f}%" for k, x in rank[:3]))
    print("  wrote results/tables/channel_subdomains.csv")


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
    write_subdomains(s, chan)
    print("\nwrote results/tables/channel_residues.csv")


if __name__ == "__main__":
    main()
