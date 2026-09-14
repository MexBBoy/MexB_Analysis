#!/usr/bin/env python3
"""Where the two binding pockets fall along each traced tunnel.

P21 and P22 draw the path and the substrates on it, but not the landmarks the
whole mechanism is described in terms of: the proximal pocket, where a
substrate first binds, and the distal pocket across the switch loop from it.
Naming them in the caption is not the same as showing where they are.

A pocket is a chamber, not a point, and the traced line skirts through it
rather than passing down its middle: the centroids here sit 1.3 to 10.9 A off
their own trace. Worse, the two centroids are close enough together that the
nearest-point test can put the distal mark EARLIER along the path than the
proximal one, inverting the order the mechanism runs in. A tick drawn at a
centroid would assert a precision the geometry does not support.

So each pocket is measured as the extent it actually occupies: every one of
its lining residues is projected onto the trace, and the band runs from the
10th to the 90th percentile of those positions, which keeps a couple of
outlying residues from stretching it across the whole panel. The centroid
position and its offset are kept alongside, so the band can be checked against
the cruder measure.

No grid work is needed - the traces are already on disk, so this is a cheap
second pass over them.

Writes results/tables/pocket_marks.csv
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import per_structure_tunnels as pst
from published_pockets import PDBDIR
from mexb_common import (CXDIR, DBP, PBP, STRUCT_DIR, TABLES, Structure,
                         centroid, fmt, write_csv)


def main():
    out = []
    for (tab, kind) in (("full_tunnels.csv", "fixed"),
                        ("own_route_tunnels.csv", "seeded")):
        path = os.path.join(TABLES, tab)
        if not os.path.exists(path):
            continue
        for r in pst.rows_of(path):
            f = os.path.join(CXDIR, r["trace_file"])
            if not os.path.exists(f):
                continue
            P, rad = pst.read_trace(f)
            arc = np.concatenate([[0.0], np.cumsum(
                np.linalg.norm(np.diff(P, axis=0), axis=1))])
            sp = os.path.join(STRUCT_DIR, f"{r['pdb']}.pdb")
            if not os.path.exists(sp):
                sp = os.path.join(PDBDIR, f"{r['pdb']}.pdb")
            if not os.path.exists(sp):
                continue
            ca = Structure(sp).ca(r["chain"])
            for (name, res) in (("proximal", PBP), ("distal", DBP)):
                c = centroid(ca, res)
                d = np.linalg.norm(P - c, axis=1)
                i = int(d.argmin())
                X = np.array([ca[q] for q in res if q in ca])
                if len(X) < 4:
                    continue
                j = np.linalg.norm(X[:, None, :] - P[None, :, :],
                                   axis=2).argmin(1)
                a = np.sort(arc[j])
                lo, hi = (float(np.percentile(a, 10)),
                          float(np.percentile(a, 90)))
                out.append([r["pdb"], r["chain"], r["ligand"], kind, name,
                            len(X), fmt(100.0 * lo / arc[-1]),
                            fmt(100.0 * hi / arc[-1]),
                            fmt(100.0 * arc[i] / arc[-1]), fmt(float(d[i])),
                            fmt(float(arc[-1]))])
                print(f"  {r['ligand']:16} {kind:6} {name:8} band "
                      f"{100*lo/arc[-1]:5.1f}-{100*hi/arc[-1]:5.1f}%  "
                      f"(centroid {100*arc[i]/arc[-1]:5.1f}%, "
                      f"{d[i]:4.1f} A off)")
    write_csv(os.path.join(TABLES, "pocket_marks.csv"),
              ["pdb", "chain", "ligand", "route", "pocket", "lining_residues",
               "band_lo_pct", "band_hi_pct", "centroid_pct",
               "centroid_offset_A", "tunnel_length_A"], out)
    print(f"\nwrote results/tables/pocket_marks.csv ({len(out)})")


if __name__ == "__main__":
    main()
