#!/usr/bin/env python3
"""Diagnostic: how open is the switch-loop (F615) gate in each protomer?

PROTOCOL section 6 expects the ampicillin chain-E tunnel to be constricted at
F615 with a 2.01 A bottleneck. This pipeline's widest-path search instead
routes chain E around F615 and pinches at the PC2/DC gate. This script
measures the F615 gate itself, so the two statements can be compared.

For each chain it reports:
  - the widest clearance anywhere in the F615/F617/F136 gate neighbourhood
  - the bottleneck of a path forced to pass through that gate
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import DBP, PBP, TABLES, centroid, coords, fmt, \
    load_structures, write_csv
from tunnels import (ClearanceGrid, ExactClearance, bulk_region, densify,
                     refine_path, smooth_path, widest_path)

GATE_RES = [615, 617, 136]


def main():
    rows = []
    for s in load_structures():
        atoms = [a for a in s.protein_atoms if not a.is_hydrogen]
        grid = ClearanceGrid(atoms, step=0.6, verbose=False)
        cf = ExactClearance(atoms)
        bulk = bulk_region(grid.clearance)
        for ch in s.chains:
            sc = []
            for r in GATE_RES:
                sc += [a for a in s.residue_atoms(ch, r)
                       if a.name not in ("N", "CA", "C", "O")]
            if not sc:
                continue
            gate = coords(sc).mean(axis=0)
            gidx, gclr = grid.free_seed(gate, min_clear=0.4, search=10.0)
            if gidx is None:
                rows.append([s.name, ch, "", "", "no free point in the gate"])
                continue
            gate_pt = grid.point_of(gidx)
            # widest clearance anywhere within 7 A of the gate centre
            c = grid.index_of(gate)
            rad = int(np.ceil(7.0 / grid.step))
            sl = tuple(slice(max(0, c[i] - rad),
                             min(grid.shape[i], c[i] + rad + 1))
                       for i in range(3))
            widest = float(grid.clearance[sl].max())
            # path forced through the gate: site -> gate, gate -> bulk
            ca = s.ca(ch)
            site = 0.5 * (centroid(ca, DBP) + centroid(ca, PBP))
            sidx, _ = grid.free_seed(site)
            seg = []
            ok = True
            for a, b in ((sidx, gidx), (gidx, None)):
                if b is None:
                    R, path = widest_path(grid.clearance, a, bulk)
                else:
                    target = np.zeros(grid.shape, dtype=bool)
                    target[b] = True
                    R, path = widest_path(grid.clearance, a, target)
                if R is None:
                    ok = False
                    break
                raw = np.array([grid.point_of(q) for q in path])
                if len(raw) >= 5:
                    raw = densify(refine_path(smooth_path(raw, n=1), cf))
                seg.append(float(cf(raw).min()))
            via = min(seg) if ok and seg else None
            print(f"  {s.name} {ch}: gate widest {widest:.2f} A, "
                  f"forced-through-gate bottleneck "
                  f"{('%.2f A' % via) if via else 'no path'}")
            rows.append([s.name, ch, fmt(widest), fmt(via), ""])
    write_csv(os.path.join(TABLES, "switch_gate.csv"),
              ["structure", "chain", "gate_widest_clearance_A",
               "bottleneck_forced_through_gate_A", "note"], rows)
    print("\nwrote results/tables/switch_gate.csv")


if __name__ == "__main__":
    main()
