#!/usr/bin/env python3
"""Every tunnel cut at the same exit, so the lengths can be compared.

The tunnel search stops at the first voxel of the bulk-connected open region -
the moment a route joins solvent that reaches the outside of the box. That is
a uniform rule, but it does not stop the seven routes at comparable places:
measured on their own trimer, ampicillin's route ends with 31 protein heavy
atoms within 12 A of its last point while chloramphenicol's still has 203. One
route has come out into the open; the other has only just broken the surface.
Raw path lengths therefore mix two things - how far the route runs inside the
protein, and how far past the surface the search happened to carry it.

So this re-cuts every trace at a matched enclosure. Enclosure is the number of
protein heavy atoms within 12 A of a trace point, counted on that protomer's
own trimer (the chain plus its two nearest neighbours, so a second trimer in
the asymmetric unit does not count as burial). Hydrogens are excluded: our
models carry them and the deposited ones do not, which would otherwise double
the count. The threshold is the largest end-of-trace enclosure in the set, so
every route reaches it. Each trace is cut at its point of no return - the last
point beyond which it never re-buries above the threshold - because two routes
(3W9J, 21FO) dip into an open vestibule and then run back into the protein.

Writes:
  results/tables/common_exit_profiles.csv   enclosure and radius along each
  results/tables/common_exit_summary.csv    one row per protomer
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import per_structure_tunnels as pst
from published_pockets import PDBDIR, load_channel
from mexb_common import STRUCT_DIR, TABLES, Structure, coords, fmt, write_csv

SHELL = 12.0        # enclosure counting radius, A
END_TRIM = 3.0      # same convention as the route bottleneck elsewhere


def enclosure(s, chain, P):
    """Protein heavy atoms within SHELL of each trace point, own trimer only."""
    heavy = [a for a in s.protein_atoms
             if (a.element or "").strip().upper() != "H"]
    cen = {c: coords([a for a in heavy
                      if a.chain == c and a.name.strip() == "CA"]).mean(0)
           for c in s.chains}
    near = sorted(cen, key=lambda c: float(np.linalg.norm(cen[c] - cen[chain])))[:3]
    tree = cKDTree(coords([a for a in heavy if a.chain in near]))
    return np.array([len(tree.query_ball_point(q, SHELL)) for q in P])


def neck(rad, depth):
    """Narrowest radius away from either end, the panel's bottleneck rule."""
    mid = (depth >= END_TRIM) & (depth <= depth[0] - END_TRIM)
    return float(rad[mid].min()) if mid.any() else float(rad.min())


def main():
    chan = load_channel()
    if chan is None:
        print("  reference channel missing - run tunnels.py first")
        return
    want = pst.panel_protomers()
    print("=== every tunnel cut at a matched enclosure ===")

    traces = {}
    for (pid, ch), nm in sorted(want.items(), key=lambda x: x[1]):
        path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            path = os.path.join(PDBDIR, f"{pid}.pdb")
        if not os.path.exists(path):
            print(f"  {nm}: {pid} not found - skipped")
            continue
        s = Structure(path)
        tr = pst.find_trace(s, ch, chan)
        if tr is None:
            print(f"  {nm}: no trace for chain {ch} - skipped")
            continue
        P, rad = pst.read_trace(tr)
        arc = np.concatenate([[0.0], np.cumsum(
            np.linalg.norm(np.diff(P, axis=0), axis=1))])
        traces[(pid, ch)] = (nm, arc, rad, enclosure(s, ch, P))

    if not traces:
        return
    thresh = max(int(e[-1]) for _, _, _, e in traces.values())
    print(f"    enclosure = protein heavy atoms within {SHELL:.0f} A, own "
          f"trimer\n    threshold {thresh}, the largest end-of-trace "
          f"enclosure in the set\n")

    prof, summ = [], []
    for (pid, ch), (nm, arc, rad, enc) in sorted(traces.items(),
                                                 key=lambda x: x[1][0]):
        # point of no return: first index from which enclosure never again
        # rises above the threshold
        i = int(np.argmax(np.maximum.accumulate(enc[::-1])[::-1] <= thresh))
        first = int(np.argmax(enc <= thresh))
        d = arc[i] - arc[:i + 1]                 # depth from the cut mouth
        for j in range(len(arc)):
            prof.append([pid, ch, nm, fmt(arc[j]), fmt(rad[j]), int(enc[j]),
                         "yes" if j <= i else "no"])
        summ.append([pid, ch, nm, fmt(arc[-1]), fmt(arc[i]), fmt(arc[first]),
                     thresh, int(enc[0]), int(enc[-1]),
                     fmt(neck(rad[:i + 1], d))])
        print(f"  {nm:16} {pid} {ch}: full {arc[-1]:6.1f} A, cut at "
              f"{arc[i]:6.1f} A (first crossing {arc[first]:6.1f}), "
              f"enclosure {enc[0]:4d} at the ligand, {enc[-1]:4d} at its "
              f"own end")

    write_csv(os.path.join(TABLES, "common_exit_profiles.csv"),
              ["pdb", "chain", "ligand", "depth_from_ligand_A", "radius_A",
               f"enclosure_atoms_{SHELL:.0f}A", "kept"], prof)
    write_csv(os.path.join(TABLES, "common_exit_summary.csv"),
              ["pdb", "chain", "ligand", "full_length_A", "matched_length_A",
               "first_crossing_A", "enclosure_threshold",
               "enclosure_at_ligand", "enclosure_at_own_end",
               "route_bottleneck_matched_A"], summ)
    print("\nwrote results/tables/common_exit_profiles.csv and "
          "common_exit_summary.csv")


if __name__ == "__main__":
    main()
