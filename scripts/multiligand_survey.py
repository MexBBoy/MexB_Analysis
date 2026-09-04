#!/usr/bin/env python3
"""Is more than one ligand per protomer actually unprecedented in RND pumps?

The poster's conclusion claims three ligands in one binding pocket have never
been seen in an RND inner-membrane component. That claim is checkable, and it
does not survive: AcrB 4DX7 holds two doxorubicin molecules 5.7 A apart in one
protomer. This script measures multi-ligand occupancy in the AcrB entries that
carry drugs, so the distinction the DDM x3 model can defend is quantified
rather than asserted.

What separates our structure is not the count but the arrangement. AcrB's pair
is a stack at one site. The three DDM molecules form a contiguous chain, each
within van der Waals contact of the next, spanning the proximal-to-distal path.

Writes results/tables/multiligand_survey.csv.
"""
from __future__ import annotations

import os
import sys
import urllib.request
from itertools import combinations

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from published_pockets import CRYO, PDBDIR, pocket_ligands
from mexb_common import STRUCT_DIR, TABLES, Structure, coords, fmt, write_csv

# AcrB entries carrying drugs or inhibitors, plus the classic apo references
ACRB = ["1IWG", "1T9U", "2DRD", "2DR6", "3AOA", "3AOB", "3AOC", "3AOD",
        "4DX5", "4DX7", "4ZLJ", "5ENO", "5ENP", "5ENQ", "5ENR", "5ENS"]
MIN_HEAVY = 10


def fetch(pid, into):
    os.makedirs(into, exist_ok=True)
    p = os.path.join(into, f"{pid}.pdb")
    if not os.path.exists(p) or os.path.getsize(p) < 10000:
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{pid}.pdb", p)
        except Exception as exc:
            print(f"  could not fetch {pid}: {exc}")
            return None
    return p


def groups(s):
    """{chain: [(resname, coords)]} for non-cryoprotectant ligands."""
    out = {}
    for (ch, rs, rn, ats) in s.ligands():
        h = [a for a in ats if not a.is_hydrogen]
        if rn in CRYO or len(h) < MIN_HEAVY:
            continue
        out.setdefault(ch, []).append((rn, coords(h)))
    return out


def main():
    acr = os.path.join(os.path.dirname(PDBDIR), "acrb")
    rows = []
    print("=== multi-ligand occupancy of one protomer ===")
    print("    'chain' = contiguous if every molecule is within 4.5 A of "
          "another in the set\n")

    targets = [(os.path.join(STRUCT_DIR, "MexB_DDM_3_20260730.pdb"),
                "MexB_DDM_3 (this work)", "E")]
    for pid in ACRB:
        p = fetch(pid, acr)
        if p:
            targets.append((p, f"AcrB {pid}", None))

    for path, lab, only in targets:
        s = Structure(path)
        for ch, items in sorted(groups(s).items()):
            if only and ch != only:
                continue
            if only:                       # our model: porter-pocket ligands
                items = [(rn, coords(h))
                         for c, rn, h in pocket_ligands(s) if c == ch]
            if len(items) < 2:
                continue
            n = len(items)
            cen = [x[1].mean(0) for x in items]
            trees = [cKDTree(x[1]) for x in items]
            pairs = []
            adj = {i: set() for i in range(n)}
            for i, j in combinations(range(n), 2):
                d_cen = float(np.linalg.norm(cen[i] - cen[j]))
                d_surf = float(trees[i].query(items[j][1])[0].min())
                pairs.append((i, j, d_cen, d_surf))
                if d_surf <= 4.5:
                    adj[i].add(j); adj[j].add(i)
            # is the whole set one contiguous run of touching molecules?
            seen, stack = {0}, [0]
            while stack:
                k = stack.pop()
                for m in adj[k] - seen:
                    seen.add(m); stack.append(m)
            contiguous = len(seen) == n
            spread = max(p[2] for p in pairs)
            nearest = min(p[3] for p in pairs)
            names = "+".join(sorted({x[0] for x in items}))
            print(f"  {lab:22} {ch}  {n} ligands ({names:12})  "
                  f"centroids {min(p[2] for p in pairs):5.1f}-{spread:5.1f} A"
                  f"  closest atoms {nearest:4.1f} A  "
                  + ("CONTIGUOUS CHAIN" if contiguous else "separate sites"))
            rows.append([lab, ch, n, names, fmt(spread), fmt(nearest),
                         "yes" if contiguous else "no"])

    write_csv(os.path.join(TABLES, "multiligand_survey.csv"),
              ["structure", "chain", "n_ligands", "ligands",
               "max_centroid_separation_A", "closest_approach_A",
               "one_contiguous_chain"], rows)
    print("\nwrote results/tables/multiligand_survey.csv")


if __name__ == "__main__":
    main()
