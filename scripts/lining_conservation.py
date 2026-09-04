#!/usr/bin/env python3
"""Are the pocket-lining residues conserved scaffold, or specificity tuning?

P4/P5/P7 say the MexB pocket is pre-organised: it does not resize or rotate
its lining for a different ligand. P6 says different ligands dock at
different stations on that one surface. This asks what those stations are
built from - residues shared across the RND family, or residues that differ
between transporters with different substrate ranges.

Method. Each homologue is aligned to MexB (P52002) by global Needleman-Wunsch
with BLOSUM62 and affine gaps, and the residue aligned to each MexB pocket
position is read off. Pairwise-to-MexB rather than a multiple alignment
because only the MexB column matters here and pairwise alignment introduces
no ordering artefacts. The panel is chosen to contrast substrate ranges:
MexY and AcrD prefer aminoglycosides, which MexB and AcrB export poorly.

Writes results/tables/lining_conservation.csv.
"""
from __future__ import annotations

import json
import os
import sys
import urllib.request

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import (DBP, PBP, SWITCH_LOOP, TABLES, WORK_DIR, fmt,
                         write_csv)

PANEL = [("MexB", "P52002", "broad; beta-lactams, fluoroquinolones"),
         ("MexD", "Q9HVI9", "broad, P. aeruginosa"),
         ("MexF", "Q51396", "chloramphenicol, fluoroquinolones"),
         ("MexY", "G3XCW2", "AMINOGLYCOSIDES"),
         ("AcrB", "P31224", "broad, E. coli"),
         ("AcrF", "P24181", "broad, E. coli"),
         ("MdtF", "P37637", "broad; the Lawrence et al. transporter"),
         ("AcrD", "P24177", "AMINOGLYCOSIDES")]
AROM = {"F", "Y", "W", "H"}


def aligner():
    """Global BLOSUM62 aligner with affine gaps, from Biopython.

    Biopython's validated matrix and C implementation rather than a
    hand-entered matrix and a hand-rolled Needleman-Wunsch: the alignment is
    the one thing here that must not be subtly wrong.
    """
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    al.mode = "global"
    return al


AL = aligner()


def align(ref, other):
    """Index in `other` aligned to each index of `ref`; -1 for a gap."""
    aln = AL.align(ref, other)[0]
    m = np.full(len(ref), -1, dtype=int)
    for (r0, r1), (o0, o1) in zip(*aln.aligned):
        for k in range(r1 - r0):
            m[r0 + k] = o0 + k
    return m


def fetch(acc):
    """UniProt sequence, cached under work/."""
    os.makedirs(os.path.join(WORK_DIR, "seqs"), exist_ok=True)
    p = os.path.join(WORK_DIR, "seqs", f"{acc}.txt")
    if os.path.exists(p):
        return open(p).read().strip()
    d = json.load(urllib.request.urlopen(
        f"https://rest.uniprot.org/uniprotkb/{acc}.json"))
    s = d["sequence"]["value"]
    open(p, "w").write(s)
    return s


def main():
    ref = fetch(PANEL[0][1])
    print(f"=== conservation of the pocket lining across {len(PANEL)} RND "
          f"transporters ===")
    print(f"    MexB {len(ref)} aa; zero numbering offset, so residue n is "
          f"sequence index n-1\n")

    maps = {}
    for name, acc, note in PANEL[1:]:
        s = fetch(acc)
        m = align(ref, s)
        ident = sum(1 for i, j in enumerate(m)
                    if j >= 0 and ref[i] == s[j]) / len(ref) * 100
        maps[name] = (s, m)
        print(f"  {name:5} {acc}  {len(s):4} aa  {ident:4.1f}% identity to "
              f"MexB   {note}")

    # F664 and F666 are contacted by the outermost DDM but are not in the
    # PROTOCOL lining lists, so they are added here as "contacted"
    EXTRA = [664, 666]
    sites = {"DBP": DBP, "PBP": PBP, "switch": SWITCH_LOOP,
             "contacted": EXTRA}
    lining = sorted(set(DBP) | set(PBP) | set(SWITCH_LOOP) | set(EXTRA))
    rows = []
    print(f"\n  {'res':>7} {'site':<8} " +
          " ".join(f"{n:<5}" for n, _, _ in PANEL[1:]) + "  identity")
    for r in lining:
        i = r - 1
        if not (0 <= i < len(ref)):
            continue
        wt = ref[i]
        others, same = [], 0
        for name, _, _ in PANEL[1:]:
            s, m = maps[name]
            j = m[i]
            c = s[j] if j >= 0 else "-"
            others.append(c)
            if c == wt:
                same += 1
        pct = same / len(others) * 100
        where = ";".join(k for k, v in sites.items() if r in v)
        print(f"  {wt}{r:<6} {where:<8} " +
              " ".join(f"{c:<5}" for c in others) + f"  {pct:5.0f}%")
        rows.append([r, wt, where, "".join(others), fmt(pct, 0),
                     "yes" if wt in AROM else "no",
                     sum(1 for c in others if c in AROM)])

    write_csv(os.path.join(TABLES, "lining_conservation.csv"),
              ["resseq", "mexb_residue", "site", "homologues_" +
               "".join(n[-1] for n, _, _ in PANEL[1:]), "percent_identical",
               "mexb_aromatic", "homologues_aromatic"], rows)

    print("\n  --- summary ---")
    for lab, ids in (("DBP", DBP), ("PBP", PBP), ("switch loop", SWITCH_LOOP)):
        v = [float(r[4]) for r in rows if r[0] in ids]
        print(f"    {lab:12} mean identity {np.mean(v):5.1f}%  "
              f"({sum(1 for x in v if x == 100)}/{len(v)} invariant)")
    arom = [r for r in rows if r[5] == "yes"]
    v = [float(r[4]) for r in arom]
    print(f"    aromatics    mean identity {np.mean(v):5.1f}%  "
          f"({sum(1 for x in v if x == 100)}/{len(v)} invariant)")
    nonarom = [float(r[4]) for r in rows if r[5] == "no"]
    print(f"    non-aromatic mean identity {np.mean(nonarom):5.1f}%")

    # how many of MexB's pocket aromatics does each homologue keep aromatic?
    print("\n  --- aromatic scaffold retained, per transporter ---")
    idx = [i for i, r in enumerate(rows) if r[5] == "yes"]
    for k, (name, _, note) in enumerate(PANEL[1:]):
        keep = sum(1 for i in idx if rows[i][3][k] in AROM)
        exact = sum(1 for i in idx if rows[i][3][k] == rows[i][1])
        print(f"    {name:5} {exact}/{len(idx)} identical, "
              f"{keep}/{len(idx)} still aromatic   {note}")

    # the three stations of the DDM x3 structure (from P6)
    print("\n  --- conservation of the three DDM stations ---")
    STATION = {"outermost (33.7 A)": [136, 573, 617, 628, 664, 666, 327],
               "middle (50.7 A)": [615, 617],
               "deepest (62.1 A)": [178, 610, 615, 628]}
    look = {r[0]: r for r in rows}
    for lab, ids in STATION.items():
        v = [float(look[i][4]) for i in ids if i in look]
        got = [f"{look[i][1]}{i} {float(look[i][4]):.0f}%"
               for i in ids if i in look]
        print(f"    {lab:20} mean identity {np.mean(v):5.1f}%   "
              + ", ".join(got))
    print("\nwrote results/tables/lining_conservation.csv")


if __name__ == "__main__":
    main()
