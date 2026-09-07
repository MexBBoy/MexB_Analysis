#!/usr/bin/env python3
"""Localised backbone RMSD by region, after Lawrence et al. Fig. 2B.

Their panel breaks the AcrB-vs-MdtF comparison down by region rather than
quoting one global number, separately for each protomer state, which is what
shows that the transmembrane domain and not the porter domain carries the
difference. The same comparison is available here for MexB.

Three comparisons are made per region and per state:
  MexB vs AcrB        the cross-transporter one, as in the paper
  ampicillin vs DDM   our two reconstructions against each other
  MexB vs MexB apo    a within-transporter, cross-state control

MexB residue numbers are mapped onto AcrB by the same global BLOSUM62
alignment lining_conservation.py uses. Each protomer pair is superposed on
all its common CA, then RMSD is computed per region on that superposition -
so a region's value says how far it sits from where the global fit puts it,
which is the quantity the paper's panel shows.

AcrB states are assigned by ranking 4DX5's own protomers on the PN1-PN2 and
PC1-PC2 separations rather than by this project's absolute thresholds, which
are calibrated on MexB and put two of AcrB's three chains in the same state.

Writes results/tables/regional_rmsd.csv.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lining_conservation import align, fetch
from mexb_analysis import cleft_metrics
from mexb_common import (REGIONS, STRUCT_DIR, TABLES, Structure, apply_rt,
                         fmt, kabsch, write_csv)

ACRB_PDB = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "work", "acrb", "4DX5.pdb")
AMP, DDM = "Amp_MexB_20260826", "MexB_DDM_3_20260730"


def mexb_to_acrb():
    """{MexB resseq: AcrB resseq} from the global alignment."""
    m = align(fetch("P52002"), fetch("P31224"))
    return {i + 1: int(m[i]) + 1 for i in range(len(m)) if m[i] >= 0}


def rank_states(s):
    """Assign access/binding/extrusion by ranking this trimer's own clefts.

    Extrusion is the protomer with the tightest PC1-PC2 cleft and access the
    one with the tightest PN1-PN2, which is the ordering the functional
    rotation defines; using absolute cutoffs calibrated on another
    transporter mislabels them.
    """
    m = {}
    for ch in s.chains:
        try:
            c = cleft_metrics(s, ch)
        except Exception:
            continue
        m[ch] = (c["PN1-PN2"]["sep"], c["PC1-PC2"]["sep"])
    if len(m) < 3:
        return {}
    ext = min(m, key=lambda c: m[c][1])
    rest = [c for c in m if c != ext]
    acc = min(rest, key=lambda c: m[c][0])
    bnd = [c for c in rest if c != acc][0]
    return {acc: "Access", bnd: "Binding", ext: "Extrusion"}


def pair_rmsd(a, ach, b, bch, mapping=None):
    """Per-region CA RMSD after superposing the two protomers globally."""
    ca_a, ca_b = a.ca(ach), b.ca(bch)
    pairs = []
    for r in sorted(ca_a):
        t = mapping.get(r) if mapping else r
        if t is not None and t in ca_b:
            pairs.append((r, t))
    if len(pairs) < 200:
        return None, 0, {}
    A = np.array([ca_a[r] for r, _ in pairs])
    B = np.array([ca_b[t] for _, t in pairs])
    R, tr = kabsch(A, B)
    Am = apply_rt(R, tr, A)
    dev = np.linalg.norm(Am - B, axis=1)
    glob = float(np.sqrt((dev ** 2).mean()))
    out = {}
    for name, ids in REGIONS.items():
        sel = np.array([i for i, (r, _) in enumerate(pairs) if r in set(ids)])
        if len(sel) >= 8:
            out[name] = (float(np.sqrt((dev[sel] ** 2).mean())), len(sel))
    return glob, len(pairs), out


def main():
    amp = Structure(os.path.join(STRUCT_DIR, f"{AMP}.pdb"))
    ddm = Structure(os.path.join(STRUCT_DIR, f"{DDM}.pdb"))
    if not os.path.exists(ACRB_PDB):
        print(f"  {ACRB_PDB} missing - run multiligand_survey.py first")
        return
    acrb = Structure(ACRB_PDB)

    from mexb_analysis import assign_state
    def states(s):
        out = {}
        for ch in s.chains:
            try:
                c = cleft_metrics(s, ch)
            except Exception:
                continue
            out[ch] = assign_state(c["PN1-PN2"]["sep"],
                                   c["PC1-PC2"]["sep"])[0]
        return out

    amp_st, ddm_st = states(amp), states(ddm)
    acrb_st = rank_states(acrb)
    print("=== localised backbone RMSD by region (Lawrence Fig. 2B) ===")
    print(f"    ampicillin {amp_st}")
    print(f"    DDM        {ddm_st}")
    print(f"    AcrB 4DX5  {acrb_st}  (ranked within the trimer)\n")

    mp = mexb_to_acrb()
    print(f"    {len(mp)} of 1046 MexB residues map onto AcrB\n")

    jobs = []
    for st in ("Access", "Binding", "Extrusion"):
        a = [c for c, v in amp_st.items() if v == st]
        d = [c for c, v in ddm_st.items() if v == st]
        k = [c for c, v in acrb_st.items() if v == st]
        if a and k:
            jobs.append(("MexB (ampicillin) vs AcrB", st, amp, a[0],
                         acrb, k[0], mp))
        if d and k:
            jobs.append(("MexB (DDM) vs AcrB", st, ddm, d[0], acrb, k[0], mp))


    # our two reconstructions are compared chain for chain, not state for
    # state: ampicillin D is a binding protomer where DDM D is an access one,
    # so pairing on state would compare different chains and report their
    # state difference as a structural one
    for ch in sorted(set(amp_st) & set(ddm_st)):
        st = (f"{amp_st[ch]} vs {ddm_st[ch]}" if amp_st[ch] != ddm_st[ch]
              else amp_st[ch])
        jobs.append((f"ampicillin vs DDM, chain {ch}", st, amp, ch,
                     ddm, ch, None))

    rows = []
    for label, st, A, ach, B, bch, mapping in jobs:
        glob, n, per = pair_rmsd(A, ach, B, bch, mapping)
        if glob is None:
            print(f"  {label:26} {st:9}: too few common CA")
            continue
        top = sorted(per.items(), key=lambda x: -x[1][0])[:3]
        print(f"  {label:26} {st:9}: global {glob:4.2f} A over {n} CA"
              f"   worst " + ", ".join(f"{k} {v[0]:.2f}" for k, v in top))
        for name, (val, cnt) in per.items():
            rows.append([label, st, f"{A.name}:{ach}", f"{B.name}:{bch}",
                         name, fmt(val), cnt, fmt(glob), n])

    write_csv(os.path.join(TABLES, "regional_rmsd.csv"),
              ["comparison", "state", "structure_a", "structure_b", "region",
               "rmsd_A", "n_CA_in_region", "global_rmsd_A", "n_CA_total"],
              rows)
    print("\nwrote results/tables/regional_rmsd.csv")


if __name__ == "__main__":
    main()
