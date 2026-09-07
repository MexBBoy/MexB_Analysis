#!/usr/bin/env python3
"""Per-helix movement in the transmembrane domain.

Lawrence et al. Fig. 2C annotates specific helices as swinging out or staying
fixed. This measures which ones actually do, so a label on the equivalent
MexB panel rests on a number rather than on eyeballing the cartoon.

Helix boundaries are not hard-coded: PyMOL's dss assigns secondary structure
on the MexB model itself and every helical run of at least MIN_LEN residues
inside the transmembrane windows is taken. Canonical TM names are attached
only where a segment matches the published AcrB topology closely; the rest
are identified by residue range, since inventing a number would be worse
than leaving it out.

Two movements are reported per helix:
  vs AcrB     after superposing the whole TM domain of the same state, how
              far that helix sits from its AcrB counterpart
  vs Access   within MexB, how far it moves between protomer states - the
              swing the functional rotation produces

Each is given as a centroid displacement and as the angle between helix axes,
because a helix can slide without tilting or tilt without sliding.

Writes results/tables/helix_displacement.csv.
"""
from __future__ import annotations

import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from regional_rmsd import ACRB_PDB, DDM, mexb_to_acrb, rank_states
from tm_overlay import TM, states_of
from mexb_common import (STRUCT_DIR, TABLES, WORK_DIR, Structure, apply_rt,
                         fmt, kabsch, write_csv)

WINDOWS = [(5, 40), (330, 575), (855, 1040)]
MIN_LEN = 12
# AcrB topology, used only to name segments that clearly correspond
CANON = {(5, 40): "TM1", (330, 361): "TM2", (362, 390): "Iα",
         (391, 425): "TM3", (426, 460): "TM4", (461, 500): "TM5",
         (505, 545): "TM6", (546, 575): "TM6b", (855, 893): "TM7",
         (894, 921): "TM8", (922, 956): "TM9", (957, 991): "TM10",
         (992, 1040): "TM11/12"}


def helices(pdb, chain):
    """[(first, last)] helical runs inside the TM windows, from PyMOL dss."""
    pml = os.path.join(WORK_DIR, "dss.pml")
    with open(pml, "w") as fh:
        fh.write(f"load {pdb}, m\ncreate p, m and chain {chain} and polymer\n"
                 f"delete m\ndss p\n"
                 f"iterate p and name CA and ss H, print('H', resi)\nquit\n")
    r = subprocess.run(["pymol", "-cq", pml], capture_output=True, text=True)
    res = sorted({int(l.split()[1]) for l in r.stdout.splitlines()
                  if l.startswith("H ") and l.split()[1].lstrip("-").isdigit()})
    if not res:
        return []
    segs, s, p = [], res[0], res[0]
    for x in res[1:]:
        if x == p + 1:
            p = x
        else:
            segs.append((s, p)); s = p = x
    segs.append((s, p))
    return [(a, b) for a, b in segs if b - a + 1 >= MIN_LEN
            and any(a >= lo and b <= hi for lo, hi in WINDOWS)]


def name_of(a, b):
    mid = (a + b) / 2
    for (lo, hi), nm in CANON.items():
        if lo <= mid <= hi:
            return nm
    return ""


def axis(pts):
    """Unit vector along a helix, from the principal axis of its CA."""
    c = pts.mean(0)
    _, _, vt = np.linalg.svd(pts - c)
    v = vt[0]
    return c, v / np.linalg.norm(v)


def compare(pa, pb):
    """(centroid shift, axis angle) between two sets of matched CA."""
    ca, va = axis(pa)
    cb, vb = axis(pb)
    ang = np.degrees(np.arccos(min(1.0, abs(float(np.dot(va, vb))))))
    return float(np.linalg.norm(ca - cb)), float(ang)


def main():
    ddm_path = os.path.join(STRUCT_DIR, f"{DDM}.pdb")
    ddm = Structure(ddm_path)
    acrb = Structure(ACRB_PDB)
    ddm_st, acrb_st = states_of(ddm), rank_states(acrb)
    mp = mexb_to_acrb()

    ref_ch = [c for c, v in ddm_st.items() if v == "Binding"][0]
    hel = helices(ddm_path, ref_ch)
    print(f"=== per-helix movement, {len(hel)} TM helices from dss ===\n")
    for a, b in hel:
        print(f"    {name_of(a, b) or '(unnamed)':8} {a}-{b}")

    rows = []

    def cas(s, ch, ids, mapping=None):
        c = s.ca(ch)
        out = []
        for r in ids:
            t = mapping.get(r) if mapping else r
            if t is not None and t in c:
                out.append((r, c[t]))
        return out

    # --- MexB vs AcrB, per state, after superposing the whole TM domain
    print("\n  --- MexB vs AcrB, after whole-TM superposition ---")
    for st in ("Access", "Binding", "Extrusion"):
        m = [c for c, v in ddm_st.items() if v == st]
        a = [c for c, v in acrb_st.items() if v == st]
        if not m or not a:
            continue
        mch, ach = m[0], a[0]
        pr = [(r, mp[r]) for r in TM
              if r in ddm.ca(mch) and mp.get(r) in acrb.ca(ach)]
        R, t = kabsch(np.array([ddm.ca(mch)[r] for r, _ in pr]),
                      np.array([acrb.ca(ach)[u] for _, u in pr]))
        best = []
        for (h0, h1) in hel:
            ids = range(h0, h1 + 1)
            A = cas(ddm, mch, ids)
            B = cas(acrb, ach, ids, mp)
            common = sorted(set(r for r, _ in A) & set(r for r, _ in B))
            if len(common) < 8:
                continue
            PA = apply_rt(R, t, np.array([dict(A)[r] for r in common]))
            PB = np.array([dict(B)[r] for r in common])
            shift, ang = compare(PA, PB)
            rmsd = float(np.sqrt(((PA - PB) ** 2).sum(1).mean()))
            rows.append(["MexB vs AcrB", st, name_of(h0, h1), h0, h1,
                         len(common), fmt(shift), fmt(ang, 1), fmt(rmsd)])
            best.append((shift, ang, name_of(h0, h1) or f"{h0}-{h1}"))
        best.sort(key=lambda x: -x[0])
        print(f"    {st:9} " + ",  ".join(
            f"{n} {s:.1f} A/{a:.0f}°" for s, a, n in best[:3]))

    # --- within MexB, how far each helix swings between states
    print("\n  --- within MexB, movement from the access protomer ---")
    acc = [c for c, v in ddm_st.items() if v == "Access"]
    if acc:
        a0 = acc[0]
        for st in ("Binding", "Extrusion"):
            m = [c for c, v in ddm_st.items() if v == st]
            if not m:
                continue
            mch = m[0]
            pr = [r for r in TM if r in ddm.ca(a0) and r in ddm.ca(mch)]
            R, t = kabsch(np.array([ddm.ca(mch)[r] for r in pr]),
                          np.array([ddm.ca(a0)[r] for r in pr]))
            best = []
            for (h0, h1) in hel:
                common = [r for r in range(h0, h1 + 1)
                          if r in ddm.ca(a0) and r in ddm.ca(mch)]
                if len(common) < 8:
                    continue
                PA = apply_rt(R, t, np.array([ddm.ca(mch)[r] for r in common]))
                PB = np.array([ddm.ca(a0)[r] for r in common])
                shift, ang = compare(PA, PB)
                rmsd = float(np.sqrt(((PA - PB) ** 2).sum(1).mean()))
                rows.append([f"MexB {st} vs Access", st, name_of(h0, h1),
                             h0, h1, len(common), fmt(shift), fmt(ang, 1),
                             fmt(rmsd)])
                best.append((shift, ang, name_of(h0, h1) or f"{h0}-{h1}"))
            best.sort(key=lambda x: -x[0])
            print(f"    {st:9} " + ",  ".join(
                f"{n} {s:.1f} A/{a:.0f}°" for s, a, n in best[:3]))

    write_csv(os.path.join(TABLES, "helix_displacement.csv"),
              ["comparison", "state", "helix", "first_res", "last_res",
               "n_CA", "centroid_shift_A", "axis_angle_deg", "rmsd_A"], rows)
    print("\nwrote results/tables/helix_displacement.csv")


if __name__ == "__main__":
    main()
