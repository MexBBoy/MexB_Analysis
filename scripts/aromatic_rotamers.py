#!/usr/bin/env python3
"""Side-chain rotamers of the pocket aromatics, across every MexB protomer.

After Lawrence et al. (Nat Commun 2025;16:10601) Fig. 3C/3F, where a chi1
rotation of the DBP phenylalanines opens a 'hydrophobic nook' and, in the
wild type, contracts the pocket through steric clashes between F613, F626,
F610 and Y327.

P4 and P5 showed that pocket *volume* does not respond to the ligand, only to
the conformational state. This asks the other half of the question: does the
lining rearrange instead? A pocket can accommodate a different ligand by
rotating its side chains without changing size at all, and chi1 is where that
would show.

chi1 is N-CA-CB-CG, chi2 CA-CB-CG-CD1 (equivalent atoms for Trp/His). Both
come straight from the deposited coordinates, so this needs no superposition
and no sampling.

Writes results/tables/aromatic_rotamers.csv.
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_analysis import assign_state, cleft_metrics
from published_pockets import PDBDIR, ensure_pdbs, numbering_ok
from mexb_common import STRUCT_DIR, TABLES, Structure, fmt, write_csv

# the aromatic cluster that lines the porter pocket, distal set first
AROM = [136, 178, 327, 573, 610, 615, 617, 626, 628, 664, 666]
CHI1 = {"PHE": ("N", "CA", "CB", "CG"), "TYR": ("N", "CA", "CB", "CG"),
        "TRP": ("N", "CA", "CB", "CG"), "HIS": ("N", "CA", "CB", "CG")}
CHI2 = {"PHE": ("CA", "CB", "CG", "CD1"), "TYR": ("CA", "CB", "CG", "CD1"),
        "TRP": ("CA", "CB", "CG", "CD1"), "HIS": ("CA", "CB", "CG", "ND1")}


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w),
                                       np.dot(v, w))))


def chi(atoms, names):
    by = {a.name.strip(): np.asarray(a.xyz, float) for a in atoms}
    if not all(n in by for n in names):
        return None
    return dihedral(*(by[n] for n in names))


def rotamer(x):
    """chi1 well: the three staggered rotamers of a chi1 angle."""
    if x is None:
        return ""
    x = (x + 360) % 360
    return "g-" if 0 <= x < 120 else "t" if 120 <= x < 240 else "g+"


def circ_mean(a):
    a = np.radians(np.asarray(a, float))
    return float(np.degrees(np.arctan2(np.sin(a).mean(), np.cos(a).mean())))


def circ_sd(a):
    a = np.radians(np.asarray(a, float))
    r = np.hypot(np.sin(a).mean(), np.cos(a).mean())
    return float(np.degrees(np.sqrt(-2 * np.log(max(r, 1e-12)))))


def circ_diff(a, b):
    return float(abs((a - b + 180) % 360 - 180))


def main():
    ensure_pdbs()
    targets = [(os.path.join(STRUCT_DIR, f"{n}.pdb"), n)
               for n in ("Amp_MexB_20260826", "MexB_DDM_3_20260730")]
    targets += [(f, os.path.basename(f)[:-4])
                for f in sorted(glob.glob(os.path.join(PDBDIR, "*.pdb")))]
    OURS = ("Amp_MexB_20260826", "MexB_DDM_3_20260730")

    print("=== chi1 of the pocket aromatics, every protomer ===\n")
    rows = []
    for path, pid in targets:
        s = Structure(path)
        for ch in s.chains:
            good, bad = numbering_ok(s, ch)
            if bad or good < 4:
                continue
            m = cleft_metrics(s, ch)
            call, _, _, _ = assign_state(m["PN1-PN2"]["sep"],
                                         m["PC1-PC2"]["sep"])
            for rid in AROM:
                at = s.residue_atoms(ch, rid)
                if not at:
                    continue
                rn = at[0].resname
                if rn not in CHI1:
                    rows.append([pid, ch, call, rid, rn, "", "", "",
                                 "not aromatic"])
                    continue
                c1, c2 = chi(at, CHI1[rn]), chi(at, CHI2[rn])
                rows.append([pid, ch, call, rid, rn, fmt(c1), fmt(c2),
                             rotamer(c1), ""])

    write_csv(os.path.join(TABLES, "aromatic_rotamers.csv"),
              ["pdb", "chain", "state_call", "resseq", "resname",
               "chi1_deg", "chi2_deg", "chi1_rotamer", "note"], rows)

    # --- does the multi-ligand protomer sit outside the published spread?
    print("  Binding protomers only. Reference spread = published "
          "single-ligand structures.\n")
    print(f"  {'res':>6}  {'published chi1':>18}  {'rotamers':<12}"
          f"  {'DDM x3 E':>10}  {'dev':>6}   {'ampicillin E':>12}  {'dev':>6}")
    hits = []
    for rid in AROM:
        pub = [r for r in rows if r[2] == "Binding" and r[3] == rid
               and r[0] not in OURS and r[5]]
        if len(pub) < 4:
            continue
        ang = [float(r[5]) for r in pub]
        mu, sd = circ_mean(ang), circ_sd(ang)
        rots = sorted({r[7] for r in pub})
        line = f"  {pub[0][4]}{rid:<4}  {mu:>8.0f} +- {sd:<5.0f} (n={len(ang)})" \
               f"  {'/'.join(rots):<12}"
        for who in ("MexB_DDM_3_20260730", "Amp_MexB_20260826"):
            mine = [r for r in rows if r[0] == who and r[1] == "E"
                    and r[3] == rid and r[5]]
            if not mine:
                line += f"  {'-':>10}  {'-':>6}"
                continue
            v = float(mine[0][5])
            d = circ_diff(v, mu)
            line += f"  {v:>10.0f}  {d:>5.0f}"
            if who == "MexB_DDM_3_20260730" and d > max(2 * sd, 30):
                hits.append((rid, pub[0][4], v, mu, sd, d, mine[0][7],
                             "/".join(rots)))
        print(line)

    print("\n  --- where the DDM x3 protomer's lining differs ---")
    if not hits:
        print("    nothing outside 2 SD (or 30 deg) of the published "
              "Binding spread: the aromatics sit in the same rotamers.")
    for rid, rn, v, mu, sd, d, rot, rots in hits:
        print(f"    {rn}{rid}: chi1 {v:.0f} deg vs published {mu:.0f} "
              f"+- {sd:.0f} ({d:.0f} deg away); rotamer {rot}, "
              f"published {rots}")
    print("\nwrote results/tables/aromatic_rotamers.csv")


if __name__ == "__main__":
    main()
