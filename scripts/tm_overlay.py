#!/usr/bin/env python3
"""Cartoon overlays of the MexB and AcrB transmembrane domains, per state.

The equivalent of Lawrence et al. Fig. 2C: the two transporters' TM domains
superposed and drawn as cartoons for each protomer state, split into the two
pseudo-symmetric repeats R1 (TM1-6) and R2 (TM7-12), so the helical
rearrangement through the functional rotation can be read off.

Every panel is rendered from one camera. Each MexB protomer is first
superposed on the AcrB protomer of the same state over their shared TM CA,
then that rigid pair is moved onto the AcrB access protomer, so all six
images share a frame and differences between panels are real.

Needs PyMOL (pip install pymol-open-source). Writes the panels to
results/figures/tm_overlay/ and a manifest to results/tables/tm_overlay.csv.
"""
from __future__ import annotations

import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from regional_rmsd import ACRB_PDB, AMP, DDM, mexb_to_acrb, rank_states
from mexb_analysis import assign_state, cleft_metrics
from mexb_common import (FIGURES, REGIONS, STRUCT_DIR, TABLES, WORK_DIR,
                         Structure, apply_rt, fmt, kabsch, write_csv)

R1 = sorted(set(REGIONS["TM1"]) | set(REGIONS["TM2"]) | set(REGIONS["Ialpha"])
            | set(REGIONS["TM3-6"]) | set(REGIONS["TM6b"]))
# TM7 starts at ~861, before the TM7-12 region boundary, so R2 is extended
# back through the junction to keep that helix whole
R2 = sorted(set(REGIONS["junction859-875"]) | set(REGIONS["TM7-12"]))
TM = sorted(set(R1) | set(R2))
OUT = os.path.join(FIGURES, "tm_overlay")
WORK = os.path.join(WORK_DIR, "tmoverlay")
MEXB_COL, ACRB_COL = "0x2E5FE8", "0xE59BD8"   # blue / pink, as the paper
SWING_COL = "0xD11149"                        # the helix that moves most


def states_of(s):
    out = {}
    for ch in s.chains:
        try:
            c = cleft_metrics(s, ch)
        except Exception:
            continue
        out[ch] = assign_state(c["PN1-PN2"]["sep"], c["PC1-PC2"]["sep"])[0]
    return out


def write_pdb(path, atoms, xyz):
    with open(path, "w") as fh:
        for i, (a, p) in enumerate(zip(atoms, xyz), 1):
            nm = a.name if len(a.name) >= 4 else f" {a.name:<3.3s}"
            fh.write(f"ATOM  {i:5d} {nm:<4.4s} {a.resname:>3.3s} "
                     f"A{a.resseq:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}"
                     f"  1.00  0.00\n")
        fh.write("END\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    if not os.path.exists(ACRB_PDB):
        print("  AcrB 4DX5 missing - run multiligand_survey.py first")
        return
    acrb = Structure(ACRB_PDB)
    ddm = Structure(os.path.join(STRUCT_DIR, f"{DDM}.pdb"))
    acrb_st = rank_states(acrb)
    ddm_st = states_of(ddm)
    mp = mexb_to_acrb()
    print("=== TM cartoon overlays, MexB (DDM x3 model) on AcrB 4DX5 ===")
    print(f"    MexB {ddm_st}\n    AcrB {acrb_st}\n")

    def ca_pairs(mch, ach, ids):
        mca, aca = ddm.ca(mch), acrb.ca(ach)
        out = []
        for r in ids:
            t = mp.get(r)
            if r in mca and t in aca:
                out.append((r, t))
        return out

    # common view frame: the AcrB access protomer
    view_ch = [c for c, v in acrb_st.items() if v == "Access"][0]
    rows, objs = [], []
    for st in ("Access", "Binding", "Extrusion"):
        m = [c for c, v in ddm_st.items() if v == st]
        a = [c for c, v in acrb_st.items() if v == st]
        if not m or not a:
            print(f"  {st}: no matching protomer pair"); continue
        mch, ach = m[0], a[0]

        # MexB onto its AcrB partner, over the shared TM
        pr = ca_pairs(mch, ach, TM)
        Rm, tm_ = kabsch(np.array([ddm.ca(mch)[r] for r, _ in pr]),
                         np.array([acrb.ca(ach)[t] for _, t in pr]))
        fit = float(np.sqrt(((apply_rt(Rm, tm_, np.array(
            [ddm.ca(mch)[r] for r, _ in pr]))
            - np.array([acrb.ca(ach)[t] for _, t in pr])) ** 2).sum(1).mean()))

        # that pair onto the AcrB access protomer, for a shared camera
        pv = [(t, mp_t) for t, mp_t in
              [(r, mp.get(r)) for r in TM]
              if mp_t in acrb.ca(ach) and mp_t in acrb.ca(view_ch)]
        Rv, tv = kabsch(np.array([acrb.ca(ach)[t] for _, t in pv]),
                        np.array([acrb.ca(view_ch)[t] for _, t in pv]))

        for rep, ids in (("R1", R1), ("R2", R2)):
            keep = set(ids)
            ma = [x for x in ddm.protein_atoms
                  if x.chain == mch and not x.is_hydrogen
                  and x.resseq in keep]
            aset = {mp[r] for r in ids if r in mp}
            aa = [x for x in acrb.protein_atoms
                  if x.chain == ach and not x.is_hydrogen
                  and x.resseq in aset]
            mx = apply_rt(Rv, tv, apply_rt(
                Rm, tm_, np.array([x.xyz for x in ma])))
            ax = apply_rt(Rv, tv, np.array([x.xyz for x in aa]))
            for tag, at, xyz in (("mexb", ma, mx), ("acrb", aa, ax)):
                f = os.path.join(WORK, f"{st}_{rep}_{tag}.pdb")
                write_pdb(f, at, xyz)
                objs.append((st, rep, tag, f))
            rows.append([st, rep, f"{DDM}:{mch}", f"4DX5:{ach}", len(pr),
                         fmt(fit), len(ma), len(aa)])
        print(f"  {st:9}: MexB {mch} on AcrB {ach}, TM fit {fit:.2f} A "
              f"over {len(pr)} CA")

    write_csv(os.path.join(TABLES, "tm_overlay.csv"),
              ["state", "repeat", "mexb", "acrb", "n_TM_CA_fitted",
               "tm_fit_rmsd_A", "mexb_atoms", "acrb_atoms"], rows)


    # the helix that swings furthest between states, measured rather than
    # picked by eye - helix_displacement.py writes it
    swing = None
    hd = os.path.join(TABLES, "helix_displacement.csv")
    if os.path.exists(hd):
        import csv as _csv
        cand = [r for r in _csv.DictReader(open(hd))
                if r["comparison"].startswith("MexB") and "vs Access"
                in r["comparison"] and r["helix"]]
        if cand:
            b = max(cand, key=lambda r: float(r["centroid_shift_A"]))
            swing = (b["helix"], int(b["first_res"]), int(b["last_res"]),
                     float(b["centroid_shift_A"]), b["state"])
            print(f"\n  highlighting {swing[0]} ({swing[1]}-{swing[2]}), "
                  f"which moves {swing[3]:.1f} A from access to "
                  f"{swing[4].lower()}")

    # one PyMOL session, one camera per repeat
    pml = os.path.join(WORK, "render.pml")
    with open(pml, "w") as fh:
        fh.write("set ray_opaque_background, 0\nset cartoon_transparency, 0\n"
                 "set ray_shadows, 0\nset antialias, 2\n"
                 "set cartoon_cylindrical_helices, 1\n"
                 "set cartoon_helix_radius, 2.2\nbg_color white\n")
        for st, rep, tag, f in objs:
            fh.write(f"load {f}, {st}_{rep}_{tag}\n")
        fh.write("set orthoscopic, 1\nhide everything\nshow cartoon\n")
        for st, rep, tag, f in objs:
            col = MEXB_COL if tag == "mexb" else ACRB_COL
            fh.write(f"color {col}, {st}_{rep}_{tag}\n")
        if swing:
            _, lo, hi, _, _ = swing
            for st, rep, tag, f in objs:
                if tag == "mexb":
                    fh.write(f"color {SWING_COL}, {st}_{rep}_{tag} "
                             f"and resi {lo}-{hi}\n")
        # one rotation for all six panels: orient on the whole TM domain so
        # the membrane normal is vertical, then only re-centre per repeat
        fh.write("show cartoon\norient all\nturn z, 90\n")
        for rep in ("R1", "R2"):
            fh.write(f"zoom *_{rep}_*, 6\n")
            for st in ("Access", "Binding", "Extrusion"):
                fh.write(f"hide everything\n"
                         f"show cartoon, {st}_{rep}_mexb {st}_{rep}_acrb\n"
                         f"png {os.path.join(OUT, f'{st}_{rep}.png')}, "
                         f"width=900, height=1300, dpi=300, ray=1\n")
        fh.write("quit\n")
    r = subprocess.run(["pymol", "-cq", pml], capture_output=True, text=True)
    made = sorted(f for f in os.listdir(OUT) if f.endswith(".png"))
    print(f"\n  rendered {len(made)} panels: {', '.join(made)}")
    if not made:
        print((r.stderr or r.stdout)[-500:])
    print(f"\nwrote {OUT} and results/tables/tm_overlay.csv")


if __name__ == "__main__":
    main()
