#!/usr/bin/env python3
"""What the two pockets are made of, distal against proximal.

Three independent readings of the same question, because any one of them can
be argued with:

  composition  - what the residues are, by class, straight from the sequence
  hydropathy   - Kyte-Doolittle mean over those residues
  lipophilic   - the fraction of side-chain heavy atoms that are carbon or
                 sulphur, measured on the structures rather than the sequence,
                 which is what a substrate's van der Waals surface actually
                 meets
and, as an observation rather than a property of the protein, the apolar
fraction of the contacts the bound ligands actually make in each site.

Writes results/tables/pocket_chemistry.csv (per protomer) and prints the
summary across protomers.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from per_structure_tunnels import rows_of
from published_pockets import PDBDIR
from mexb_common import (DBP, PBP, STRUCT_DIR, TABLES, Structure, coords, fmt,
                         write_csv)

# Kyte-Doolittle
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
      "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
      "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
         "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
         "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
         "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
AROM = set("FWY")
ALIPH = set("AVLIMPG")
POLAR = set("STNQCH")
CHARGED = set("DEKR")
PLUS, MINUS = set("KR"), set("DE")
BACKBONE = {"N", "CA", "C", "O", "OXT"}


def sidechain(atoms):
    return [a for a in atoms
            if a.name.strip() not in BACKBONE and not a.is_hydrogen]


def main():
    pids = []
    for base in (STRUCT_DIR, PDBDIR):
        if os.path.isdir(base):
            pids += [(f[:-4], os.path.join(base, f))
                     for f in sorted(os.listdir(base)) if f.endswith(".pdb")]
    seen, files = set(), []
    for pid, p in pids:
        if pid not in seen:
            seen.add(pid)
            files.append((pid, p))

    # a MexB chain, not any chain that happens to have those residue numbers:
    # 22XK and 22XM are complexes whose other components number into the same
    # range and would otherwise be scored as pockets
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    rres = {}
    for a in ref.protein_atoms:
        if a.chain == "E":
            rres.setdefault(a.resseq, a.resname.strip())
    want_all = DBP + PBP
    refseq = {r: THREE.get(rres[r], "X") for r in want_all if r in rres}

    print("=== what the two pockets are made of ===")
    rows, skipped = [], []
    for pid, path in files:
        s = Structure(path)
        for ch in sorted(s.chains):
            res = {}
            for a in s.protein_atoms:
                if a.chain == ch:
                    res.setdefault(a.resseq, []).append(a)
            same = [r for r in refseq if r in res
                    and THREE.get(res[r][0].resname.strip(), "X") == refseq[r]]
            here = [r for r in refseq if r in res]
            if len(here) < 10 or len(same) < 0.8 * len(here):
                skipped.append(f"{pid}/{ch}")
                continue
            for name, want in (("distal", DBP), ("proximal", PBP)):
                got = [r for r in want if r in res]
                if len(got) < 0.8 * len(want):
                    continue
                seq = [THREE.get(res[r][0].resname.strip(), "X") for r in got]
                sc = [a for r in got for a in sidechain(res[r])]
                if not sc:
                    continue
                apolar = sum(1 for a in sc
                             if (a.element or "").strip().upper() in ("C", "S"))
                rows.append([
                    pid, ch, name, len(got),
                    sum(1 for c in seq if c in AROM),
                    sum(1 for c in seq if c in ALIPH),
                    sum(1 for c in seq if c in POLAR),
                    sum(1 for c in seq if c in CHARGED),
                    sum(1 for c in seq if c in PLUS)
                    - sum(1 for c in seq if c in MINUS),
                    fmt(float(np.mean([KD.get(c, 0.0) for c in seq]))),
                    len(sc), fmt(100.0 * apolar / len(sc), 1),
                    "".join(seq)])

    write_csv(os.path.join(TABLES, "pocket_chemistry.csv"),
              ["pdb", "chain", "pocket", "residues_present", "aromatic",
               "aliphatic", "polar", "charged", "net_charge",
               "mean_kyte_doolittle", "sidechain_atoms",
               "apolar_sidechain_atoms_pct", "sequence"], rows)

    if skipped:
        print(f"    not MexB at these positions, skipped: "
              f"{', '.join(sorted(set(skipped)))}")

    for name in ("distal", "proximal"):
        mine = [r for r in rows if r[2] == name]
        if not mine:
            continue
        lip = np.array([float(r[11]) for r in mine])
        kd = np.array([float(r[9]) for r in mine])
        r0 = mine[0]
        n = r0[3]
        print(f"\n  {name} pocket ({n} residues, {len(mine)} protomers)")
        print(f"    {r0[4]} aromatic, {r0[5]} aliphatic, {r0[6]} polar, "
              f"{r0[7]} charged (net {r0[8]:+d})")
        print(f"    Kyte-Doolittle mean {kd.mean():+.2f} "
              f"(SD {kd.std(ddof=1):.2f} across protomers)")
        print(f"    apolar side-chain atoms {lip.mean():.1f}% "
              f"(SD {lip.std(ddof=1):.1f})")
        print(f"    sequence {r0[12]}")

    # what the ligands actually touch, by site
    env = rows_of(os.path.join(TABLES, "ligand_environment.csv"))
    if env:
        print("\n  contacts the bound ligands actually make (4.5 A):")
        for site, label in (("DBP", "distal only"), ("PBP", "proximal only"),
                            ("both", "spanning both")):
            got = [float(r["percent_apolar"]) for r in env
                   if r.get("site") == site and r.get("percent_apolar")]
            if got:
                print(f"    {label:15} {np.mean(got):5.1f}% apolar "
                      f"(n = {len(got)} ligand copies)")
    print("\nwrote results/tables/pocket_chemistry.csv")


if __name__ == "__main__":
    main()
