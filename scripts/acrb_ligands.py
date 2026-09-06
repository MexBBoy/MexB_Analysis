#!/usr/bin/env python3
"""Where do AcrB's bound drugs sit on the MexB transport path?

P10 drew AcrB's sequential cycle as a schematic because AcrB coordinates were
not measured in this project's frame. They can be. AcrB is 70% identical to
MexB, so its pocket-lining residues map onto MexB's by alignment, and its
protomers superpose on the same reference used everywhere else in this
project. Its bound drugs can then be placed on the same entry channel and
given the same depth coordinate as the MexB ligands.

That turns the comparison into a measurement: whether AcrB's drugs occupy the
same stations, and how many of them a single protomer holds at once.

Detergents and alkanes are excluded - the question is where drugs bind.

Writes results/tables/acrb_ligands.csv.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lining_conservation import align, fetch
from multiligand_survey import ACRB, fetch as fetch_pdb
from published_pockets import LINING, channel_depth, load_channel
from mexb_common import (DBP, PBP, STRUCT_DIR, TABLES, WORK_DIR, Structure,
                         apply_rt, centroid, coords, fmt, kabsch, write_csv)

# detergents, lipids and alkanes seen in the AcrB entries; we want drugs
NOT_DRUG = {"LMT", "LMU", "D10", "D12", "C14", "UND", "DDQ", "PGE", "1PE",
            "P6G", "PG4", "HP6", "OCT", "DAO", "MPD", "GOL", "EDO", "SO4",
            "PO4", "CL", "NA", "MG", "K", "ACT", "HOH", "DMS", "TRS"}
DRUG = {"RFP": "rifampicin", "ERY": "erythromycin", "MIY": "minocycline",
        "DM2": "doxorubicin", "CPF": "ciprofloxacin",
        "5QG": "MBX inhibitor", "5QE": "MBX inhibitor",
        "5QF": "MBX inhibitor", "5QH": "MBX inhibitor"}
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
         "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
         "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
         "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def numbering_ok(s, ch, seq, need=200):
    """AcrB PDB numbering must match P31224 with zero offset."""
    ca = s.ca(ch)
    good = bad = 0
    for r in sorted(ca):
        at = s.residue_atoms(ch, r)
        if not at or not (1 <= r <= len(seq)):
            continue
        aa = THREE.get(at[0].resname)
        if aa is None:
            continue
        (good := good + 1) if seq[r - 1] == aa else (bad := bad + 1)
    return good >= need and bad == 0


def main():
    chan = load_channel()
    if chan is None:
        print("  entry channel trace missing - run tunnels.py first")
        return
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    rca = ref.ca("E")

    mexb, acrb_seq = fetch("P52002"), fetch("P31224")
    m = align(mexb, acrb_seq)          # MexB index -> AcrB index
    # MexB residue n is sequence index n-1 (zero offset, PROTOCOL section 0);
    # AcrB PDB numbering matches P31224 with zero offset, checked per chain
    mapping = {}
    for r in LINING:
        j = m[r - 1] if 0 <= r - 1 < len(m) else -1
        if j >= 0:
            mapping[r] = j + 1
    print("=== AcrB drugs on the MexB transport path ===")
    print(f"    {len(mapping)} of {len(LINING)} pocket-lining residues map "
          f"from MexB onto AcrB by alignment\n")

    acr = os.path.join(WORK_DIR, "acrb")
    rows = []
    for pid in ACRB:
        path = fetch_pdb(pid, acr)
        if not path:
            continue
        s = Structure(path)
        # drugs in this entry, by chain
        bych = {}
        for (ch, rs, rn, ats) in s.ligands():
            h = [a for a in ats if not a.is_hydrogen]
            if rn in NOT_DRUG or len(h) < 10:
                continue
            bych.setdefault(ch, []).append((rn, coords(h)))
        if not bych:
            continue
        for ch, items in sorted(bych.items()):
            if not numbering_ok(s, ch, acrb_seq):
                print(f"  {pid} {ch}: numbering does not match P31224 - "
                      f"skipped")
                continue
            mca = s.ca(ch)
            pairs = [(r, mapping[r]) for r in LINING
                     if r in mapping and r in rca and mapping[r] in mca]
            if len(pairs) < 25:
                print(f"  {pid} {ch}: only {len(pairs)} lining CA in "
                      f"common - skipped")
                continue
            M = np.array([mca[a] for _, a in pairs])
            T = np.array([rca[r] for r, _ in pairs])
            R, t = kabsch(M, T)
            fit = float(np.sqrt(((apply_rt(R, t, M) - T) ** 2).sum(1).mean()))

            for rn, L in items:
                moved = apply_rt(R, t, L)
                depth, off = channel_depth(moved.mean(0), chan)
                if depth is None:
                    continue
                # which MexB lining set does it land in, in the common frame?
                near = []
                for r, a in pairs:
                    d = float(np.linalg.norm(rca[r] - moved.mean(0)))
                    if d < 14:
                        near.append(r)
                n_dbp = len(set(near) & set(DBP))
                n_pbp = len(set(near) & set(PBP))
                site = ("DBP" if n_dbp > n_pbp else
                        "PBP" if n_pbp > n_dbp else
                        "both" if n_dbp else "outside")
                print(f"  {pid} {ch}  {DRUG.get(rn, rn):16} depth "
                      f"{depth:5.1f} A  offset {off:4.1f} A  {site:7} "
                      f"fit {fit:4.2f} A")
                rows.append([pid, ch, rn, DRUG.get(rn, rn), len(L),
                             fmt(depth), fmt(off), site, fmt(fit),
                             len(pairs), len(items)])

    write_csv(os.path.join(TABLES, "acrb_ligands.csv"),
              ["pdb", "chain", "ligand", "name", "heavy_atoms",
               "depth_from_entrance_A", "offset_from_channel_A", "site",
               "fit_rmsd_A", "n_lining_CA_mapped", "drugs_in_this_protomer"],
              rows)

    if rows:
        d = [float(r[5]) for r in rows]
        print(f"\n  {len(rows)} drug copies in {len({r[0] for r in rows})} "
              f"entries, depths {min(d):.1f}-{max(d):.1f} A")
        most = max(int(r[10]) for r in rows)
        print(f"  most drugs in one AcrB protomer: {most}")
    print("\nwrote results/tables/acrb_ligands.csv")


if __name__ == "__main__":
    main()
