#!/usr/bin/env python3
"""Export viewable structures: ChimeraX session scripts and B-factor-painted
PDBs, so the analysis can be inspected in a molecular viewer.

Writes into results/chimerax/:
  <structure>_tunnels.cxc        one-command ChimeraX session: trimer,
                                 pocket residues, ligands, tunnel traces
                                 coloured by local radius
  <structure>_chainE_deviation.pdb
                                 chain E with the per-residue Binding-vs-
                                 Extrusion deviation written into the
                                 B-factor column, so any viewer can colour
                                 by it
"""
from __future__ import annotations

import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import CXDIR, DBP, PBP, SWITCH_LOOP, STRUCT_DIR, TABLES, \
    Structure, load_structures

CHAIN_COLOR = {"D": "#8a8984", "E": "#2a78d6", "F": "#eb6834"}


def spec(resids):
    return ",".join(str(r) for r in sorted(resids))


def write_cxc(s):
    traces = sorted(glob.glob(os.path.join(CXDIR, f"{s.name}_protein_*_tunnel.pdb")))
    lig_names = sorted({n for (_, _, n, _) in s.ligands()})
    out = os.path.join(CXDIR, f"{s.name}_tunnels.cxc")
    L = []
    A = L.append
    A(f"# ChimeraX session for {s.name}")
    A("# run with:  chimerax " + os.path.basename(out))
    A("#   or from inside ChimeraX:  open <this file>")
    A("# Tunnel traces carry the local tunnel radius in the B-factor "
      "column.")
    A("")
    A("set bgColor white")
    A(f"open {os.path.relpath(s.path, CXDIR)}")
    A("hide atoms")
    A("show cartoon")
    A("cartoon style protein modeHelix tube sides 12")
    for ch, col in CHAIN_COLOR.items():
        A(f"color /{ch} {col} cartoons")
    A("transparency /D,F 55 cartoons")
    A("")
    A("# substrate site: distal pocket (red) and proximal pocket (teal)")
    A(f"show /E:{spec(DBP)} atoms")
    A(f"color /E:{spec(DBP)} #e34948 atoms")
    A(f"show /E:{spec(PBP)} atoms")
    A(f"color /E:{spec(PBP)} #1baf7a atoms")
    A(f"show /E:{spec(SWITCH_LOOP)} atoms")
    A(f"color /E:{spec(SWITCH_LOOP)} #eda100 atoms")
    A("style /E sidechain stick")
    A("hide H")
    A("")
    if lig_names:
        A("# bound ligands")
        for n in lig_names:
            A(f"show :{n} atoms")
            A(f"style :{n} ball")
            A(f"color :{n} #4a3aa7 atoms")
        A("")
    A("# tunnel traces, coloured by local radius")
    for i, t in enumerate(traces, start=2):
        A(f"open {os.path.basename(t)}")
    if traces:
        first, last = 2, len(traces) + 1
        rng = f"#{first}" if first == last else f"#{first}-{last}"
        A(f"style {rng} sphere")
        A(f"size {rng} atomRadius 0.35")
        A(f"color byattribute bfactor {rng} "
          f"palette 1.0,#e34948:2.0,#eda100:3.0,#1baf7a:4.5,#2a78d6")
        A(f"# key: red ~1 A (impassable) -> blue ~4.5 A (open)")
    A("")
    A("view")
    A("lighting soft")
    A("graphics silhouettes true width 1.5")
    with open(out, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"  wrote {os.path.relpath(out)}")
    return out


def write_painted(s, chain="E"):
    """Chain `chain` with the E-vs-F sliding deviation in the B column."""
    src = os.path.join(TABLES,
                       f"per_residue_{s.name}_EvsF_TM_trimmed.csv")
    if not os.path.exists(src):
        return None
    with open(src) as fh:
        dev = {int(r["resseq"]): float(r["sliding_rms_A"])
               for r in csv.DictReader(fh)}
    out = os.path.join(CXDIR, f"{s.name}_chain{chain}_deviation.pdb")
    n = 0
    with open(out, "w") as fh:
        fh.write("REMARK  B-factor column = sliding-window Ca deviation "
                 "(A),\n")
        fh.write("REMARK  chain E (Binding) versus chain F (Extrusion), "
                 "fitted on the trimmed TM domain.\n")
        fh.write("REMARK  Colour by B-factor to see where the two states "
                 "differ.\n")
        for line in open(s.path):
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[21] != chain:
                continue
            if line[76:78].strip().upper() == "H":
                continue
            r = int(line[22:26])
            b = dev.get(r, 0.0)
            fh.write(line[:60] + f"{b:6.2f}" + line[66:])
            n += 1
        fh.write("END\n")
    print(f"  wrote {os.path.relpath(out)} ({n} atoms)")
    return out


def main():
    print("=== viewer exports ===")
    for s in load_structures():
        write_cxc(s)
        write_painted(s)


if __name__ == "__main__":
    main()
