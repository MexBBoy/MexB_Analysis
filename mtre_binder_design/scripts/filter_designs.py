#!/usr/bin/env python3
"""Geometry filters for MtrE binder designs, run before any confidence metric.

The false positives you get from co-folding models on an outer-membrane protein
are mostly not subtle. They are designs docked onto surfaces that do not exist
in a cell: the lipid-facing belt of the beta-barrel, or the periplasmic half of
the channel. No amount of ipTM tightening removes them, because the model is
confident and self-consistent - it simply has no bilayer. These are cheap,
deterministic geometric rejections, so apply them first and only spend
AF2/Boltz/Rosetta compute on what survives.

Filters, in order:
  1. MEMBRANE   reject any binder atom below z = +58 (outer-leaflet boundary,
                set by the F326 aromatic girdle at z = +59.5).
  2. PERIPLASM  reject any binder atom below z = +34. Redundant with (1) but
                reported separately because it is the catastrophic failure.
  3. LOS        warn when the binder reaches below z = +62. The LOS inner core
                extends roughly 10-15 A above the lipid headgroups, so contacts
                in the z = 58-62 shell are real protein contacts but may be
                sterically contested on a live cell.
  4. DEPTH      for plug designs only, require that the binder actually inserts:
                at least `--min-insertion` atoms below the mouth plane z = +62.
  5. BURIAL     report buried surface area so you can drop the kissers.

Reads design PDBs where the binder is a single chain not present in the target.

Usage:
    python3 scripts/filter_designs.py designs/mtre_crown/*.pdb
    python3 scripts/filter_designs.py --mode plug --min-insertion 25 \
        designs/mtre_lumen_plug/*.pdb --out accepted_plug.tsv
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, ShrakeRupley

SITES = Path(__file__).resolve().parent.parent / "targets" / "mtre_sites.json"

TARGET_CHAINS = {"A", "B", "C"}


def load_planes() -> dict:
    if SITES.exists():
        s = json.loads(SITES.read_text())
        return {"extracellular": s["extracellular_z"],
                "lumen_floor": s["lumen_z_floor"],
                "membrane": tuple(s["membrane_z"])}
    return {"extracellular": 58.0, "lumen_floor": 34.0, "membrane": (35.0, 60.0)}


def binder_chain(model):
    for ch in model:
        if ch.id not in TARGET_CHAINS:
            return ch
    return None


def evaluate(path: Path, planes: dict, mode: str, min_insertion: int) -> dict | None:
    model = PDBParser(QUIET=True).get_structure(path.stem, str(path))[0]
    ch = binder_chain(model)
    if ch is None:
        return None

    z = np.array([a.coord[2] for r in ch for a in r if a.element != "H"])
    if z.size == 0:
        return None

    reasons = []
    if z.min() < planes["extracellular"]:
        reasons.append(f"membrane_clash(z_min={z.min():.1f})")
    if z.min() < planes["lumen_floor"]:
        reasons.append("periplasmic_reach")

    insertion = int((z < 62.0).sum())
    if mode == "plug" and insertion < min_insertion:
        reasons.append(f"no_insertion({insertion}<{min_insertion})")

    warn = []
    if planes["extracellular"] <= z.min() < 62.0:
        warn.append("LOS_contested")

    sr = ShrakeRupley()
    sr.compute(model, level="C")
    complex_sasa = {c.id: c.sasa for c in model}
    for other in list(model):
        if other.id != ch.id:
            model.detach_child(other.id)
    sr.compute(model, level="C")
    bsa = complex_sasa[ch.id]
    free = model[ch.id].sasa
    buried = free - bsa

    if buried < 600:
        reasons.append(f"low_burial({buried:.0f})")

    return {
        "design": path.stem,
        "n_res": len(list(ch)),
        "z_min": round(float(z.min()), 1),
        "z_max": round(float(z.max()), 1),
        "insertion_atoms": insertion,
        "buried_A2": round(float(buried)),
        "verdict": "PASS" if not reasons else "REJECT",
        "reasons": ";".join(reasons) or "-",
        "warnings": ";".join(warn) or "-",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdbs", nargs="+")
    ap.add_argument("--mode", choices=["crown", "plug", "loop2"], default="crown")
    ap.add_argument("--min-insertion", type=int, default=25,
                    help="plug mode: minimum binder atoms below the mouth plane")
    ap.add_argument("--out", help="write TSV here as well as stdout")
    args = ap.parse_args()

    planes = load_planes()
    rows = []
    for p in args.pdbs:
        try:
            row = evaluate(Path(p), planes, args.mode, args.min_insertion)
        except Exception as exc:                       # noqa: BLE001
            print(f"skipped {p}: {exc}", file=sys.stderr)
            continue
        if row:
            rows.append(row)

    if not rows:
        raise SystemExit("no designs evaluated")

    cols = list(rows[0])
    lines = ["\t".join(cols)] + ["\t".join(str(r[c]) for c in cols) for r in rows]
    print("\n".join(lines))
    if args.out:
        Path(args.out).write_text("\n".join(lines) + "\n")

    n_pass = sum(r["verdict"] == "PASS" for r in rows)
    print(f"\n{n_pass}/{len(rows)} passed geometry "
          f"({100 * n_pass / len(rows):.0f}%)", file=sys.stderr)
    print("Surviving designs go to AF2 initial-guess, then Rosetta, then Boltz-2 "
          "as a third vote, then the decoy panel (TolC 1EK9, OprM 1WP1).",
          file=sys.stderr)


if __name__ == "__main__":
    main()
