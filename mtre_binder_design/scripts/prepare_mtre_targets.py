#!/usr/bin/env python3
"""Build extracellular-only binder-design targets from the MtrE open-state structure.

Downloads PDB 4MT0 (biological assembly 1, the trimer), places it in a membrane
frame, and writes three trimmed target PDBs plus a JSON site definition:

  mtre_trimer_framed.pdb   full trimer, membrane normal on +z, extracellular at +z
  target_crown.pdb         extracellular crown only (loops + rim)
  target_lumen.pdb         crown + beta-barrel lumen wall (the plug site)
  target_loop2.pdb         Loop 2 epitope + immediate rim, one face of the trimer
  mtre_sites.json          residue sets, hotspots, membrane planes

Numbering throughout is 4MT0 AUTH numbering (precursor, includes the 20-residue
signal peptide). Mature numbering = auth - 20. See PROTOCOL_MtrE.md.

Usage:
    python3 scripts/prepare_mtre_targets.py [--outdir targets]
"""

from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select, ShrakeRupley
from Bio.PDB.Polypeptide import index_to_one, three_to_index

ASSEMBLY_URL = "https://files.rcsb.org/download/4MT0-assembly1.cif"

# Ground truth, verified against the coordinates by this script's self-checks.
SIGNAL_PEPTIDE_OFFSET = 20          # auth = mature + 20
ASP_RING = (422, 425)               # auth; mature 402/405. PERIPLASMIC - not a target.
BARREL_STRANDS = ((106, 112), (119, 129), (308, 319), (331, 340))
LOOP1 = (113, 118)                  # core turn, mature 93-98
LOOP2 = (319, 331)                  # 13-aa epitope SVELGGLFKSGTG, mature 299-311

# Membrane slab in the framed coordinate system, set by the two aromatic
# girdles (lower: F137/F305/F306/F344/W346 ~ z+35.8; upper: F326 ~ z+59.5).
MEMBRANE_Z = (35.0, 60.0)
# Everything at or above this z is solvent-accessible from outside the cell.
EXTRACELLULAR_Z = 58.0
# Depth of barrel lumen reachable from outside without passing any constriction.
LUMEN_Z_FLOOR = 34.0

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}


def fetch(outdir: Path) -> Path:
    cif = outdir / "4mt0_assembly1.cif"
    if not cif.exists():
        urllib.request.urlretrieve(ASSEMBLY_URL, cif)
    return cif


def one(res) -> str:
    try:
        return index_to_one(three_to_index(res.get_resname()))
    except Exception:
        return "X"


def load_trimer(cif: Path):
    model = MMCIFParser(QUIET=True).get_structure("mtre", cif)[0]
    for ch in list(model):
        for res in list(ch):
            if res.id[0] != " ":
                ch.detach_child(res.id)
    # The assembly file labels symmetry copies A / A-2 / A-3; PDB format allows
    # one character, and every downstream tool expects A / B / C.
    for src, dst in (("A", "A"), ("A-2", "B"), ("A-3", "C")):
        model[src].id = dst
    chains = [model[c] for c in ("A", "B", "C")]
    return model, chains


def frame(model, chains) -> np.ndarray:
    """Centre on the trimer centroid; 4MT0's crystallographic 3-fold is already
    on z. Flip if needed so that the beta-barrel (extracellular) is at +z."""
    ca = np.array([r["CA"].coord for c in chains for r in c])
    centre = ca.mean(0)
    for atom in model.get_atoms():
        atom.coord = atom.coord - centre

    # The barrel strands must sit at positive z; the Asp ring at negative z.
    barrel_z = np.mean(
        [r["CA"].coord[2] for c in chains for r in c
         if any(lo <= r.id[1] <= hi for lo, hi in BARREL_STRANDS)]
    )
    if barrel_z < 0:
        for atom in model.get_atoms():
            atom.coord = atom.coord * np.array([1.0, 1.0, -1.0])
    return centre


def check(model, chains) -> dict:
    """Self-checks. Any failure means the assumptions above no longer hold and
    every downstream residue list is suspect."""
    res_a = {r.id[1]: r for r in chains[0]}
    problems = []

    for n in ASP_RING:
        if one(res_a[n]) != "D":
            problems.append(f"auth {n} is {one(res_a[n])}, expected D (Asp ring)")

    asp_z = np.mean([res_a[n]["CA"].coord[2] for n in ASP_RING])
    if asp_z > 0:
        problems.append(f"Asp ring at z={asp_z:.1f}; expected periplasmic (negative)")

    barrel_z = np.mean(
        [r["CA"].coord[2] for r in chains[0]
         if any(lo <= r.id[1] <= hi for lo, hi in BARREL_STRANDS)]
    )
    if barrel_z < 0:
        problems.append(f"barrel at z={barrel_z:.1f}; expected extracellular (positive)")

    loop2_seq = "".join(one(res_a[i]) for i in range(LOOP2[0], LOOP2[1] + 1))
    if loop2_seq != "SVELGGLFKSGTG":
        problems.append(f"Loop 2 reads {loop2_seq}, expected SVELGGLFKSGTG")

    if problems:
        raise SystemExit("TARGET PREP FAILED:\n  " + "\n  ".join(problems))

    return {"asp_ring_z": float(asp_z), "barrel_z": float(barrel_z),
            "loop2_seq": loop2_seq}


def pore_profile(chains, z_lo=-60, z_hi=68, step=2.0) -> list[dict]:
    """Radius of the largest sphere that fits on the channel axis at each z."""
    atoms = [a for c in chains for r in c for a in r]
    xyz = np.array([a.coord for a in atoms])
    rad = np.array([VDW.get(a.element, 1.7) for a in atoms])
    owner = [a.get_parent().id[1] for a in atoms]
    out = []
    for z in np.arange(z_lo, z_hi + step, step):
        sel = np.abs(xyz[:, 2] - z) < step
        if sel.sum() < 6:
            continue
        d = np.hypot(xyz[sel, 0], xyz[sel, 1]) - rad[sel]
        idx = np.where(sel)[0]
        lining = sorted({owner[i] for i in idx[np.argsort(d)[:8]]})
        out.append({"z": round(float(z), 1),
                    "pore_radius": round(float(d.min()), 2),
                    "lining": lining})
    return out


def spans(resnums: list[int], present: set[int], max_gap: int = 6) -> list[tuple[int, int]]:
    """Collapse a residue selection into contiguous spans.

    A selection picked per-residue (by SASA, by radius) comes out full of
    one- and two-residue holes. Those holes are not harmless: every one is an
    artificial chain break, which RFdiffusion has to be told about in the
    contig string and which AF2 reads as a real terminus. Filling gaps up to
    max_gap gives contiguous spans that are cleaner to design against and
    simpler to write contigs for, at the cost of a few extra residues."""
    if not resnums:
        return []
    out, start, prev = [], resnums[0], resnums[0]
    for n in resnums[1:]:
        if n - prev <= max_gap + 1:
            prev = n
        else:
            out.append((start, prev))
            start = prev = n
    out.append((start, prev))
    # Never invent residues the structure does not model.
    return [(lo, hi) for lo, hi in out if any(i in present for i in range(lo, hi + 1))]


def expand(resnums: list[int], present: set[int], max_gap: int = 6) -> list[int]:
    keep = []
    for lo, hi in spans(resnums, present, max_gap):
        keep.extend(i for i in range(lo, hi + 1) if i in present)
    return sorted(keep)


def contig(chain_ids: list[str], resnums: list[int], present: set[int]) -> str:
    """RFdiffusion contig fragment for a trimmed target, e.g. 'A105-133/0 A310-340/0'."""
    parts = []
    for ch in chain_ids:
        for lo, hi in spans(resnums, present):
            parts.append(f"{ch}{lo}-{hi}/0")
    return " ".join(parts)


def sasa_by_residue(model) -> dict:
    ShrakeRupley().compute(model, level="R")
    return {(r.get_parent().id, r.id[1]): float(r.sasa)
            for ch in model for r in ch}


def lumen_wall(chains, sasa) -> list[int]:
    """Residues lining the barrel lumen between the extracellular mouth and the
    deepest point reachable without passing a constriction."""
    keep = set()
    for ch in chains:
        for r in ch:
            for a in r:
                z = a.coord[2]
                rr = float(np.hypot(a.coord[0], a.coord[1]))
                if LUMEN_Z_FLOOR <= z <= 68.0 and rr < 14.0:
                    keep.add(r.id[1])
                    break
    return sorted(keep)


def crown(chains, sasa) -> list[int]:
    """Solvent-exposed residues above the outer leaflet."""
    keep = set()
    for ch in chains:
        for r in ch:
            if sasa.get((ch.id, r.id[1]), 0.0) < 10.0:
                continue
            if max(a.coord[2] for a in r) >= EXTRACELLULAR_Z:
                keep.add(r.id[1])
    return sorted(keep)


class Keep(Select):
    def __init__(self, chains, resnums):
        self.chains = set(chains)
        self.resnums = set(resnums)

    def accept_chain(self, chain):
        return chain.id in self.chains

    def accept_residue(self, residue):
        return residue.id[1] in self.resnums

    def accept_atom(self, atom):
        return atom.element != "H"


def write(model, path: Path, chains, resnums):
    io = PDBIO()
    io.set_structure(model)
    io.save(str(path), Keep(chains, resnums))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="targets")
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    cif = fetch(out)
    model, chains = load_trimer(cif)
    frame(model, chains)
    checks = check(model, chains)

    profile = pore_profile(chains)
    sasa = sasa_by_residue(model)
    lumen = lumen_wall(chains, sasa)
    crown_res = crown(chains, sasa)

    all_ids = [c.id for c in chains]
    present = {r.id[1] for r in chains[0]}

    # Collapse each per-residue selection into contiguous spans before writing,
    # so the targets carry no artificial chain breaks (see spans()).
    crown_span = expand(crown_res, present)
    lumen_span = expand(sorted(set(lumen) | set(crown_res)), present)
    loop2_span = expand(sorted(i for i in range(LOOP2[0] - 4, LOOP2[1] + 5)
                               if i in present), present)

    write(model, out / "mtre_trimer_framed.pdb", all_ids, present)
    write(model, out / "target_crown.pdb", all_ids, crown_span)
    write(model, out / "target_lumen.pdb", all_ids, lumen_span)
    # Loop 2 site: one inter-protomer wedge, Loop 2 of A plus the rim of B.
    write(model, out / "target_loop2.pdb", ["A", "B"], loop2_span)

    extracellular_mouth = [p for p in profile if p["z"] >= EXTRACELLULAR_Z]
    sites = {
        "source": "PDB 4MT0 biological assembly 1 (MtrE open state, 3.29 A)",
        "numbering": "4MT0 auth (precursor). mature = auth - %d" % SIGNAL_PEPTIDE_OFFSET,
        "checks": checks,
        "membrane_z": MEMBRANE_Z,
        "extracellular_z": EXTRACELLULAR_Z,
        "lumen_z_floor": LUMEN_Z_FLOOR,
        "asp_ring_auth": list(ASP_RING),
        "asp_ring_note": "PERIPLASMIC. Not reachable from outside. Do not target.",
        "loop1_auth": list(LOOP1),
        "loop2_auth": list(LOOP2),
        "loop2_seq": checks["loop2_seq"],
        "crown_residues_auth": crown_res,
        "lumen_wall_auth": lumen,
        "contigs": {
            "crown": contig(["A", "B", "C"], crown_res, present),
            "lumen": contig(["A", "B", "C"],
                            sorted(set(lumen) | set(crown_res)), present),
            "loop2": contig(["A", "B"],
                            sorted(i for i in range(LOOP2[0] - 4, LOOP2[1] + 5)
                                   if i in present), present),
        },
        "pore_profile": profile,
        "min_extracellular_pore_radius": min(p["pore_radius"] for p in extracellular_mouth),
    }
    (out / "mtre_sites.json").write_text(json.dumps(sites, indent=2))

    print(f"checks passed: {checks}")
    print(f"crown residues ({len(crown_res)}): {crown_res}")
    print(f"lumen wall  ({len(lumen)}): {lumen}")
    print(f"narrowest extracellular pore radius: "
          f"{sites['min_extracellular_pore_radius']:.1f} A")
    print(f"wrote {out}/")


if __name__ == "__main__":
    main()
