#!/usr/bin/env python3
"""Validate the models against their cryo-EM density.

Every conclusion in this pipeline so far rests on coordinates alone. This
stage adds the map, and targets the specific flags the coordinate analysis
could not settle:

  * R971 sits 12-17 A from D407/D408 in all six protomers (known issue 1).
    Is the modelled rotamer supported by density, and is there unexplained
    density near the aspartates that a different rotamer would fit?
  * The DDM occlusion result depends on three detergent molecules being
    real. Are all three supported?
  * The ampicillin contact list depends on one pose, which changed between
    model versions.
  * The tunnel bottleneck was confirmed by two independent tools, but both
    trusted the coordinates. How well resolved are the constriction residues?

Metrics per residue / ligand:
  rscc        real-space correlation between the map and density simulated
              from the model, inside a mask around the group
  mean_z      mean map value at atom centres, in sigma above the map mean
  min_z       weakest single atom, which is what catches a bad side chain

Usage (maps are not in this repository - point at wherever they live):

    python3 map_validation.py --map sharpened.mrc --structure Amp_MexB_20260826 \\
                              --resolution 3.0

Needs only numpy, scipy and mrcfile, so it can be run on the machine that
holds the full-size maps and the resulting CSV brought back here.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from map_tools import Map
from mexb_common import DBP, PBP, RELAY, SWITCH_LOOP, TABLES, coords, fmt, \
    load_structures, write_csv

# electron scattering falls off roughly as a Gaussian at these resolutions;
# width is tied to the map resolution rather than fitted per atom
W = {"C": 6.0, "N": 6.5, "O": 7.0, "S": 12.0, "P": 12.0}


def simulate(atoms, pts, resolution):
    """Density simulated from a set of atoms at the given points."""
    sigma = max(resolution / 2.4, 0.8)
    X = coords(atoms)
    w = np.array([W.get(a.element, 6.0) for a in atoms])
    d2 = ((pts[:, None, :] - X[None, :, :]) ** 2).sum(axis=2)
    return (w[None, :] * np.exp(-d2 / (2 * sigma ** 2))).sum(axis=1)


def group_metrics(m, atoms, resolution, mu, sd, pad=2.2):
    """RSCC and per-atom z for one residue or ligand."""
    heavy = [a for a in atoms if not a.is_hydrogen]
    if not heavy:
        return None
    X = coords(heavy)
    lo = m.world_to_grid(X.min(axis=0) - pad).astype(int)
    hi = m.world_to_grid(X.max(axis=0) + pad).astype(int) + 1
    lo = np.maximum(lo, 0)
    hi = np.minimum(hi, m.shape)
    if np.any(hi <= lo):
        return None
    gi = np.stack(np.meshgrid(*[np.arange(lo[i], hi[i]) for i in range(3)],
                              indexing="ij"), axis=-1).reshape(-1, 3)
    world = m.grid_to_world(gi)
    # restrict to a shell around the atoms so the correlation is not
    # dominated by empty space
    dmin = np.sqrt(((world[:, None, :] - X[None, :, :]) ** 2)
                   .sum(axis=2)).min(axis=1)
    keep = dmin <= pad
    if keep.sum() < 12:
        return None
    world, gi = world[keep], gi[keep]
    obs = m.data[gi[:, 0], gi[:, 1], gi[:, 2]].astype(float)
    calc = simulate(heavy, world, resolution)
    if obs.std() < 1e-9 or calc.std() < 1e-9:
        return None
    rscc = float(np.corrcoef(obs, calc)[0, 1])
    at = (m.sample(X) - mu) / (sd or 1.0)
    return {"rscc": rscc, "mean_z": float(at.mean()),
            "min_z": float(at.min()), "n_atoms": len(heavy),
            "n_voxels": int(keep.sum()),
            "worst_atom": heavy[int(np.argmin(at))].name}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--map", required=True)
    ap.add_argument("--structure", required=True)
    ap.add_argument("--resolution", type=float, default=3.0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    s = next((x for x in load_structures() if x.name == a.structure), None)
    if s is None:
        raise SystemExit(f"structure {a.structure} not found in structures/")
    m = Map(a.map)
    mu, sd = m.stats()
    print(f"=== map validation: {a.structure} vs {os.path.basename(a.map)} ===")
    print(f"  {tuple(m.shape)} voxels, {m.voxel[0]:.3f} A/voxel, "
          f"mean {mu:.4f}, sigma {sd:.4f}")

    # which groups to score, and why
    targets = []
    for ch in s.chains:
        for r in sorted(RELAY):
            targets.append((ch, r, "proton relay (known issue 1/2)"))
        for r in DBP:
            targets.append((ch, r, "distal pocket"))
        for r in PBP:
            targets.append((ch, r, "proximal pocket"))
        for r in SWITCH_LOOP:
            targets.append((ch, r, "switch loop"))
    # constriction residues from the tunnel stage
    p = os.path.join(TABLES, "tunnels.csv")
    if os.path.exists(p):
        with open(p) as fh:
            for row in csv.DictReader(fh):
                if (row["structure"] != a.structure or row["mode"] != "protein"
                        or row["tunnel_rank"] != "1"):
                    continue
                for e in (row["constriction_lining_clearance_A"] or "").split(";"):
                    e = e.strip()
                    if not e:
                        continue
                    import re
                    mm = re.match(r"^[A-Z]{2,3}(\d+)([A-Za-z])", e)
                    if mm:
                        targets.append((mm.group(2), int(mm.group(1)),
                                        "tunnel constriction"))

    seen, rows = set(), []
    for ch, r, why in targets:
        if (ch, r) in seen:
            continue
        seen.add((ch, r))
        at = s.residue_atoms(ch, r)
        if not at:
            continue
        g = group_metrics(m, at, a.resolution, mu, sd)
        if g is None:
            continue
        rows.append([a.structure, ch, r, at[0].resname, why,
                     fmt(g["rscc"], 3), fmt(g["mean_z"], 2),
                     fmt(g["min_z"], 2), g["worst_atom"], g["n_atoms"]])

    for (lch, lres, lname, ats) in s.ligands():
        g = group_metrics(m, ats, a.resolution, mu, sd)
        if g is None:
            continue
        rows.append([a.structure, lch, lres, lname, "bound ligand",
                     fmt(g["rscc"], 3), fmt(g["mean_z"], 2),
                     fmt(g["min_z"], 2), g["worst_atom"], g["n_atoms"],
                     g["n_voxels"]])

    out = a.out or os.path.join(TABLES, f"map_validation_{a.structure}.csv")
    write_csv(out, ["structure", "chain", "resseq", "resname", "why",
                    "rscc", "mean_z", "min_z", "weakest_atom", "n_heavy",
                    "n_mask_voxels"], rows)

    def show(label, sel):
        sub = [r for r in rows if sel(r)]
        if not sub:
            return
        v = [float(r[5]) for r in sub if r[5]]
        print(f"\n  {label}: {len(sub)} groups, median RSCC "
              f"{np.median(v):.3f}")
        for r in sorted(sub, key=lambda x: float(x[5] or 0))[:6]:
            print(f"    {r[3]}{r[2]}{r[1]:<2} RSCC {r[5]:>6} mean_z {r[6]:>6} "
                  f"min_z {r[7]:>6} (weakest {r[8]})  [{r[4]}]")

    show("proton relay", lambda r: "relay" in r[4])
    show("bound ligands", lambda r: r[4] == "bound ligand")
    show("tunnel constriction", lambda r: r[4] == "tunnel constriction")
    show("switch loop", lambda r: r[4] == "switch loop")
    print(f"\n  wrote {out} ({len(rows)} rows)")
    print("\n  Reading these: RSCC below ~0.5 or min_z below ~0.5 sigma "
          "marks a group the density does not support.")
    print("  RSCC is unreliable for very small groups - a glycine with a "
          "handful of atoms in a small mask scores low even against a "
          "perfect map, because there is little density contrast to "
          "correlate. Check n_mask_voxels, and prefer mean_z/min_z for "
          "residues with fewer than ~6 heavy atoms.")


if __name__ == "__main__":
    main()
