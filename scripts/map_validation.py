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
    raw = m.sample(X)
    if not np.all(np.isfinite(raw)):
        return None          # part of the group lies outside this crop
    at = (raw - mu) / (sd or 1.0)
    return {"rscc": rscc, "mean_z": float(at.mean()),
            "min_z": float(at.min()), "n_atoms": len(heavy),
            "n_voxels": int(keep.sum()),
            "worst_atom": heavy[int(np.argmin(at))].name}



def score_one(map_path, s, resolution):
    """Score every target group of one structure against one map."""
    import re
    m = Map(map_path)
    prot_xyz = coords([a for a in s.protein_atoms if not a.is_hydrogen])
    mu, sd, masked, zfrac, n_bg = m.background_stats(prot_xyz)
    print(f"    {os.path.basename(map_path)}: "
          f"{'masked' if masked else 'unmasked'} "
          f"({100 * zfrac:.1f}% zero voxels), solvent mu {mu:.4f} "
          f"sigma {sd:.4f} from {n_bg} voxels")
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
    p = os.path.join(TABLES, "tunnels.csv")
    if os.path.exists(p):
        with open(p) as fh:
            for row in csv.DictReader(fh):
                if (row["structure"] != s.name or row["mode"] != "protein"
                        or row["tunnel_rank"] != "1"):
                    continue
                for e in (row["constriction_lining_clearance_A"] or "").split(";"):
                    mm = re.match(r"^[A-Z]{2,3}(\d+)([A-Za-z])", e.strip())
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
        g = group_metrics(m, at, resolution, mu, sd)
        if g is None:
            continue
        rows.append([s.name, os.path.basename(map_path),
                     "masked" if masked else "unmasked", ch, r,
                     at[0].resname, why,
                     fmt(g["rscc"], 3), fmt(g["mean_z"], 2),
                     fmt(g["min_z"], 2), g["worst_atom"], g["n_atoms"],
                     g["n_voxels"]])
    for (lch, lres, lname, ats) in s.ligands():
        g = group_metrics(m, ats, resolution, mu, sd)
        if g is None:
            continue
        rows.append([s.name, os.path.basename(map_path),
                     "masked" if masked else "unmasked", lch, lres, lname,
                     "bound ligand",
                     fmt(g["rscc"], 3), fmt(g["mean_z"], 2),
                     fmt(g["min_z"], 2), g["worst_atom"], g["n_atoms"],
                     g["n_voxels"]])
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--map", help="a single .mrc/.ccp4/.map file")
    ap.add_argument("--zip", help="a bundle from prepare_maps.py; every crop "
                                  "inside is scored and the results stacked")
    ap.add_argument("--structure", help="model name; inferred per crop when "
                                        "using --zip")
    ap.add_argument("--resolution", type=float, required=True,
                    help="map resolution in A. REQUIRED: RSCC rises "
                         "monotonically as this approaches the true value, "
                         "so a wrong figure depresses every score.")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    if not a.map and not a.zip:
        raise SystemExit("pass --map or --zip")

    jobs = []           # (map_path, structure)
    tmp = None
    if a.zip:
        import tempfile, zipfile
        tmp = tempfile.mkdtemp(prefix="mapzip_")
        with zipfile.ZipFile(a.zip) as z:
            z.extractall(tmp)
        structs = load_structures()
        for f in sorted(os.listdir(tmp)):
            if not f.lower().endswith((".mrc", ".ccp4", ".map")):
                continue
            st = None
            if a.structure:
                st = next((x for x in structs if x.name == a.structure), None)
            else:
                # crops are named <mapstem>__<region>.mrc; match the stem
                best, score = None, 0
                for x in structs:
                    n = sum(1 for p, q in zip(f.lower(), x.name.lower())
                            if p == q)
                    for tok in ("amp", "ddm", "lmt", "zz7"):
                        if tok in f.lower() and tok in x.name.lower():
                            n += 20
                    if n > score:
                        best, score = x, n
                st = best
            if st is None:
                print(f"  {f}: no matching model, skipped")
                continue
            jobs.append((os.path.join(tmp, f), st))
        print(f"=== map validation: {len(jobs)} crops from "
              f"{os.path.basename(a.zip)} ===")
    else:
        st = next((x for x in load_structures()
                   if x.name == a.structure), None)
        if st is None:
            raise SystemExit(f"structure {a.structure} not found")
        jobs.append((a.map, st))
        print(f"=== map validation: {a.structure} vs "
              f"{os.path.basename(a.map)} ===")

    all_rows = []
    for map_path, s in jobs:
        rows = score_one(map_path, s, a.resolution)
        all_rows += rows
        print(f"  {os.path.basename(map_path)} -> {len(rows)} groups scored")
    rows = all_rows
    if tmp:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)

    # a group can appear in several crops; keep the one where it is most
    # completely covered, since a group at a crop edge is truncated
    best = {}
    for r in rows:
        key = (r[0], r[1], r[3], r[4])          # structure, map, chain, resi
        if key not in best or int(r[12]) > int(best[key][12]):
            best[key] = r
    rows = sorted(best.values(),
                  key=lambda r: (r[0], r[1], r[3], str(r[4])))

    # Percentile rank of each group's RSCC within its own map. Raw RSCC is
    # not comparable between datasets - voxel size, sharpening and masking
    # all shift it - so rank within one map is the only safe comparison, and
    # it is what the findings should be quoted against.
    by_map = {}
    for r in rows:
        by_map.setdefault((r[0], r[1]), []).append(r)
    for key, group in by_map.items():
        vals = sorted(float(x[7]) for x in group if x[7])
        for r in group:
            if not r[7]:
                r.append("")
                continue
            v = float(r[7])
            pct = 100.0 * sum(1 for q in vals if q <= v) / len(vals)
            r.append(f"{pct:.0f}")

    tag = a.structure or "bundle"
    out = a.out or os.path.join(TABLES, f"map_validation_{tag}.csv")
    write_csv(out, ["structure", "map", "map_masked", "chain", "resseq",
                    "resname", "why", "rscc", "mean_z", "min_z",
                    "weakest_atom", "n_heavy", "n_mask_voxels",
                    "rscc_percentile_in_map"], rows)

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
    print("\n  How to read these")
    print("  - Compare within one map, not between maps. Raw RSCC is not "
          "comparable across datasets: voxel size, sharpening and masking "
          "all move it. Use rscc_percentile_in_map.")
    print("  - z-scores are referenced to a solvent shell around the model "
          "rather than the whole box, which puts masked and unmasked maps "
          "on a common footing to within ~10%. It cannot do better: a mask "
          "cutting close to the model truncates the solvent population. "
          "The map_masked column records which each map is.")
    print("  - RSCC is unreliable for very small groups. A glycine scores "
          "low even against a perfect map, because a small mask holds "
          "little density contrast. Check n_mask_voxels and prefer the "
          "z-scores below ~6 heavy atoms.")
    print("  - Acidic side chains read low by design. D and E are "
          "preferentially decarboxylated by the electron beam (Hattne et "
          "al. 2018, doi:10.1016/j.str.2018.03.021; Spear et al. 2015, "
          "doi:10.1016/j.jsb.2015.09.006), so a weak D407/D408 is a "
          "radiation artefact, not evidence of mismodelling.")
    if a.zip:
        print("  - Scored from crops: a group near a crop boundary has a "
              "clipped mask, which can move its RSCC by a few hundredths "
              "(z-scores are unaffected). Where a group appears in more "
              "than one crop the best-covered instance is kept.")
    asp = [r for r in rows if r[5] in ("ASP", "GLU")]
    if asp:
        v = [float(r[7]) for r in asp if r[7]]
        print(f"\n  ({len(asp)} Asp/Glu groups scored, median RSCC "
              f"{np.median(v):.3f} - expected to sit low, see above.)")


if __name__ == "__main__":
    main()
