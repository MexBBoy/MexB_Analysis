#!/usr/bin/env python3
"""Cut a small, uploadable sub-map around the region that matters.

A full MexB map is far too large to move around; the questions this analysis
needs answered are all local. This writes a cube around a chosen target, which
is typically a few MB instead of several hundred.

Run it where the full maps live:

    python3 crop_map.py --map sharpened.mrc --pdb Amp_MexB_20260826.pdb \\
                        --target ligand --out amp_ligand_crop.mrc

Targets:
  ligand    the bound ligand(s) in the given chain      (pose validation)
  relay     D407 / D408 / K939 / R971 / T976            (known issue 1)
  site      the DBP/PBP midpoint                        (pocket density)
  neck      the tunnel constriction from tunnels.csv    (bottleneck support)
  resi      an explicit residue list via --resi
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from map_tools import Map, crop
from mexb_common import DBP, PBP, TABLES, Structure, centroid, coords


def target_centre(s, chain, kind, resi=None):
    if kind == "ligand":
        at = [a for a in s.het_atoms if a.chain == chain
              and not a.is_hydrogen]
        if not at:
            raise SystemExit(f"no ligand in chain {chain}")
        return coords(at).mean(axis=0), f"{len(at)} ligand atoms"
    if kind == "relay":
        at = []
        for r in (407, 408, 939, 971, 976):
            at += [a for a in s.residue_atoms(chain, r)]
        return coords(at).mean(axis=0), "proton relay residues"
    if kind == "site":
        ca = s.ca(chain)
        return 0.5 * (centroid(ca, DBP) + centroid(ca, PBP)), "DBP/PBP midpoint"
    if kind == "neck":
        p = os.path.join(TABLES, "tunnels.csv")
        with open(p) as fh:
            for r in csv.DictReader(fh):
                if (r["chain"] == chain and r["mode"] == "protein"
                        and r["tunnel_rank"] == "1" and r["narrowest_point_xyz"]):
                    xyz = np.array([float(x) for x in
                                    r["narrowest_point_xyz"].split(",")])
                    return xyz, "tunnel constriction"
        raise SystemExit("no tunnel found in results/tables/tunnels.csv")
    if kind == "resi":
        if not resi:
            raise SystemExit("--resi required for target resi")
        at = []
        for r in [int(x) for x in resi.split(",")]:
            at += s.residue_atoms(chain, r)
        return coords(at).mean(axis=0), f"residues {resi}"
    raise SystemExit(f"unknown target {kind}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--map", required=True)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--chain", default="E")
    ap.add_argument("--target", default="ligand",
                    choices=["ligand", "relay", "site", "neck", "resi"])
    ap.add_argument("--resi")
    ap.add_argument("--radius", type=float, default=15.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    s = Structure(a.pdb)
    m = Map(a.map)
    cen, what = target_centre(s, a.chain, a.target, a.resi)
    print(f"  map {os.path.basename(a.map)}: {tuple(m.shape)} voxels, "
          f"{m.voxel[0]:.3f} A/voxel")
    print(f"  target: {what} at ({cen[0]:.1f}, {cen[1]:.1f}, {cen[2]:.1f})")
    shape, out = crop(m, cen, a.radius, a.out)
    mb = os.path.getsize(out) / 1e6
    print(f"  wrote {out}: {shape} voxels, {mb:.1f} MB")
    if mb > 25:
        print("  (still large - reduce --radius)")


if __name__ == "__main__":
    main()
