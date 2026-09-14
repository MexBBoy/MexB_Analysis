#!/usr/bin/env python3
"""How much room each binding pocket actually encloses, in our two structures.

The tunnel panels reduce each pocket to where it falls along a line, and that
projection loses most of what a pocket is: the centroids sit up to 11 A off
their own trace, and the two pockets separate across the path as much as along
it - 10.4 A apart in every structure, but running along the route by 10.1 A in
ampicillin and 0.0 A in CYMAL-7. A one-dimensional mark cannot express that.
Volume can.

Method. Clearance is measured to van der Waals surfaces on a 0.6 A grid over
the trimer, so a neighbouring protomer walls the pocket properly. A voxel
counts as pocket space if a probe of radius `probe` fits in it and it is NOT
connected to bulk solvent - buried room, not surface. The buried space around
the porter domain is one connected system, so proximal and distal are
separated by assigning each voxel to whichever pocket's lining residues are
nearer, a Voronoi split on the two residue sets.

Two things are reported that a single number would hide. Volumes are given at
three probe radii, because a pocket's volume is a function of the probe and
quoting one figure invites it to be read as an intrinsic constant. And the
ligand is measured against the pocket it sits in, since a pocket that looks
large may be mostly occupied.

Ligands are stripped before the grid is built, so these are the volumes of the
empty site.

Writes results/tables/porter_pocket_volumes.csv
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tunnels as T
from published_pockets import PDBDIR, pocket_ligands
from mexb_common import (DBP, PBP, STRUCT_DIR, TABLES, Structure, coords,
                         centroid, fmt, vdw, write_csv)

PROBES = (1.4, 1.8, 2.2)        # water, then two larger spheres
REACH = 8.0                     # how far from the lining residues counts as pocket
OURS = {("Amp_MexB_20260826", "E"): "Ampicillin",
        ("MexB_DDM_3_20260730", "E"): "DDM x3"}


def main():
    print("=== pocket volumes, ligands stripped ===")
    rows = []
    for (pid, ch), nm in sorted(OURS.items(), key=lambda x: x[1]):
        t0 = time.time()
        path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            path = os.path.join(PDBDIR, f"{pid}.pdb")
        if not os.path.exists(path):
            print(f"  {nm}: no structure"); continue
        s = Structure(path)

        cen = {c: coords([a for a in s.protein_atoms if a.chain == c
                          and a.name.strip() == "CA"]).mean(0) for c in s.chains}
        near = sorted(cen, key=lambda c: float(
            np.linalg.norm(cen[c] - cen[ch])))[:3]
        atoms = [a for a in s.protein_atoms
                 if not a.is_hydrogen and a.chain in set(near)]
        grid = T.ClearanceGrid(atoms, step=T.STEP, verbose=False)
        vox = grid.step ** 3

        ca = s.ca(ch)
        lin = {"proximal": np.array([ca[r] for r in PBP if r in ca]),
               "distal": np.array([ca[r] for r in DBP if r in ca])}
        gx = [grid.origin[i] + grid.step * np.arange(grid.shape[i])
              for i in range(3)]
        G = np.stack(np.meshgrid(*gx, indexing="ij"), -1)

        # nearest lining residue of each set, for the Voronoi split and to
        # keep only voxels that are actually in the pockets' neighbourhood
        dist = {}
        for k, X in lin.items():
            d = np.full(grid.shape, np.inf, np.float32)
            for p in X:
                np.minimum(d, np.linalg.norm(G - p, axis=-1), out=d)
            dist[k] = d
        owner = dist["proximal"] <= dist["distal"]
        nearest = np.minimum(dist["proximal"], dist["distal"])

        lig = [h for (c, rn, h) in pocket_ligands(s) if c == ch]
        for probe in PROBES:
            free = grid.clearance >= probe
            bulk = T.bulk_region(grid.clearance)
            # buried room only: free space not connected to the outside
            buried = free & ~bulk
            lab, n = ndimage.label(buried)
            keep = np.zeros_like(buried)
            for k, X in lin.items():
                c = X.mean(0)
                gi = tuple(int(round((c[i] - grid.origin[i]) / grid.step))
                           for i in range(3))
                sl = tuple(slice(max(0, gi[i] - 25),
                                 min(grid.shape[i], gi[i] + 26))
                           for i in range(3))
                ids = set(np.unique(lab[sl][buried[sl]]))
                for q in ids - {0}:
                    keep |= lab == q
            pocket = keep & (nearest <= REACH)
            for k in ("proximal", "distal"):
                m = pocket & (owner if k == "proximal" else ~owner)
                vol = float(m.sum()) * vox
                occ = np.nan
                if lig:
                    # how much of that space the bound ligand fills
                    inside = np.zeros_like(m)
                    for h in lig:
                        for a in h:
                            p = a.xyz
                            gi = [int(round((p[i] - grid.origin[i])
                                            / grid.step)) for i in range(3)]
                            r = int(np.ceil(vdw(a.element) / grid.step))
                            sl = tuple(slice(max(0, gi[i] - r),
                                             min(grid.shape[i], gi[i] + r + 1))
                                       for i in range(3))
                            d = np.linalg.norm(G[sl] - p, axis=-1)
                            inside[sl] |= d <= vdw(a.element)
                    occ = 100.0 * float((m & inside).sum()) / max(1, m.sum())
                rows.append([pid, ch, nm, k, fmt(probe), fmt(vol),
                             fmt(vol / 1000.0), fmt(occ)])
                print(f"  {nm:10} {k:9} probe {probe:.1f} A: "
                      f"{vol:7.0f} A^3   ligand fills {occ:4.1f}%")
        print(f"    [{nm} {time.time()-t0:.0f}s]")

    write_csv(os.path.join(TABLES, "porter_pocket_volumes.csv"),
              ["pdb", "chain", "structure", "pocket", "probe_A", "volume_A3",
               "volume_nm3", "ligand_fills_pct"], rows)
    print("\nwrote results/tables/porter_pocket_volumes.csv")


if __name__ == "__main__":
    main()
