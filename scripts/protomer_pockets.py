#!/usr/bin/env python3
"""Pocket volume in every protomer of every MexB structure, not just the
ligand-bound one.

Each trimer carries three protomers in different points of the functional
rotation, so each structure contributes three measurements rather than one.
That answers a different question from published_pockets.py: not "does the
pocket enlarge for bigger ligands", but "how much does the pocket change
around the cycle, and do our two structures sit where the published ones do".

Method, and its one important caveat. Every protomer is superposed on the
pocket-lining CA of a single reference (Amp chain E, a Binding protomer) and
the free volume measured with a sphere fixed in that reference frame. So the
sphere sits at the same anatomical position in all of them - which is what
makes the numbers comparable - but for an Access or Extrusion protomer it is
NOT that protomer's own pocket as it would be defined in isolation. The fit
RMSD column says how far each protomer's lining has moved from the reference,
and is the honest measure of how much to trust each row.

Writes results/tables/protomer_pockets.csv.
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.spatial import cKDTree

from mexb_analysis import assign_state, cleft_metrics
from pockets import site_volume, vdw
from published_pockets import (ENTRIES, Frozen, LABEL, LINING, PDBDIR, RADII,
                               ensure_pdbs, numbering_ok, pocket_ligands)
from mexb_common import (DBP, PBP, STRUCT_DIR, TABLES, Structure, apply_rt,
                         centroid, coords, fmt, kabsch, write_csv)


def main():
    ensure_pdbs()
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    rca = ref.ca("E")
    ref_centre = 0.5 * (centroid(rca, DBP) + centroid(rca, PBP))

    targets = [(os.path.join(STRUCT_DIR, f"{n}.pdb"), n)
               for n in ("Amp_MexB_20260826", "MexB_DDM_3_20260730")]
    targets += [(f, os.path.basename(f)[:-4])
                for f in sorted(glob.glob(os.path.join(PDBDIR, "*.pdb")))]

    print("=== pocket volume in every protomer ===")
    print(f"    sphere fixed in the frame of {ref.name} chain E; every "
          f"protomer superposed on its {len(LINING)} pocket-lining CA\n")
    rows = []
    for path, pid in targets:
        s = Structure(path)
        # which protomers hold a ligand in the porter pocket, and how big
        lig_by_chain = {}
        for (ch, rn, heavy) in pocket_ligands(s):
            n, names = lig_by_chain.get(ch, (0, set()))
            lig_by_chain[ch] = (n + len(heavy), names | {rn})

        for ch in s.chains:
            good, bad = numbering_ok(s, ch)
            if bad or good < 4:
                print(f"  {pid} {ch}: numbering check failed - skipped")
                continue
            mca = s.ca(ch)
            common = [r for r in LINING if r in mca and r in rca]
            if len(common) < 25:
                print(f"  {pid} {ch}: only {len(common)} lining CA in "
                      f"common - skipped")
                continue
            M = np.array([mca[r] for r in common])
            T = np.array([rca[r] for r in common])
            R, t = kabsch(M, T)
            fit = float(np.sqrt(((apply_rt(R, t, M) - T) ** 2).sum(1).mean()))

            m = cleft_metrics(s, ch)
            pn, pc = m["PN1-PN2"]["sep"], m["PC1-PC2"]["sep"]
            call, dev, pn_best, pc_best = assign_state(pn, pc)
            agree = "yes" if pn_best == pc_best else "no"

            prot = [a for a in s.protein_atoms
                    if not a.is_hydrogen and a.chain == ch]
            moved = apply_rt(R, t, coords(prot))
            atoms = [Frozen(p, a.element) for p, a in zip(moved, prot)]
            # Two measures, because the first one degenerates. The
            # connected volume needs a free seed at the sphere centre; in an
            # Access or Extrusion protomer that point is often inside the
            # closed lining, and site_volume then falls back to the nearest
            # free voxel and reports whatever isolated bubble it lands in
            # (fractions of an A^3, which round to a plausible-looking 0).
            # Total free volume in the sphere needs no seed and is defined
            # for every protomer, so it is the one that is comparable across
            # the cycle. The connected number is kept where the seed landed
            # in the site rather than in a bubble - the exact midpoint is
            # inside an atom in most protomers, so clearance at the centre
            # is not the test; the size of what the seed found is.
            tree = cKDTree(moved)
            rad = np.array([vdw(a.element) for a in prot])
            dd, ii = tree.query(ref_centre, k=32)
            centre_clear = float((dd - rad[ii]).min())
            vols, tots = {}, {}
            for r in RADII:
                vols[r], tots[r] = site_volume(atoms, ref_centre, radius=r,
                                               step=0.5, probe=1.4)

            # Lipophilic index of the site surface, after the pocket
            # lipophilicity of Ramaswamy et al. (Front Microbiol 2018) but
            # computed directly rather than with the PyMOL MLP plugin: sample
            # the free space in the sphere, take the protein atom nearest each
            # free point that is close enough to line it, and report the
            # fraction of that lining that is carbon or sulfur. No borrowed
            # parameter table, and it is the same quantity the poster's
            # polarity panel shows qualitatively.
            g = np.arange(-16.0, 16.0 + 1e-9, 0.8)
            G = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
            G = G[np.linalg.norm(G, axis=1) <= 16.0] + ref_centre
            dd, ii = tree.query(G, k=1)
            free = dd - rad[ii] >= 1.4          # not inside any atom
            near = dd - rad[ii] <= 5.0          # close enough to be lined by it
            sel = free & near
            if sel.sum() >= 50:
                el = np.array([a.element.strip().upper() for a in prot])[ii[sel]]
                lipo = float(np.isin(el, ("C", "S")).mean() * 100)
            else:
                lipo = float("nan")

            n_heavy, names = lig_by_chain.get(ch, (0, set()))
            ok = (vols[16.0] is not None and tots[16.0] > 0
                  and vols[16.0] >= 0.10 * tots[16.0])
            print(f"  {pid:22} {ch}  {call:9} (PN {pn:5.2f} PC {pc:5.2f}, "
                  f"agree {agree})  lig {n_heavy:>3}  fit {fit:4.2f} A  "
                  f"free "
                  + "  ".join(f"r{int(r)}={tots[r]:.0f}" for r in RADII)
                  + ("" if ok else "   [seed found only an isolated bubble: "
                                    "connected volume not measurable]"))
            rows.append([pid, LABEL.get(pid, ""), ch, call, fmt(dev),
                         pn_best, pc_best, agree, fmt(pn), fmt(pc),
                         "+".join(sorted(names)) or "none", n_heavy,
                         fmt(fit), len(common), fmt(centre_clear),
                         "yes" if ok else "no", fmt(lipo, 1)]
                        + [fmt(tots[r], 0) for r in RADII]
                        + [fmt(vols[r], 0) if ok else "" for r in RADII])

    write_csv(os.path.join(TABLES, "protomer_pockets.csv"),
              ["pdb", "description", "chain", "state_call", "state_L1_dev",
               "PN_nearest", "PC_nearest", "diagnostics_agree",
               "PN1_PN2_sep_A", "PC1_PC2_sep_A", "pocket_ligands",
               "ligand_heavy_atoms", "fit_rmsd_A", "n_lining_CA",
               "centre_clearance_A", "connected_volume_measurable",
               "lipophilic_index_pct"]
              + [f"free_volume_r{int(r)}_A3" for r in RADII]
              + [f"connected_volume_r{int(r)}_A3" for r in RADII], rows)

    print("\n  --- free volume in a 16 A sphere, by state ---")
    for st in ("Access", "Binding", "Extrusion"):
        v = np.array([float(r[18]) for r in rows if r[3] == st])
        if len(v) < 2:
            continue
        print(f"    {st:9} n={len(v):2d}  {v.mean():.0f} +- {v.std(ddof=1):.0f}"
              f"  (range {v.min():.0f}-{v.max():.0f})")
    ours = [r for r in rows if r[0].startswith(("Amp_", "MexB_"))]
    print("\n  --- our two structures against the published spread ---")
    for r in ours:
        same = [float(x[18]) for x in rows if x[3] == r[3]
                and not x[0].startswith(("Amp_", "MexB_"))]
        v = float(r[18])
        if len(same) >= 2:
            z = (v - np.mean(same)) / np.std(same, ddof=1)
            print(f"    {r[0]:22} {r[2]} {r[3]:9} {v:6.0f} A^3  "
                  f"z = {z:+.2f} vs {len(same)} published {r[3]} protomers")
        else:
            print(f"    {r[0]:22} {r[2]} {r[3]:9} {v:6.0f} A^3  "
                  f"(only {len(same)} published {r[3]} protomer to compare)")
    print("\n  --- lipophilic index of the site surface, by state ---")
    for st in ("Access", "Binding", "Extrusion"):
        v = np.array([float(r[16]) for r in rows
                      if r[3] == st and r[16] not in ("", "NA")])
        if len(v) < 2:
            continue
        print(f"    {st:9} n={len(v):2d}  {v.mean():.1f}% +- "
              f"{v.std(ddof=1):.1f}%  (range {v.min():.1f}-{v.max():.1f})")
    mine = [r for r in rows if r[0].startswith(("Amp_", "MexB_"))]
    for r in mine:
        print(f"    {r[0]:22} {r[2]} {r[3]:9} {float(r[16]):.1f}%")
    print("\nwrote results/tables/protomer_pockets.csv")


if __name__ == "__main__":
    main()
