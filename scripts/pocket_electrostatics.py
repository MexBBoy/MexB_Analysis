#!/usr/bin/env python3
"""Electrostatic potential in the MexB substrate site, per protomer.

After Ramaswamy et al. (Front Microbiol 2018), who characterise the MexB and
MexY pockets by APBS potential as well as by volume and lipophilicity. The
charged residues around the entrance and the two pockets - K134, D174, R620,
R649 among them - are repeatedly named as substrate-binding determinants, and
chloramphenicol is the one ligand in this survey whose contacts are mostly
polar, so a potential map tests whether there is a polar sub-site the
apolar-surface picture misses.

Method. PDB2PQR (AMBER charges, pH 7) then APBS (linearised PB, 0.15 M
implied by the default sdh boundary, 298 K, pdie 2 / sdie 78.54). The site
sphere is defined in the reference frame as everywhere else and mapped back
into each protomer's own frame, so the same anatomical volume is sampled.
Potential is reported in kT/e over the free points of that sphere.

Runs on a subset - our six protomers and one protomer per state per published
structure - because each APBS solve takes ~45 s and several GB.

Writes results/tables/pocket_electrostatics.csv.
"""
from __future__ import annotations

import glob
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from published_pockets import LINING, PDBDIR, ensure_pdbs, numbering_ok
from mexb_analysis import assign_state, cleft_metrics
from mexb_common import (DBP, PBP, STRUCT_DIR, TABLES, WORK_DIR, Structure,
                         apply_rt, centroid, coords, fmt, kabsch, vdw,
                         write_csv)

OURS = ("Amp_MexB_20260826", "MexB_DDM_3_20260730")


def write_pdb(atoms, path):
    with open(path, "w") as fh:
        for i, a in enumerate(atoms, 1):
            nm = a.name if len(a.name) >= 4 else f" {a.name:<3.3s}"
            fh.write(f"ATOM  {i:5d} {nm:<4.4s} {a.resname:>3.3s} "
                     f"A{a.resseq:4d}    {a.xyz[0]:8.3f}{a.xyz[1]:8.3f}"
                     f"{a.xyz[2]:8.3f}  1.00  0.00\n")
        fh.write("END\n")


def read_dx(path):
    """(origin, spacing, values) from an APBS OpenDX potential map."""
    origin = delta = None
    shape, vals = None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("object 1"):
                shape = tuple(int(x) for x in ln.split()[-3:])
            elif ln.startswith("origin"):
                origin = np.array([float(x) for x in ln.split()[1:4]])
            elif ln.startswith("delta") and delta is None:
                delta = float(ln.split()[1])
            elif ln and ln[0].isdigit() or ln.startswith("-"):
                vals.extend(float(x) for x in ln.split())
    return origin, delta, np.asarray(vals).reshape(shape)


def sample(grid, origin, delta, pts):
    """Nearest-node potential at each point; NaN outside the grid."""
    idx = np.rint((pts - origin) / delta).astype(int)
    ok = np.all((idx >= 0) & (idx < np.array(grid.shape)), axis=1)
    out = np.full(len(pts), np.nan)
    ii = idx[ok]
    out[ok] = grid[ii[:, 0], ii[:, 1], ii[:, 2]]
    return out


def main():
    ensure_pdbs()
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    rca = ref.ca("E")
    ref_centre = 0.5 * (centroid(rca, DBP) + centroid(rca, PBP))
    g = np.arange(-16.0, 16.0 + 1e-9, 1.0)
    BALL = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    BALL = BALL[np.linalg.norm(BALL, axis=1) <= 16.0]

    targets = [(os.path.join(STRUCT_DIR, f"{n}.pdb"), n) for n in OURS]
    targets += [(f, os.path.basename(f)[:-4])
                for f in sorted(glob.glob(os.path.join(PDBDIR, "*.pdb")))]

    work = os.path.join(WORK_DIR, "apbs")
    os.makedirs(work, exist_ok=True)
    rows = []
    print("=== electrostatic potential in the substrate site ===\n")
    for path, pid in targets:
        s = Structure(path)
        seen = set()
        for ch in s.chains:
            good, bad = numbering_ok(s, ch)
            if bad or good < 4:
                continue
            m = cleft_metrics(s, ch)
            call, _, _, _ = assign_state(m["PN1-PN2"]["sep"],
                                         m["PC1-PC2"]["sep"])
            # one protomer per state per published structure; all of ours
            if pid not in OURS:
                if call in seen:
                    continue
                seen.add(call)

            mca = s.ca(ch)
            common = [r for r in LINING if r in mca and r in rca]
            if len(common) < 25:
                continue
            R, t = kabsch(np.array([mca[r] for r in common]),
                          np.array([rca[r] for r in common]))
            # the site sphere, mapped back into this protomer's own frame
            local = (BALL + ref_centre - t) @ R

            prot = [a for a in s.protein_atoms
                    if not a.is_hydrogen and a.chain == ch]
            tag = f"{pid}_{ch}"
            pqr = os.path.join(work, f"{tag}.pqr")
            dx = os.path.join(work, f"{tag}.pqr-PE0.dx")
            if not os.path.exists(dx):
                pdb = os.path.join(work, f"{tag}.pdb")
                inp = os.path.join(work, f"{tag}.in")
                write_pdb(prot, pdb)
                r = subprocess.run(["pdb2pqr", "--ff=AMBER", "--whitespace",
                                    f"--apbs-input={inp}", pdb, pqr],
                                   capture_output=True, text=True)
                if not os.path.exists(pqr):
                    print(f"  {tag}: pdb2pqr failed - {r.stderr[-120:]}")
                    continue
                r = subprocess.run(["apbs", os.path.basename(inp)], cwd=work,
                                   capture_output=True, text=True)
                if not os.path.exists(dx):
                    print(f"  {tag}: apbs failed - {r.stdout[-200:]}")
                    continue

            origin, delta, grid = read_dx(dx)
            X = coords(prot)
            tree = cKDTree(X)
            rad = np.array([vdw(a.element) for a in prot])
            dd, ii = tree.query(local, k=1)
            free = (dd - rad[ii]) >= 1.4
            phi = sample(grid, origin, delta, local[free])
            phi = phi[np.isfinite(phi)]
            if len(phi) < 50:
                print(f"  {tag}: only {len(phi)} sampled points - skipped")
                continue
            neg = float((phi < -1).mean() * 100)
            pos = float((phi > 1).mean() * 100)
            print(f"  {pid:22} {ch} {call:9}  mean {phi.mean():+6.2f} kT/e  "
                  f"median {np.median(phi):+6.2f}  "
                  f"{neg:4.1f}% below -1, {pos:4.1f}% above +1  "
                  f"(n={len(phi)})")
            rows.append([pid, ch, call, fmt(phi.mean()), fmt(np.median(phi)),
                         fmt(phi.std()), fmt(neg, 1), fmt(pos, 1), len(phi)])

    write_csv(os.path.join(TABLES, "pocket_electrostatics.csv"),
              ["pdb", "chain", "state_call", "mean_potential_kT_e",
               "median_potential_kT_e", "sd_potential_kT_e",
               "percent_below_minus1", "percent_above_plus1",
               "n_points_sampled"], rows)

    print("\n  --- by state ---")
    for st in ("Access", "Binding", "Extrusion"):
        v = [float(r[3]) for r in rows if r[2] == st]
        if len(v) >= 2:
            print(f"    {st:9} n={len(v):2d}  mean potential "
                  f"{np.mean(v):+.2f} +- {np.std(v, ddof=1):.2f} kT/e")
    print("\nwrote results/tables/pocket_electrostatics.csv")


if __name__ == "__main__":
    main()
