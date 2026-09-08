#!/usr/bin/env python3
"""The widest route out of every protomer of every MexB structure.

The per-structure tunnel work covered seven protomers, one per ligand
chemistry, because that is what the panel needed. Nothing about the method
was limited to those seven: the search only needs a protomer and a seed. This
runs it over every MexB chain of every structure in hand - ligand-bound and
empty alike - so the tunnels can be compared across the whole set rather than
across a chosen subset.

Method, unchanged from tunnels.py: clearance grid over the protomer's own
trimer, widest-path (max-min) search from the seed to bulk solvent, then the
path re-centred on the medial axis and its local radius measured exactly.

Three things this has to get right that the seven-protomer version could take
for granted:

  the seed     a chain with a pocket ligand is seeded on it; an empty one on
               the transferred DBP/PBP midpoint, which lands inside the
               protein in a closed protomer, so it is nudged to the roomiest
               point within 11 A
  the chain    only chains whose residues match MexB at the pocket positions,
               which drops the MexA and OprM chains of the complexes and the
               MexBYB chimera
  the grid     built once per trimer and reused for its three chains, since
               the grid is most of the cost

Writes results/tables/all_channels.csv, and each trace as a PDB under
results/chimerax/ with the local radius in the B-factor column.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tunnels as T
from per_structure_tunnels import rows_of
from published_pockets import PDBDIR, load_channel, pocket_ligands
from mexb_common import (CXDIR, DBP, PBP, STRUCT_DIR, SUBDOMAINS, TABLES,
                         Structure, centroid, coords, fmt, write_csv)

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
         "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
         "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
         "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
END_TRIM = 3.0      # the route bottleneck ignores the terminal cap at each end
REACH = 11.0        # how far a seed may be nudged to find open space
MOUTH = 12.0        # the stretch of path, at the bulk end, that names the exit

DOCK = set(range(181, 278)) | set(range(718, 814))     # DN and DC


def exit_call(atoms, pts, rad):
    """Name the exit by what lines it, not by which way it points.

    Height along the trimer axis is a poor test: the periplasmic cleft mouth
    sits about 18 A *below* the porter pocket, so a route out of the cleft
    reads as heading down towards the membrane. What distinguishes the exits
    is the subdomain that lines them, which is also how CH1/CH2/CH3 are
    defined in the literature.
    """
    arc = np.concatenate([[0.0], np.cumsum(
        np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    near = arc >= arc[-1] - MOUTH
    P, R = pts[near], rad[near]
    tree = cKDTree(coords(atoms))
    seen = set()
    for q, r in zip(P, R):
        for j in tree.query_ball_point(q, float(r) + 3.0):
            seen.add(atoms[j].resseq)
    if not seen:
        return "unassigned", 0.0, 0.0, 0.0, 0.0
    g = {"PC": 0, "PN": 0, "TM": 0, "dock": 0}
    for r in seen:
        if r in SUBDOMAINS["PC1"] or r in SUBDOMAINS["PC2"]:
            g["PC"] += 1
        elif r in SUBDOMAINS["PN1"] or r in SUBDOMAINS["PN2"]:
            g["PN"] += 1
        elif r in DOCK:
            g["dock"] += 1
        elif r <= 35 or 337 <= r <= 565 or 876 <= r <= 1030:
            g["TM"] += 1
    tot = max(1, sum(g.values()))
    f = {k: v / tot for k, v in g.items()}
    top = max(f, key=lambda k: f[k])
    name = {"PC": "CH1 - PC1/PC2 periplasmic cleft",
            "PN": "CH3 - PN1/PN2 groove",
            "TM": "CH2 - membrane and central cavity",
            "dock": "funnel - docking domain"}[top]
    if f[top] < 0.35:
        name += " (mixed)"
    return name, f["PC"], f["PN"], f["TM"], f["dock"]


def states_of():
    out = {}
    for r in rows_of(os.path.join(TABLES, "protomer_pockets.csv")):
        out[(r["pdb"], r["chain"])] = r["state_call"]
    return out


def mexb_reference():
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    out = {}
    for a in ref.protein_atoms:
        if a.chain == "E" and a.name == "CA":
            out[a.resseq] = THREE.get(a.resname.strip(), "X")
    return {r: out[r] for r in DBP + PBP if r in out}


def is_mexb(s, ch, refseq):
    here = {}
    for a in s.protein_atoms:
        if a.chain == ch and a.name == "CA" and a.resseq in refseq:
            here[a.resseq] = THREE.get(a.resname.strip(), "X")
    if len(here) < 10:
        return False
    same = sum(1 for r in here if here[r] == refseq[r])
    return same >= 0.8 * len(here)


def seed_of(s, ch):
    ca = s.ca(ch)
    dbp = centroid(ca, DBP)
    best, bd, nm = None, np.inf, ""
    for (lch, rn, heavy) in pocket_ligands(s):
        if lch != ch:
            continue
        d = float(np.linalg.norm(coords(heavy).mean(0) - dbp))
        if d < bd:
            best, bd, nm = coords(heavy).mean(0), d, rn
    if best is not None:
        return best, nm, bd
    return 0.5 * (dbp + centroid(ca, PBP)), "", float("nan")


def free_point(clear_fn, p, reach=REACH, step=0.8, floor=1.4):
    g = np.arange(-reach, reach + 1e-9, step)
    off = np.array([[x, y, z] for x in g for y in g for z in g])
    off = off[np.linalg.norm(off, axis=1) <= reach]
    r = clear_fn(p + off)
    k = int(np.argmax(r - 0.02 * np.linalg.norm(off, axis=1)))
    if r[k] < floor:
        return None, 0.0, 0.0
    return p + off[k], float(r[k]), float(np.linalg.norm(off[k]))


def trimer_of(s, ch):
    cen = {c: coords([a for a in s.protein_atoms
                      if a.chain == c and a.name.strip() == "CA"]).mean(0)
           for c in s.chains}
    return tuple(sorted(sorted(cen, key=lambda c: float(
        np.linalg.norm(cen[c] - cen[ch])))[:3]))


def main():
    only = sys.argv[1:] or None
    refseq = mexb_reference()
    state = states_of()
    files = []
    for base in (STRUCT_DIR, PDBDIR):
        if os.path.isdir(base):
            for f in sorted(os.listdir(base)):
                if f.endswith(".pdb") and (only is None or f[:-4] in only):
                    files.append((f[:-4], os.path.join(base, f)))
    seen, todo = set(), []
    for pid, p in files:
        if pid not in seen:
            seen.add(pid)
            todo.append((pid, p))

    print(f"=== the widest route out of every MexB protomer ===")
    rows = []
    for pid, path in todo:
        s = Structure(path)
        chains = [c for c in sorted(s.chains) if is_mexb(s, c, refseq)]
        if not chains:
            print(f"  {pid}: no MexB chain")
            continue
        axis_c, axis = T.trimer_axis(s)
        grids = {}
        for ch in chains:
            t0 = time.time()
            key = trimer_of(s, ch)
            if key not in grids:
                atoms = [a for a in s.protein_atoms
                         if not a.is_hydrogen and a.chain in set(key)]
                grids[key] = (T.ClearanceGrid(atoms, step=T.STEP,
                                              verbose=False),
                              T.ExactClearance(atoms))
            grid, clear_fn = grids[key]
            bulk = T.bulk_region(grid.clearance)

            raw_seed, lig, ligd = seed_of(s, ch)
            seed, clr, moved = free_point(clear_fn, raw_seed)
            if seed is None:
                print(f"  {pid} {ch}: nothing water-wide near the seed")
                continue
            sidx = grid.free_seed(seed)[0]
            if sidx is None:
                print(f"  {pid} {ch}: seed off the grid")
                continue
            R_grid, p = T.widest_path(grid.clearance, sidx, bulk)
            if R_grid is None:
                print(f"  {pid} {ch}: no route to bulk")
                continue
            pts = T.densify(T.refine_path(
                T.smooth_path(np.array([grid.point_of(q) for q in p]), n=1),
                clear_fn), spacing=0.15)
            rad = clear_fn(pts)
            arc = np.concatenate([[0.0], np.cumsum(
                np.linalg.norm(np.diff(pts, axis=0), axis=1))])
            total = float(arc[-1])
            depth = total - arc               # from the bulk end, as always
            mid = (depth >= END_TRIM) & (depth <= total - END_TRIM)
            neck = float(rad[mid].min()) if mid.any() else float(rad.min())
            neck_d = float(depth[int(np.argmin(np.where(mid, rad, np.inf)))]
                           ) if mid.any() else total

            h0 = float(np.dot(pts[0] - axis_c, axis))
            h1 = float(np.dot(pts[-1] - axis_c, axis))
            r0 = float(np.linalg.norm((pts[0] - axis_c) - h0 * axis))
            r1 = float(np.linalg.norm((pts[-1] - axis_c) - h1 * axis))
            dh = h1 - h0
            call, fPC, fPN, fTM, fDK = exit_call(
                [a for a in s.protein_atoms
                 if not a.is_hydrogen and a.chain in set(key)], pts, rad)

            out = os.path.join(CXDIR, f"all_{pid}_{ch}_tunnel.pdb")
            T.write_trace(out, pts, rad)
            rows.append([pid, ch, state.get((pid, ch), ""), lig,
                         fmt(ligd) if lig else "", fmt(moved), fmt(clr),
                         fmt(total), fmt(float(rad[0])), fmt(neck),
                         fmt(neck_d), fmt(float(rad.min())), fmt(dh),
                         fmt(r1 - r0), call, fmt(fPC, 2), fmt(fPN, 2),
                         fmt(fTM, 2), fmt(fDK, 2), os.path.basename(out)])
            print(f"  {pid:20} {ch} [{state.get((pid, ch), '?'):9}] "
                  f"{'on ' + lig if lig else 'empty':>9}: {total:6.1f} A, "
                  f"neck {neck:.2f} A, rises {dh:+6.1f} A -> "
                  f"{call:38} ({time.time() - t0:.0f}s)")

    write_csv(os.path.join(TABLES, "all_channels.csv"),
              ["pdb", "chain", "state", "seed_ligand",
               "ligand_to_distal_pocket_A", "seed_moved_A",
               "seed_clearance_A", "length_A", "radius_at_seed_A",
               "route_bottleneck_A", "route_bottleneck_depth_A",
               "whole_trace_min_A", "rise_along_axis_A", "outward_step_A",
               "exit", "mouth_frac_PC", "mouth_frac_PN", "mouth_frac_TM",
               "mouth_frac_docking", "trace_file"], rows)
    print(f"\nwrote results/tables/all_channels.csv ({len(rows)} protomers)")


if __name__ == "__main__":
    main()
