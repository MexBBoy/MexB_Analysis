#!/usr/bin/env python3
"""Stage 6 - pocket and cavity analysis.

pyKVFinder for cavity detection (ligand-guided and unguided) and fpocket
where available, plus a grid-based substrate-site volume computed with and
without ligands so that occlusion can be quantified.

The grid volumes here are internally comparable across these structures.
They are NOT drop-in replacements for fpocket or CASTp values and are
labelled as such wherever reported.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import (
    AA3to1, DBP, DETERGENTS, KD, PBP, REPO, TABLES, WORK_DIR, Structure,
    centroid, coords, fmt, load_structures, vdw, write_csv,
)

FPOCKET = os.path.join(WORK_DIR, "fpocket", "bin", "fpocket")
SITE_RADIUS = 16.0     # sphere about the substrate site for grid volumes
GRID_STEP = 0.5
PROBE = 1.4


# ------------------------------------------------------------- pdb writing

def write_pdb(path, atoms):
    with open(path, "w") as fh:
        for i, a in enumerate(atoms, start=1):
            rec = "HETATM" if a.hetatm else "ATOM  "
            # cols 13-16 name, 17 altLoc, 18-20 resname, 22 chain,
            # 23-26 resseq. The altLoc space is required, or every residue
            # name is read one column short (ASN -> "SN").
            nm = a.name if len(a.name) >= 4 else f" {a.name:<3.3s}"
            fh.write(f"{rec}{i:5d} {nm:<4.4s} {a.resname:>3.3s} "
                     f"{a.chain}{a.resseq:4d}    "
                     f"{a.xyz[0]:8.3f}{a.xyz[1]:8.3f}{a.xyz[2]:8.3f}"
                     f"  1.00{a.bfac:6.2f}          {a.element:>2.2s}\n")
        fh.write("END\n")
    return path


# --------------------------------------------------------- grid site volume

def site_volume(atoms, centre, radius=SITE_RADIUS, step=GRID_STEP,
                probe=PROBE):
    """Free volume within `radius` of `centre` that is connected to the
    centre. Free = no atom vdW surface within `probe` of the point.

    Internally comparable only; not equivalent to fpocket/CASTp volumes.
    """
    X = coords(atoms)
    r = np.array([vdw(a.element) for a in atoms])
    keep = np.linalg.norm(X - centre, axis=1) < radius + 12.0
    X, r = X[keep], r[keep]
    n = int(np.ceil(2 * radius / step)) + 1
    ax = centre[0] - radius + step * np.arange(n)
    ay = centre[1] - radius + step * np.arange(n)
    az = centre[2] - radius + step * np.arange(n)
    G = np.stack(np.meshgrid(ax, ay, az, indexing="ij"), axis=-1)
    pts = G.reshape(-1, 3)
    inside = np.linalg.norm(pts - centre, axis=1) <= radius
    tree = cKDTree(X)
    d, idx = tree.query(pts, k=8, workers=-1)
    clear = (d - r[idx]).min(axis=1)
    free = (clear >= probe) & inside
    free = free.reshape(n, n, n)
    if not free.any():
        return 0.0, 0.0
    lab, _ = ndimage.label(free, structure=ndimage.generate_binary_structure(3, 1))
    c = tuple(int(round((centre[i] - (centre[i] - radius)) / step))
              for i in range(3))
    sid = lab[c]
    if sid == 0:
        # centre itself is occluded: take the largest component touching it
        near = lab[max(0, c[0]-6):c[0]+7, max(0, c[1]-6):c[1]+7,
                   max(0, c[2]-6):c[2]+7]
        ids, cnt = np.unique(near[near > 0], return_counts=True)
        if len(ids) == 0:
            return 0.0, float(free.sum() * step ** 3)
        sid = ids[int(np.argmax(cnt))]
    vox = step ** 3
    return float((lab == sid).sum() * vox), float(free.sum() * vox)


# ------------------------------------------------------------- pyKVFinder

def kv_cavities(pdb_path, ligand_path=None, ligand_cutoff=5.0):
    import pyKVFinder
    atomic = pyKVFinder.read_pdb(pdb_path)
    vertices = pyKVFinder.get_vertices(atomic)
    latomic = pyKVFinder.read_pdb(ligand_path) if ligand_path else None
    ncav, cavities = pyKVFinder.detect(
        atomic, vertices, step=0.6, probe_in=1.4, probe_out=4.0,
        removal_distance=2.4, volume_cutoff=5.0,
        latomic=latomic, ligand_cutoff=ligand_cutoff)
    if ncav == 0:
        return None
    surface, volume, area = pyKVFinder.spatial(cavities, step=0.6)
    residues = pyKVFinder.constitutional(cavities, atomic, vertices,
                                         step=0.6)
    depths, max_depth, avg_depth = pyKVFinder.depth(cavities, step=0.6)
    scales, avg_hyd = pyKVFinder.hydropathy(surface, atomic, vertices,
                                            step=0.6)
    return dict(n=ncav, volume=volume, area=area, residues=residues,
                max_depth=max_depth, avg_depth=avg_depth, avg_hyd=avg_hyd)


def res_set(residue_list):
    """pyKVFinder residue entries -> {(chain, resseq)}"""
    out = set()
    for item in residue_list:
        try:
            resseq, chain = int(item[0]), item[1]
        except (ValueError, IndexError):
            continue
        out.add((chain, resseq))
    return out


# ---------------------------------------------------------------- fpocket

def run_fpocket(pdb_path):
    if not os.path.exists(FPOCKET):
        return None
    out = pdb_path[:-4] + "_out"
    shutil.rmtree(out, ignore_errors=True)
    try:
        subprocess.run([FPOCKET, "-f", pdb_path], check=True,
                       capture_output=True, timeout=1800)
    except Exception as e:
        print(f"    fpocket failed: {e}")
        return None
    info = os.path.join(out, os.path.basename(pdb_path)[:-4] + "_info.txt")
    if not os.path.exists(info):
        return None
    pockets, cur = [], None
    for line in open(info):
        line = line.strip()
        if line.startswith("Pocket"):
            cur = {"id": int(line.split()[1])}
            pockets.append(cur)
        elif cur is not None and ":" in line:
            k, v = line.split(":", 1)
            try:
                cur[k.strip()] = float(v.strip())
            except ValueError:
                pass
    return pockets, out


def fpocket_atoms(out_dir, pdb_stem, pid):
    """Residues lining one fpocket pocket."""
    f = os.path.join(out_dir, "pockets", f"pocket{pid}_atm.pdb")
    res = set()
    if not os.path.exists(f):
        return res
    for line in open(f):
        if line.startswith(("ATOM", "HETATM")):
            res.add((line[21], int(line[22:26])))
    return res


# ------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-fpocket", action="store_true")
    a = ap.parse_args()
    structs = load_structures()
    cav_rows, vol_rows, hyd_rows, fp_rows = [], [], [], []
    have_fp = os.path.exists(FPOCKET) and not a.skip_fpocket
    print("=== Stage 6: pockets and cavities ===")
    print(f"  pyKVFinder: available;  fpocket: "
          f"{'available' if have_fp else 'NOT USED'}")

    for s in structs:
        prot = [x for x in s.protein_atoms if not x.is_hydrogen]
        lig_all = [x for x in s.het_atoms if not x.is_hydrogen]
        # Ligands are stripped for the reference cavity/volume calculation so
        # that sites are comparable across structures. For the DDM structure
        # this means removing the three detergent molecules, which otherwise
        # fill the very grooves being measured.
        p_only = os.path.join(WORK_DIR, f"{s.name}_protein.pdb")
        write_pdb(p_only, prot)
        p_with = os.path.join(WORK_DIR, f"{s.name}_withlig.pdb")
        write_pdb(p_with, prot + lig_all)
        print(f"\n  {s.name}")

        # ---- unguided detection on the trimer, ligands stripped
        kv = kv_cavities(p_only)
        site_cav = {}
        if kv:
            print(f"    unguided (ligands stripped): {kv['n']} cavities")
            for ch in s.chains:
                dbp_keys = {(ch, r) for r in DBP}
                pbp_keys = {(ch, r) for r in PBP}
                best, bestov = None, 0
                for cid, reslist in kv["residues"].items():
                    rs = res_set(reslist)
                    ov = len(rs & dbp_keys)
                    if ov > bestov:
                        best, bestov = cid, ov
                if best is None:
                    print(f"      chain {ch}: no cavity overlaps the DBP")
                    continue
                rs = res_set(kv["residues"][best])
                n_pbp = len(rs & pbp_keys)
                cont = ("one continuous cavity" if n_pbp >= 3
                        else "PBP not in the same cavity")
                site_cav[ch] = best
                print(f"      chain {ch}: cavity {best} vol "
                      f"{kv['volume'][best]:.0f} A^3, area "
                      f"{kv['area'][best]:.0f} A^2, max depth "
                      f"{kv['max_depth'][best]:.1f} A, avg hydropathy "
                      f"{kv['avg_hyd'][best]:+.2f}; DBP overlap {bestov}/"
                      f"{len(DBP)}, PBP overlap {n_pbp}/{len(PBP)} "
                      f"-> {cont}")
                cav_rows.append([s.name, ch, "unguided_ligands_stripped",
                                 best, fmt(kv["volume"][best], 0),
                                 fmt(kv["area"][best], 0),
                                 fmt(kv["max_depth"][best], 1),
                                 fmt(kv["avg_depth"][best], 2),
                                 fmt(kv["avg_hyd"][best]), bestov,
                                 n_pbp, cont, len(rs)])

        # ---- ligand-guided detection, per ligand
        for (lch, lres, lname, ats) in s.ligands():
            heavy = [x for x in ats if not x.is_hydrogen]
            lp = os.path.join(WORK_DIR,
                              f"{s.name}_{lname}{lch}{lres}_lig.pdb")
            write_pdb(lp, heavy)
            kvl = kv_cavities(p_only, ligand_path=lp)
            if not kvl:
                print(f"    {lname}{lch}{lres}: ligand-guided detection "
                      f"found no cavity")
                continue
            cid = max(kvl["volume"], key=kvl["volume"].get)
            rs = res_set(kvl["residues"][cid])
            n_dbp = len(rs & {(lch, r) for r in DBP})
            n_pbp = len(rs & {(lch, r) for r in PBP})
            print(f"    {lname}{lch}{lres} ligand-guided: cavity {cid} "
                  f"vol {kvl['volume'][cid]:.0f} A^3, area "
                  f"{kvl['area'][cid]:.0f} A^2, max depth "
                  f"{kvl['max_depth'][cid]:.1f} A, avg hydropathy "
                  f"{kvl['avg_hyd'][cid]:+.2f}, {len(rs)} lining residues "
                  f"(DBP {n_dbp}, PBP {n_pbp})")
            cav_rows.append([s.name, lch, f"ligand_guided_{lname}{lres}",
                             cid, fmt(kvl["volume"][cid], 0),
                             fmt(kvl["area"][cid], 0),
                             fmt(kvl["max_depth"][cid], 1),
                             fmt(kvl["avg_depth"][cid], 2),
                             fmt(kvl["avg_hyd"][cid]), n_dbp, n_pbp,
                             "", len(rs)])

        # ---- grid site volumes, ligands stripped vs present
        for ch in s.chains:
            ca = s.ca(ch)
            cen = 0.5 * (centroid(ca, DBP) + centroid(ca, PBP))
            v_free, _ = site_volume(prot, cen)
            v_occ, _ = site_volume(prot + lig_all, cen)
            occl = v_free - v_occ
            pct = 100.0 * occl / v_free if v_free else float("nan")
            print(f"    chain {ch} site volume (grid, r={SITE_RADIUS} A): "
                  f"stripped {v_free:.0f} A^3, with ligands {v_occ:.0f} "
                  f"A^3, occluded {occl:.0f} A^3 ({pct:.1f}%)")
            vol_rows.append([s.name, ch, fmt(SITE_RADIUS, 0), fmt(GRID_STEP, 2),
                             fmt(v_free, 0), fmt(v_occ, 0), fmt(occl, 0),
                             fmt(pct, 1)])

        # ---- hydropathy of the PBP-lining and DBP-lining subsets
        for ch in s.chains:
            for pname, resids in (("PBP", PBP), ("DBP", DBP)):
                vals = []
                for r in resids:
                    at = s.residue_atoms(ch, r)
                    if at:
                        vals.append(KD.get(AA3to1.get(at[0].resname, "X"), 0))
                hyd_rows.append([s.name, ch, pname, len(vals),
                                 fmt(float(np.mean(vals)))])

        # ---- fpocket on the ligand-stripped trimer
        if have_fp:
            res = run_fpocket(p_only)
            if res:
                pockets, out = res
                stem = os.path.basename(p_only)[:-4]
                ranked = []
                for pk in pockets:
                    lining = fpocket_atoms(out, stem, pk["id"])
                    for ch in s.chains:
                        ov = len(lining & {(ch, r) for r in DBP})
                        if ov >= 3:
                            ranked.append((ov, ch, pk, lining))
                ranked.sort(key=lambda x: -x[0])
                print(f"    fpocket: {len(pockets)} pockets; "
                      f"{len(ranked)} overlap a DBP by >=3 residues")
                for ov, ch, pk, lining in ranked[:6]:
                    n_pbp = len(lining & {(ch, r) for r in PBP})
                    print(f"      chain {ch} pocket {pk['id']}: vol "
                          f"{pk.get('Volume', float('nan')):.0f} A^3, "
                          f"score {pk.get('Score', float('nan')):.3f}, "
                          f"druggability "
                          f"{pk.get('Druggability Score', float('nan')):.3f}"
                          f", DBP {ov}, PBP {n_pbp}")
                    fp_rows.append([
                        s.name, ch, pk["id"],
                        fmt(pk.get("Volume"), 0), fmt(pk.get("Score"), 3),
                        fmt(pk.get("Druggability Score"), 3),
                        fmt(pk.get("Mean local hydrophobic density"), 2),
                        fmt(pk.get("Apolar alpha sphere proportion"), 3),
                        ov, n_pbp, len(lining)])

    write_csv(os.path.join(TABLES, "cavities.csv"),
              ["structure", "chain", "detection", "cavity_id", "volume_A3",
               "area_A2", "max_depth_A", "avg_depth_A", "avg_hydropathy",
               "DBP_overlap", "PBP_overlap", "continuity", "n_lining_res"],
              cav_rows)
    write_csv(os.path.join(TABLES, "pocket_volumes.csv"),
              ["structure", "chain", "site_sphere_radius_A", "grid_step_A",
               "volume_ligands_stripped_A3", "volume_with_ligands_A3",
               "occluded_volume_A3", "occluded_pct"], vol_rows)
    write_csv(os.path.join(TABLES, "pocket_hydropathy.csv"),
              ["structure", "chain", "pocket", "n_res", "mean_KD"], hyd_rows)
    if fp_rows:
        write_csv(os.path.join(TABLES, "fpocket.csv"),
                  ["structure", "chain", "pocket_id", "volume_A3", "score",
                   "druggability", "mean_local_hydrophobic_density",
                   "apolar_alpha_sphere_prop", "DBP_overlap", "PBP_overlap",
                   "n_lining_res"], fp_rows)
    print("\nwrote cavities.csv, pocket_volumes.csv, pocket_hydropathy.csv"
          + (", fpocket.csv" if fp_rows else ""))


if __name__ == "__main__":
    main()
