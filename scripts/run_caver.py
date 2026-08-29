#!/usr/bin/env python3
"""Cross-check the tunnel numbers with CAVER 3.0.3.

CAVER is the tool reviewers expect for tunnel analysis, so the pipeline's own
bottleneck radii mean much more once a second, independent implementation
agrees with them. This runs CAVER on the same trimers, seeded on the same
points, and writes results/tables/caver.csv alongside our own numbers.

CAVER is not redistributed here. Set CAVER_HOME, or let this script download
caver_3.0.3.zip from caver.cz into work/ on first use.
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import DBP, PBP, TABLES, WORK_DIR, centroid, coords, fmt, \
    load_structures, write_csv

CAVER_URL = "https://www.caver.cz/fil/download/caver30/303/caver_3.0.3.zip"
CAVER_HOME = os.environ.get(
    "CAVER_HOME", os.path.join(WORK_DIR, "caver", "caver_3.0.3", "caver"))

# CAVER defaults from its own example config; probe_radius must sit below the
# bottleneck we are trying to measure or the tunnel is simply not reported.
CONFIG = """load_tunnels no
load_cluster_tree no
time_sparsity 1
first_frame 1
last_frame 1

starting_point_coordinates {x:.3f} {y:.3f} {z:.3f}
starting_point_protection_radius 0
probe_radius {probe}
shell_radius {shell_radius}
shell_depth {shell_depth}

clustering average_link
weighting_coefficient 1
clustering_threshold 3.5

one_tunnel_in_snapshot cheapest
save_dynamics_visualization no
generate_summary yes
generate_tunnel_characteristics yes
generate_tunnel_profiles yes
generate_histograms no
generate_bottleneck_heat_map no
generate_profile_heat_map no

compute_tunnel_residues yes
residue_contact_distance 3.0
compute_bottleneck_residues yes
bottleneck_contact_distance 3.0

number_of_approximating_balls 12
compute_errors no
save_error_profiles no
generate_trajectory no
swap no
seed 1
"""


def ensure_caver():
    jar = os.path.join(CAVER_HOME, "caver.jar")
    if os.path.exists(jar):
        return jar
    zp = os.path.join(WORK_DIR, "caver_3.0.3.zip")
    if not os.path.exists(zp):
        print(f"  downloading CAVER from {CAVER_URL}")
        urllib.request.urlretrieve(CAVER_URL, zp)
    with zipfile.ZipFile(zp) as z:
        z.extractall(os.path.join(WORK_DIR, "caver"))
    if not os.path.exists(jar):
        raise SystemExit(f"CAVER not found at {jar}; set CAVER_HOME")
    return jar


def write_pdb(path, atoms, het_as_atom=False):
    """CAVER 3.0 ignores HETATM records outright - there is no config switch
    for it - so a ligand written as HETATM does not obstruct anything and the
    run silently returns the ligand-free answer. For the 'ligand in place'
    comparison the ligand atoms must be emitted as ATOM records.
    """
    with open(path, "w") as fh:
        for i, a in enumerate(atoms, start=1):
            rec = "ATOM  " if (het_as_atom or not a.hetatm) else "HETATM"
            nm = a.name if len(a.name) >= 4 else f" {a.name:<3.3s}"
            fh.write(f"{rec}{i:5d} {nm:<4.4s} {a.resname:>3.3s} "
                     f"{a.chain}{a.resseq:4d}    "
                     f"{a.xyz[0]:8.3f}{a.xyz[1]:8.3f}{a.xyz[2]:8.3f}"
                     f"  1.00{a.bfac:6.2f}          {a.element:>2.2s}\n")
        fh.write("END\n")


def run_one(jar, tag, atoms, seed, probe, shell_radius, shell_depth,
            het_as_atom=False):
    d = os.path.join(WORK_DIR, "caver_runs", tag)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, "pdbs"), exist_ok=True)
    write_pdb(os.path.join(d, "pdbs", "1.pdb"), atoms,
              het_as_atom=het_as_atom)
    with open(os.path.join(d, "config.txt"), "w") as fh:
        fh.write(CONFIG.format(x=seed[0], y=seed[1], z=seed[2], probe=probe,
                               shell_radius=shell_radius,
                               shell_depth=shell_depth))
    out = os.path.join(d, "out")
    cmd = ["java", "-Xmx6g", "-cp", os.path.join(CAVER_HOME, "lib"),
           "-jar", jar, "-home", CAVER_HOME, "-pdb",
           os.path.join(d, "pdbs"), "-conf", os.path.join(d, "config.txt"),
           "-out", out]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=5400)
    log = os.path.join(out, "log.txt")
    if r.returncode != 0:
        return {"error": (r.stderr or r.stdout or "")[-400:]}
    # CAVER assigns radii from its own atom table and silently discards
    # atoms it cannot place - including most of an arbitrary ligand's atoms,
    # whether written as HETATM or ATOM. Record what it actually loaded so a
    # ligand-in-place comparison is never quoted as if the ligand was there.
    n_in = sum(1 for l in open(os.path.join(d, "pdbs", "1.pdb"))
               if l.startswith(("ATOM", "HETATM")))
    n_loaded = None
    if os.path.exists(log):
        for line in open(log):
            if line.strip().startswith("Atoms:"):
                n_loaded = int(line.split(":")[1].strip())
                break
    tc = os.path.join(out, "analysis", "tunnel_characteristics.csv")
    if not os.path.exists(tc):
        msg = open(log).read()[-400:] if os.path.exists(log) else ""
        return {"error": "no tunnels reported. " + msg}
    rows = []
    with open(tc) as fh:
        for row in csv.DictReader(fh, skipinitialspace=True):
            rows.append(row)
    if not rows:
        return {"error": "no tunnels reported"}
    rows.sort(key=lambda x: -float(x["Bottleneck radius"]))
    best = rows[0]
    # bottleneck-lining residues for the winning cluster
    res = []
    bres = os.path.join(out, "analysis", "residues.txt")
    btl = os.path.join(out, "analysis", "bottlenecks.csv")
    if os.path.exists(btl):
        with open(btl) as fh:
            for row in csv.DictReader(fh, skipinitialspace=True):
                if row.get("Tunnel cluster") == best["Tunnel cluster"]:
                    res = [v for k, v in row.items()
                           if k and "residue" in k.lower() and v]
                    break
    return {"n_atoms_in": n_in, "n_atoms_loaded": n_loaded,
            "n_tunnels": len(rows),
            "cluster": best["Tunnel cluster"],
            "bottleneck": float(best["Bottleneck radius"]),
            "length": float(best["Length"]),
            "curvature": float(best["Curvature"]),
            "throughput": float(best["Throughput"]),
            "residues": ";".join(res)[:200],
            "outdir": out}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", type=float, default=0.9)
    ap.add_argument("--shell-radius", type=float, default=3.0)
    ap.add_argument("--shell-depth", type=float, default=4.0)
    ap.add_argument("--chains", default="E")
    ap.add_argument("--modes", default="protein,withlig")
    a = ap.parse_args()
    jar = ensure_caver()
    print(f"=== CAVER cross-check (jar {jar}) ===")
    print(f"    probe_radius {a.probe}, shell_radius {a.shell_radius}, "
          f"shell_depth {a.shell_depth}")

    ours = {}
    p = os.path.join(TABLES, "tunnels.csv")
    if os.path.exists(p):
        with open(p) as fh:
            for r in csv.DictReader(fh):
                if r["tunnel_rank"] == "1" and not r["note"]:
                    ours.setdefault((r["structure"], r["mode"],
                                     r["chain"]), r)

    rows = []
    for s in load_structures():
        prot = [x for x in s.protein_atoms if not x.is_hydrogen]
        lig = [x for x in s.het_atoms if not x.is_hydrogen]
        ligs_by_chain = {}
        for (lch, lres, lname, ats) in s.ligands():
            h = [x for x in ats if not x.is_hydrogen]
            ligs_by_chain.setdefault(lch, []).append(coords(h).mean(axis=0))
        for ch in a.chains.split(","):
            ca = s.ca(ch)
            if ch in ligs_by_chain:
                seed = ligs_by_chain[ch][0]
                origin = "ligand centroid"
            else:
                seed = 0.5 * (centroid(ca, DBP) + centroid(ca, PBP))
                origin = "transferred DBP/PBP midpoint"
            todo = [m for m in (("protein", prot), ("withlig", prot + lig))
                    if m[0] in a.modes.split(",")]
            for mode, atoms in todo:
                tag = f"{s.name}_{mode}_{ch}"
                print(f"\n  {tag} (seed: {origin})")
                res = run_one(jar, tag, atoms, seed, a.probe,
                              a.shell_radius, a.shell_depth,
                              het_as_atom=(mode == "withlig"))
                mine = ours.get((s.name, mode, ch), {})
                mine_bn = mine.get("bottleneck_radius_A", "")
                if "error" in res:
                    print(f"    CAVER: FAILED - {res['error'][:200]}")
                    rows.append([s.name, mode, ch, origin, a.probe, "", "",
                                 "", "", "", "", "no", mine_bn, "",
                                 res["error"][:200]])
                    continue
                diff = (res["bottleneck"] - float(mine_bn)) if mine_bn else None
                dropped = ((res["n_atoms_in"] - res["n_atoms_loaded"])
                           if res["n_atoms_loaded"] else None)
                print(f"    CAVER   bottleneck {res['bottleneck']:.2f} A, "
                      f"length {res['length']:.1f} A, "
                      f"{res['n_tunnels']} tunnel(s); loaded "
                      f"{res['n_atoms_loaded']}/{res['n_atoms_in']} atoms")
                if mode == "withlig" and dropped and dropped > 10:
                    print(f"    [WARN] CAVER discarded {dropped} atoms - most "
                          f"of the ligand. Its ligand-in-place result is NOT "
                          f"a valid cross-check.")
                print(f"    ours    bottleneck {mine_bn or '—'} A, "
                      f"length {mine.get('geodesic_path_length_A','—')} A")
                if diff is not None:
                    print(f"    difference {diff:+.2f} A")
                if res["residues"]:
                    print(f"    CAVER bottleneck residues: {res['residues']}")
                valid = "yes" if (mode == "protein" or not dropped
                                  or dropped <= 10) else "no - ligand atoms discarded"
                rows.append([s.name, mode, ch, origin, a.probe,
                             fmt(res["bottleneck"]), fmt(res["length"], 1),
                             fmt(res["curvature"]), res["n_tunnels"],
                             res["n_atoms_in"], res["n_atoms_loaded"], valid,
                             mine_bn, fmt(diff) if diff is not None else "",
                             res["residues"]])
    write_csv(os.path.join(TABLES, "caver.csv"),
              ["structure", "mode", "chain", "seed_origin", "probe_radius_A",
               "caver_bottleneck_A", "caver_length_A", "caver_curvature",
               "caver_n_tunnels", "atoms_in_input", "atoms_loaded_by_caver",
               "valid_comparison", "our_bottleneck_A", "difference_A",
               "caver_bottleneck_residues"], rows)
    print("\nwrote results/tables/caver.csv")


if __name__ == "__main__":
    main()
