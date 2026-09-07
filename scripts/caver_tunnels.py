#!/usr/bin/env python3
"""Tunnel analysis of every panel structure with CAVER 3.0.3.

Our own widest-path search and CAVER answer slightly different questions:
ours takes the single widest route from the seed to bulk solvent, CAVER
clusters many candidate tunnels and ranks them. Running CAVER on the same
protomers, seeded on the same points, gives an independent full-length
tunnel per structure - coordinates and a radius profile - from the tool
reviewers expect, and lets the two be compared directly.

Ligands are stripped: CAVER assigns radii from its own atom table and
silently discards atoms it cannot place, so a ligand-in-place run is not a
fair comparison (PROTOCOL known issue). The seed is the pocket ligand's
centroid, the same seed our own search uses.

CAVER tunnel coordinates are transformed into the project's reference frame
so they sit on the same axis as everything else, and each tunnel's own
arc length is kept so it can also be drawn at full length.

Writes results/tables/caver_tunnel_profiles.csv and
results/tables/caver_tunnels.csv.
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_caver as rc
from per_structure_tunnels import panel_protomers, rows_of
from published_pockets import (LINING, PDBDIR, load_channel, pocket_ligands)
from mexb_common import (STRUCT_DIR, TABLES, WORK_DIR, Structure, apply_rt,
                         coords, fmt, kabsch, write_csv)

PROBE, SHELL_R, SHELL_D = 0.9, 3.0, 4.0


def best_profile(outdir):
    """(points, radii) of the widest-bottleneck tunnel CAVER reported."""
    f = os.path.join(outdir, "analysis", "tunnel_profiles.csv")
    if not os.path.exists(f):
        return None
    per = {}
    for r in csv.reader(open(f)):
        if len(r) < 14 or not r[1].strip().isdigit():
            continue
        key = (r[1].strip(), r[2].strip())
        axis = r[12].strip()
        if axis in ("X", "Y", "Z", "R"):
            per.setdefault(key, {})[axis] = [
                float(x) for x in r[13:] if x.strip()]
            per[key]["bottleneck"] = float(r[5])
    good = [(v["bottleneck"], k, v) for k, v in per.items()
            if all(a in v for a in "XYZR")]
    if not good:
        return None
    good.sort(key=lambda x: -x[0])
    _, key, v = good[0]
    n = min(len(v["X"]), len(v["Y"]), len(v["Z"]), len(v["R"]))
    P = np.column_stack([v["X"][:n], v["Y"][:n], v["Z"][:n]])
    return P, np.asarray(v["R"][:n]), key


def reuse(outdir):
    """Summarise a CAVER run that already completed, without rerunning it."""
    tc = os.path.join(outdir, "analysis", "tunnel_characteristics.csv")
    rows = list(csv.DictReader(open(tc), skipinitialspace=True))
    if not rows:
        return {"error": "no tunnels in cached run"}
    rows.sort(key=lambda x: -float(x["Bottleneck radius"]))
    b = rows[0]
    return {"n_tunnels": len(rows), "cluster": b["Tunnel cluster"],
            "bottleneck": float(b["Bottleneck radius"]),
            "length": float(b["Length"]), "curvature": float(b["Curvature"]),
            "throughput": float(b["Throughput"]), "residues": "",
            "n_atoms_in": None, "n_atoms_loaded": None, "outdir": outdir}


def flush(prof, summ):
    """Write both tables now, so a killed run keeps what it already has."""
    write_csv(os.path.join(TABLES, "caver_tunnel_profiles.csv"),
              ["pdb", "chain", "ligand", "depth_from_own_mouth_A", "radius_A",
               "depth_on_reference_axis_A", "offset_from_reference_A"], prof)
    write_csv(os.path.join(TABLES, "caver_tunnels.csv"),
              ["pdb", "chain", "ligand", "n_tunnels", "cluster",
               "caver_bottleneck_A", "caver_length_A", "curvature",
               "throughput", "our_bottleneck_A", "atoms_in", "atoms_loaded"],
              summ)


def main():
    jar = rc.ensure_caver()
    if not jar:
        print("  CAVER not available")
        return
    chan = load_channel()
    RP, rarc, rtot = chan
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    rca = ref.ca("E")
    ours = {(r["pdb"], r["chain"]): r
            for r in rows_of(os.path.join(TABLES,
                                          "per_structure_tunnel_summary.csv"))}

    want = panel_protomers()
    print(f"=== CAVER 3.0.3 on {len(want)} protomers "
          f"(probe {PROBE} A, ligands stripped) ===\n")
    prof, summ = [], []
    for (pid, ch), nm in sorted(want.items(), key=lambda x: x[1]):
        path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            path = os.path.join(PDBDIR, f"{pid}.pdb")
        if not os.path.exists(path):
            print(f"  {nm}: {pid} not found"); continue
        s = Structure(path)
        pl = [h for (c, rn, h) in pocket_ligands(s) if c == ch]
        if not pl:
            print(f"  {nm}: no pocket ligand on chain {ch}"); continue
        # deepest pocket ligand, matching per_structure_tunnels
        cens = [coords(h).mean(0) for h in pl]
        depths = [rtot - rarc[int(np.argmin(np.linalg.norm(RP - c, axis=1)))]
                  for c in cens]
        seed = cens[int(np.argmax(depths))]

        # several entries hold two trimers in the asymmetric unit. Only the
        # trimer containing the target chain lines its tunnel, and feeding
        # CAVER the whole file roughly doubles a run that already takes
        # minutes, so keep that chain and its two nearest neighbours.
        cen = {c: coords([a for a in s.protein_atoms
                          if a.chain == c and a.name.strip() == "CA"]).mean(0)
               for c in s.chains}
        near = sorted(cen, key=lambda c: float(np.linalg.norm(
            cen[c] - cen[ch])))[:3]
        atoms = [a for a in s.protein_atoms
                 if not a.is_hydrogen and a.chain in set(near)]
        tag = f"{pid}_{ch}"
        # resume: a completed run leaves tunnel_profiles.csv behind, and each
        # CAVER call is minutes long, so never redo one
        cached = os.path.join(WORK_DIR, "caver_runs", tag, "out")
        if os.path.exists(os.path.join(cached, "analysis",
                                       "tunnel_profiles.csv")):
            res = reuse(cached)
            print(f"  {nm:16} {pid} {ch}: reusing cached run")
        else:
            res = rc.run_one(jar, tag, atoms, seed, PROBE, SHELL_R, SHELL_D)
        if "error" in res:
            print(f"  {nm:16} {pid} {ch}: CAVER failed - {res['error'][:90]}")
            continue
        got = best_profile(res["outdir"])
        if got is None:
            print(f"  {nm:16} {pid} {ch}: no usable profile"); continue
        P, rad, key = got

        mca = s.ca(ch)
        common = [r for r in LINING if r in mca and r in rca]
        R, t = kabsch(np.array([mca[r] for r in common]),
                      np.array([rca[r] for r in common]))
        Q = apply_rt(R, t, P)
        arc = np.concatenate([[0.0], np.cumsum(
            np.linalg.norm(np.diff(P, axis=0), axis=1))])
        # CAVER writes tunnels starting at the seed, so depth from the mouth
        # is total minus arc, the same convention used everywhere here
        own = arc[-1] - arc
        d = np.linalg.norm(Q[:, None, :] - RP[None, :, :], axis=2)
        j = d.argmin(1)
        off = d[np.arange(len(Q)), j]
        refd = rtot - rarc[j]

        for i in range(len(Q)):
            prof.append([pid, ch, nm, fmt(own[i]), fmt(rad[i]),
                         fmt(refd[i]), fmt(off[i])])

        mine = ours.get((pid, ch), {})
        print(f"  {nm:16} {pid} {ch}: {res['n_tunnels']:3d} tunnels, best "
              f"bottleneck {res['bottleneck']:.2f} A over {arc[-1]:.0f} A"
              f"   (ours {mine.get('tunnel_bottleneck_A', '-'):>5} A)")
        summ.append([pid, ch, nm, res["n_tunnels"], key[0],
                     fmt(res["bottleneck"]), fmt(arc[-1]),
                     fmt(res["curvature"]), fmt(res["throughput"]),
                     mine.get("tunnel_bottleneck_A", ""),
                     res["n_atoms_in"], res["n_atoms_loaded"]])
        flush(prof, summ)

    flush(prof, summ)
    print("\nwrote results/tables/caver_tunnels.csv and "
          "caver_tunnel_profiles.csv")


if __name__ == "__main__":
    main()
