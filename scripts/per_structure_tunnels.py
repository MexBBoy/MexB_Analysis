#!/usr/bin/env python3
"""Each structure's own tunnel, re-expressed on the reference channel's axis.

P11 draws one reference channel on every row. That keeps the ligand depths
comparable but says nothing about how each structure's own route differs.
This traces the tunnel out of every protomer in the panel, superposes it into
the reference frame, and re-parameterises it by reference depth, so each row
can be drawn at its own measured radius on the shared axis.

The re-parameterisation is the part to be careful with. A protomer's widest
route is not obliged to follow the reference channel - 21FO chain B, for
instance, leaves by a 154 A path where the reference is 63 A. Only trace
points lying within MAX_OFFSET of the reference centreline are kept, and the
fraction of the reference axis they cover is written out, so a row whose
tunnel genuinely goes elsewhere shows up as poor coverage rather than as a
confident-looking profile.

Relaxing MAX_OFFSET does not rescue the coverage - it plateaus at 40-54% by
8 A and does not move out to 20 A, while the median offset of the points that
do project is only 2-4 A. So the gaps are not a cutoff artefact: each
structure's own tunnel simply spans about half the reference axis, because
the projection is many-to-one where a route runs alongside the reference.

A full-length view therefore has to use each tunnel's own arc length, so this
also writes that: each tunnel from its own periplasmic mouth to its own deep
terminus, with its ligands placed on its own axis.

Writes, on the reference axis:
  results/tables/per_structure_tunnels.csv         one row per depth bin
  results/tables/per_structure_tunnel_summary.csv  one row per protomer
and on each structure's own axis:
  results/tables/own_axis_tunnels.csv              full-length radius profile
  results/tables/own_axis_ligands.csv              ligand depths on that axis
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tunnels
from published_pockets import (LINING, PDBDIR, load_channel,
                               pocket_ligands)
from mexb_common import (CXDIR, STRUCT_DIR, TABLES, Structure, apply_rt,
                         coords, fmt, kabsch, write_csv)

STEP = 0.8          # grid spacing for the tunnel search, A
MAX_OFFSET = 6.0    # a trace point further than this from the reference
                    # centreline is not on the same route
BIN = 1.0           # depth bin, A
END_TRIM = 3.0      # the terminal cap of a trace sits inside the seed cavity
                    # and at the bulk mouth; neither is a route constriction
OURS = ("Amp_MexB_20260826", "MexB_DDM_3_20260730")


def rows_of(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh))


def panel_protomers():
    """The (pdb, chain) set P11 draws, chosen the same way P11 chooses it."""
    env = rows_of(os.path.join(TABLES, "ligand_environment.csv"))
    NAME = {"21FP": "Chloramphenicol", "Amp_MexB_20260826": "Ampicillin",
            "2V50": "DDM", "3W9I": "DDM", "21FO": "CYMAL-7",
            "3W9J": "EPI", "6IIA": "LMNG",
            "MexB_DDM_3_20260730": "DDM x3"}
    prot = {}
    for r in env:
        prot.setdefault((r["pdb"], r["chain"]), []).append(r)
    best = {}
    for k, v in prot.items():
        nm = NAME.get(k[0], k[0])
        score = (sum(int(r["residues_contacted"]) for r in v),
                 sum(int(r["heavy_atoms"]) for r in v))
        if nm not in best or score > best[nm][0]:
            best[nm] = (score, k)
    return {k: nm for nm, (_, k) in best.items()}


def read_trace(path):
    pts, rad = [], []
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")):
            pts.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            rad.append(float(ln[60:66]))
    return np.asarray(pts, float), np.asarray(rad, float)


def find_trace(s, chain, chan):
    """Rank-1 ligand-free trace seeded at this protomer's POCKET ligand.

    tunnels.py seeds one search per ligand, so a chain carrying peripheral
    detergent as well as its substrate has several traces. Picking the first
    by filename picks an arbitrary one - for 2V50 chain B that is a surface
    detergent whose route is 10 A long, not the 51 A route from the pocket.
    Match on the pocket ligand instead, and where a protomer holds several
    (the DDM x3 one) take the deepest, whose route spans the whole site.
    """
    pre, suf = f"{s.name}_protein_{chain}_", "_t1_tunnel.pdb"
    cands = [f for f in sorted(os.listdir(CXDIR))
             if f.startswith(pre) and f.endswith(suf)]
    if not cands:
        return None
    pocket = [(rn, h) for (c, rn, h) in pocket_ligands(s) if c == chain]
    if not pocket:                       # apo protomer: the "site" seed
        for f in cands:
            if f[len(pre):-len(suf)] == "site":
                return os.path.join(CXDIR, f)
        return os.path.join(CXDIR, cands[0])

    # deepest pocket ligand first
    RP, rarc, rtot = chan
    scored = []
    for rn, h in pocket:
        cen = coords(h).mean(0)
        scored.append((rtot - rarc[int(np.argmin(
            np.linalg.norm(RP - cen, axis=1)))], cen))
    scored.sort(key=lambda x: -x[0])
    target = scored[0][1]

    best, bd = None, 1e9
    for f in cands:
        P, _ = read_trace(os.path.join(CXDIR, f))
        if not len(P):
            continue
        d = float(np.linalg.norm(P[0] - target))
        if d < bd:
            best, bd = f, d
    return os.path.join(CXDIR, best) if best else None


def write_own_axis(want, chan):
    """Each tunnel over its whole length, on its own arc-length axis.

    Depth is measured from that tunnel's own periplasmic mouth, the same
    convention the reference uses, so the shallow end of every row means the
    same thing. The deep ends do not align, because the tunnels are not the
    same length - that is the honest cost of drawing them in full.
    """
    prof, ligs, summ = [], [], []
    for (pid, ch), nm in sorted(want.items(), key=lambda x: x[1]):
        path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            path = os.path.join(PDBDIR, f"{pid}.pdb")
        if not os.path.exists(path):
            continue
        s = Structure(path)
        trace = find_trace(s, ch, chan)
        if trace is None:
            continue
        P, rad = read_trace(trace)
        if len(P) < 10:
            continue
        arc = np.concatenate([[0.0], np.cumsum(
            np.linalg.norm(np.diff(P, axis=0), axis=1))])
        # the trace runs deep terminus -> bulk exit, so depth = total - arc
        depth = arc[-1] - arc
        for d, r in zip(depth, rad):
            prof.append([pid, ch, nm, fmt(d), fmt(r)])

        # Two different numbers, both worth reporting. The clearance at the
        # deep terminus is how much room the trace has where it was seeded -
        # inside the pocket, beside the ligand. The route bottleneck is the
        # narrowest point the tunnel actually passes through on the way out,
        # taken with the terminal END_TRIM at each end removed so neither the
        # seed cavity nor the bulk mouth can masquerade as a constriction.
        seed = float(rad[0])
        mid = (depth <= depth[0] - END_TRIM) & (depth >= END_TRIM)
        if mid.any():
            k = int(np.argmin(np.where(mid, rad, np.inf)))
            neck, neck_d = float(rad[k]), float(depth[k])
        else:
            neck, neck_d = seed, float(depth[0])
        summ.append([pid, ch, nm, fmt(depth[0]), fmt(seed), fmt(neck),
                     fmt(neck_d), fmt(float(rad.min()))])
        for (c, rn, h) in pocket_ligands(s):
            if c != ch:
                continue
            cen = coords(h).mean(0)
            k = int(np.argmin(np.linalg.norm(P - cen, axis=1)))
            ligs.append([pid, ch, nm, rn, len(h), fmt(depth[k]),
                         fmt(float(np.linalg.norm(P[k] - cen))),
                         fmt(arc[-1])])

    write_csv(os.path.join(TABLES, "own_axis_tunnels.csv"),
              ["pdb", "chain", "ligand", "depth_from_own_mouth_A",
               "radius_A"], prof)
    write_csv(os.path.join(TABLES, "own_axis_ligands.csv"),
              ["pdb", "chain", "ligand", "resname", "heavy_atoms",
               "depth_from_own_mouth_A", "offset_from_own_tunnel_A",
               "tunnel_length_A"], ligs)
    write_csv(os.path.join(TABLES, "own_axis_summary.csv"),
              ["pdb", "chain", "ligand", "tunnel_length_A",
               "seed_clearance_A", "route_bottleneck_A",
               "route_bottleneck_depth_A", "whole_trace_min_A"], summ)
    print("\n  --- full length, on each tunnel's own axis ---")
    print(f"      {'ligand':16} {'length':>7} {'seed':>7} {'neck':>7} "
          f"{'at':>7}")
    for r in summ:
        print(f"      {r[2]:16} {r[3]:>7} {r[4]:>7} {r[5]:>7} {r[6]:>7}")
    print("      seed = clearance where the trace was seeded, beside the "
          "ligand;\n      neck = narrowest point on the route itself "
          f"(terminal {END_TRIM:.0f} A trimmed)")
    print("  wrote own_axis_tunnels.csv, own_axis_ligands.csv and "
          "own_axis_summary.csv")
    return {(r[0], r[1]): r for r in summ}


def main():
    chan = load_channel()
    if chan is None:
        print("  reference channel missing - run tunnels.py first")
        return
    RP, rarc, rtot = chan
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    rca = ref.ca("E")

    want = panel_protomers()
    print("=== each structure's own tunnel, on the reference axis ===")
    print(f"    {len(want)} protomers; trace points kept within "
          f"{MAX_OFFSET:.0f} A of the reference centreline\n")

    bins, summary = [], []
    for (pid, ch), nm in sorted(want.items(), key=lambda x: x[1]):
        path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            path = os.path.join(PDBDIR, f"{pid}.pdb")
        if not os.path.exists(path):
            print(f"  {nm}: {pid} not found - skipped")
            continue
        s = Structure(path)

        trace = find_trace(s, ch, chan)
        if trace is None:
            tunnels.analyse(s, "protein", STEP, max_tunnels=1, verbose=False)
            trace = find_trace(s, ch, chan)
        if trace is None:
            print(f"  {nm}: no trace produced for chain {ch} - skipped")
            continue
        P, rad = read_trace(trace)

        mca = s.ca(ch)
        common = [r for r in LINING if r in mca and r in rca]
        if len(common) < 25:
            print(f"  {nm}: only {len(common)} lining CA in common - skipped")
            continue
        R, t = kabsch(np.array([mca[r] for r in common]),
                      np.array([rca[r] for r in common]))
        Q = apply_rt(R, t, P)

        # project each trace point onto the reference channel
        d = np.linalg.norm(Q[:, None, :] - RP[None, :, :], axis=2)
        j = d.argmin(1)
        off = d[np.arange(len(Q)), j]
        depth = rtot - rarc[j]
        keep = off <= MAX_OFFSET
        if keep.sum() < 20:
            print(f"  {nm:16} {pid} {ch}: only {int(keep.sum())} of {len(Q)} "
                  f"trace points near the reference route - not projectable")
            summary.append([pid, ch, nm, len(P), int(keep.sum()), "", "",
                            fmt(rad.min()), "no"])
            continue

        dd, rr = depth[keep], rad[keep]
        edges = np.arange(0, rtot + BIN, BIN)
        idx = np.digitize(dd, edges) - 1
        got = 0
        for b in range(len(edges) - 1):
            m = idx == b
            if not m.any():
                continue
            got += 1
            bins.append([pid, ch, nm, fmt(edges[b] + BIN / 2),
                         fmt(rr[m].min()), fmt(rr[m].mean()), int(m.sum())])
        cover = 100.0 * got / (len(edges) - 1)
        print(f"  {nm:16} {pid} {ch}: {keep.sum():4d}/{len(Q):4d} points on "
              f"route, covers {cover:4.0f}% of the axis, "
              f"radius {rr.min():.2f}-{rr.max():.2f} A "
              f"(whole tunnel narrowest {rad.min():.2f} A, "
              f"length {len(P) * 0.15:.0f} A)")
        summary.append([pid, ch, nm, len(P), int(keep.sum()), fmt(cover, 0),
                        fmt(rr.min()), fmt(rad.min()),
                        "yes" if cover >= 50 else "partial"])

    own = write_own_axis(want, chan)
    for row in summary:
        o = own.get((row[0], row[1]))
        row.extend([o[4], o[5], o[6]] if o else ["", "", ""])

    write_csv(os.path.join(TABLES, "per_structure_tunnels.csv"),
              ["pdb", "chain", "ligand", "depth_from_entrance_A",
               "radius_min_A", "radius_mean_A", "n_points"], bins)
    write_csv(os.path.join(TABLES, "per_structure_tunnel_summary.csv"),
              ["pdb", "chain", "ligand", "trace_points",
               "points_on_reference_route", "axis_coverage_pct",
               "radius_min_on_route_A", "tunnel_bottleneck_A",
               "projectable", "seed_clearance_A", "route_bottleneck_A",
               "route_bottleneck_depth_A"], summary)
    print("\nwrote results/tables/per_structure_tunnels.csv and "
          "per_structure_tunnel_summary.csv")


if __name__ == "__main__":
    main()
