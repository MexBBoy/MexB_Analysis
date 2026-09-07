#!/usr/bin/env python3
"""How far into the channel each bound ligand reaches.

Every MexB ligand is projected atom by atom onto one reference channel - the
widest ligand-free route out of ampicillin chain E - after superposing its
protomer on the reference by the pocket-lining CA. That gives each ligand a
span rather than a point: the depth of its shallowest and deepest heavy atom,
which is what "how far in it goes" actually means for a molecule 15-30 A long.

Also reported, per ligand: the free radius left beside it. That is measured in
the structure's own frame, along its own tunnel trace, as the clearance to
protein *and* ligand atoms over the stretch the ligand occupies - the room a
second molecule would have to pass.

Writes results/tables/ligand_reach.csv.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import per_structure_tunnels as pst
from published_pockets import LINING, PDBDIR, load_channel
from mexb_common import (STRUCT_DIR, TABLES, Structure, apply_rt, coords, fmt,
                         kabsch, vdw, write_csv)

NEAR = 6.0          # a trace point this close to the ligand is "beside" it
WINDOW = 20.0       # arc length either side of the ligand its atoms may map to


def depths(Q, RP, rarc, rtot, win, end, u):
    """Depth from the periplasmic mouth for points already in the reference
    frame, continued past the end of the trace.

    The trace stops where its search was seeded, beside the ampicillin
    molecule, which is not the deep end of anything: the distal pocket is a
    chamber and it carries on past that point. Atoms in it would otherwise
    all pile onto the terminal point and read as the same depth. So beyond
    the terminus, depth continues as the distance past it measured along the
    trace's final direction. Nothing before the terminus is affected.
    """
    d = np.linalg.norm(Q[:, None, :] - RP[None, :, :], axis=2)
    d = np.where(win[None, :], d, np.inf)
    j = d.argmin(1)
    dep = rtot - rarc[j]
    beyond = (Q - end) @ u
    # an atom counts as past the end if it lies beyond the plane through the
    # terminus, and its nearest point on the trace is in the deepest 10 A of
    # it - so a molecule sitting 30 A out cannot be dragged past the end by a
    # plane that a winding route crosses more than once
    past = (rarc[j] <= 10.0) & (beyond > 0)
    dep[past] = np.maximum(dep[past], rtot + beyond[past])
    return dep, d[np.arange(len(Q)), j], past


def clearance(points, atoms):
    """min over atoms of (|p - x| - vdw), the same clearance used throughout."""
    X = coords(atoms)
    rad = np.array([vdw(a.element) for a in atoms])
    tree = cKDTree(X)
    out = np.full(len(points), np.inf)
    for k in (16, 48, 128):
        d, j = tree.query(points, k=k)
        out = np.minimum(out, (d - rad[j]).min(axis=1))
        if np.all(d[:, -1] - rad.max() > out):
            break
    return out


def main():
    chan = load_channel()
    if chan is None:
        print("  reference channel missing - run tunnels.py first")
        return
    RP, rarc, rtot = chan
    end = RP[0]                                   # deep terminus of the trace
    u = RP[0] - RP[20]
    u = u / float(np.linalg.norm(u))              # the direction it was going
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    rca = ref.ca("E")

    want = {}
    for r in pst.rows_of(os.path.join(TABLES,
                                      "ligand_environment.csv")):
        want.setdefault((r["pdb"], r["chain"]), []).append(r)

    print("=== how far in each ligand reaches ===")
    rows = []
    for (pid, ch), grp in sorted(want.items()):
        path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            path = os.path.join(PDBDIR, f"{pid}.pdb")
        if not os.path.exists(path):
            continue
        s = Structure(path)
        mca = s.ca(ch)
        common = [r for r in LINING if r in mca and r in rca]
        if len(common) < 25:
            continue
        R, t = kabsch(np.array([mca[r] for r in common]),
                      np.array([rca[r] for r in common]))

        trace = pst.find_trace(s, ch, chan)
        P = pst.read_trace(trace)[0] if trace is not None else None
        prot = [a for a in s.protein_atoms if not a.is_hydrogen]

        for (lch, lres, lname, ats) in s.ligands():
            if lch != ch:
                continue
            heavy = [a for a in ats if not a.is_hydrogen]
            if len(heavy) < 8:                    # ions, glycerol, water-like
                continue
            X = coords(heavy)
            Q = apply_rt(R, t, X)
            # A ligand atom is matched to the stretch of channel around its
            # own molecule, not to the whole path. Depth is arc length along
            # a route that folds back on itself, so an atom a few Angstrom
            # off the centreline can otherwise project onto a coil 30 A away
            # and stretch the span into nonsense.
            cen = X.mean(0)
            QC = apply_rt(R, t, cen[None, :])
            kc = int(np.linalg.norm(QC - RP, axis=1).argmin())
            win = np.abs(rarc - rarc[kc]) <= WINDOW
            dep, off, past = depths(Q, RP, rarc, rtot, win, end, u)
            if float(off.min()) > 12.0:           # not on this channel at all
                continue
            # the marker is the mean of the atom depths rather than the
            # depth of the centroid: past the trace's end the two coordinates
            # are not the same, and only the first is coherent with the bar
            # drawn from the shallowest atom to the deepest
            cdep = float(dep.mean())

            # the free radius left beside it, in the structure's own frame
            free = ""
            if P is not None:
                near = cKDTree(X).query(P)[0] <= NEAR
                if near.sum() >= 3:
                    free = fmt(float(clearance(P[near],
                                               prot + heavy).max()))
            rows.append([pid, ch, f"{lname}{lres}", lname, len(heavy),
                         fmt(cdep), fmt(dep.min()), fmt(dep.max()),
                         fmt(dep.max() - dep.min()), fmt(float(off.mean())),
                         free, int(past.sum()), fmt(rtot)])
            print(f"  {pid:20} {ch} {lname}{lres:>5}: mean depth "
                  f"{cdep:5.1f} A, spans {dep.min():5.1f}"
                  f"–{dep.max():5.1f} A "
                  f"({dep.max() - dep.min():4.1f} A of the path), free "
                  f"radius beside it {free or '-'} A")

    write_csv(os.path.join(TABLES, "ligand_reach.csv"),
              ["pdb", "chain", "ligand", "resname", "heavy_atoms",
               "depth_mean_A", "depth_shallowest_A", "depth_deepest_A",
               "span_A", "mean_offset_from_channel_A",
               "free_radius_beside_A", "atoms_past_the_trace_end",
               "trace_end_depth_A"], rows)
    print("\nwrote results/tables/ligand_reach.csv")


if __name__ == "__main__":
    main()
