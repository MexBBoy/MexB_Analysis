#!/usr/bin/env python3
"""The whole tunnel through each protomer, entrance to exit, same ends for all.

Earlier traces stopped at the bound ligand, so every row began at the cleft
but ended wherever that particular substrate happened to sit. This runs the
complete path instead, and fixes BOTH ends on anatomy rather than on the
ligand, so the rows are directly comparable end to end:

  entrance  the periplasmic cleft mouth of the reference channel, mapped into
            each structure's frame by superposition on the lining residues
  proximal  the proximal binding pocket, the first site on the way in
  distal    the distal binding pocket, across the switch loop from it
  exit      the funnel on the trimer's three-fold axis, above the docking
            domain, where MexB delivers into TolC

Both pockets are waypoints because that is the path MexB actually takes, and
because one of them alone does not reach every substrate: routed through the
distal pocket only, chloramphenicol sat 13.2 A off the line and dropped out of
the panel altogether. It binds the proximal pocket - 11.8 A from it against
15.4 A from the distal - and so do LMNG and one of the three DDM. Threading
both, in order, is what makes a single line comparable across structures
whichever pocket their ligand occupies.

The funnel is the one exit all three protomers share, which is what makes it
the right common terminus: the axis is taken from the docking-domain centroid
against the whole-trimer centroid, and the exit point sits just beyond the top
of the protein on that axis.

Each leg is the two-stage widest path - max-min for the bottleneck, then the
shortest route through voxels at least that wide, since a widest-path search
has no preference among routes sharing its bottleneck and wanders otherwise.

Every leg is retried rather than abandoned: the mouth and funnel targets grow,
the pocket waypoint falls back from distal to proximal to the ligand itself,
and the seed is nudged to the roomiest nearby free point. What each route
needed is recorded in `route_notes` so a relaxed one is never mistaken for a
clean one.

Writes results/tables/full_tunnels.csv, results/tables/full_tunnel_ligands.csv
and a trace per protomer under results/chimerax/whole_<pdb>_<chain>.pdb
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tunnels as T
import per_structure_tunnels as pst
from cleft_to_ligand import free_point, shortest_within
from published_pockets import LINING, PDBDIR, load_channel, pocket_ligands
from mexb_common import (CXDIR, DBP, PBP, STRUCT_DIR, TABLES, Structure,
                         apply_rt, centroid, coords, fmt, kabsch, write_csv)

DOCK = set(range(181, 278)) | set(range(718, 814))
MOUTH_R = [9.0, 12.0, 16.0, 20.0]     # grown until the target has bulk in it
FUNNEL_R = [10.0, 14.0, 18.0, 24.0]
ABOVE = [6.0, 12.0, 0.0, 18.0]        # how far past the top of the protein


def trimer_axis(s):
    """(centroid, unit axis pointing from the membrane towards the docking domain)."""
    allca, dock = [], []
    for ch in s.chains:
        ca = s.ca(ch)
        allca += list(ca.values())
        dock += [v for r, v in ca.items() if r in DOCK]
    allca, dock = np.array(allca), np.array(dock)
    cen = allca.mean(0)
    ax = dock.mean(0) - cen
    return cen, ax / np.linalg.norm(ax), allca


def ball(grid, bulk, centre, radius):
    gx = [grid.origin[i] + grid.step * np.arange(grid.shape[i])
          for i in range(3)]
    d2 = ((gx[0][:, None, None] - centre[0]) ** 2
          + (gx[1][None, :, None] - centre[1]) ** 2
          + (gx[2][None, None, :] - centre[2]) ** 2)
    return bulk & (d2 <= radius ** 2)


def leg(grid, clear_fn, src, dst):
    """Two-stage widest route between two voxel sets; None if unreachable."""
    Rg, _ = T.widest_path(grid.clearance, src, dst)
    if Rg is None:
        return None, None
    vox = shortest_within(grid, Rg - 1e-6, src, dst)
    if vox is None:
        return None, None
    return np.array([grid.point_of(q) for q in vox]), float(Rg)


def main():
    chan = load_channel()
    RP = chan[0]
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    rca = ref.ca("E")
    want = pst.panel_protomers()

    print("=== the whole tunnel: cleft mouth -> pocket -> funnel ===")
    rows, ligrows = [], []
    for (pid, ch), nm in sorted(want.items(), key=lambda x: x[1]):
        t0 = time.time()
        path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            path = os.path.join(PDBDIR, f"{pid}.pdb")
        if not os.path.exists(path):
            continue
        s = Structure(path)
        mca = s.ca(ch)
        common = [r for r in LINING if r in mca and r in rca]
        if len(common) < 25:
            print(f"  {nm}: too few lining CA in common"); continue
        R, t = kabsch(np.array([rca[r] for r in common]),
                      np.array([mca[r] for r in common]))
        mouth = apply_rt(R, t, RP[-1][None, :])[0]

        cen = {c: coords([a for a in s.protein_atoms if a.chain == c
                          and a.name.strip() == "CA"]).mean(0) for c in s.chains}
        near = sorted(cen, key=lambda c: float(
            np.linalg.norm(cen[c] - cen[ch])))[:3]
        atoms = [a for a in s.protein_atoms
                 if not a.is_hydrogen and a.chain in set(near)]
        grid = T.ClearanceGrid(atoms, step=T.STEP, verbose=False)
        clear_fn = T.ExactClearance(atoms)
        bulk = T.bulk_region(grid.clearance)
        axcen, ax, allca = trimer_axis(s)
        top = float(((allca - axcen) @ ax).max())

        notes = []
        # --- the two fixed ends, each grown until it actually contains bulk
        ent = None
        for r in MOUTH_R:
            m = ball(grid, bulk, mouth, r)
            if m.any():
                ent = m
                if r != MOUTH_R[0]:
                    notes.append(f"mouth target {r:.0f} A")
                break
        exi = None
        for up in ABOVE:
            fp = axcen + ax * (top + up)
            for r in FUNNEL_R:
                m = ball(grid, bulk, fp, r)
                if m.any():
                    exi = m
                    if up != ABOVE[0] or r != FUNNEL_R[0]:
                        notes.append(f"funnel target +{up:.0f} A, {r:.0f} A")
                    break
            if exi is not None:
                break
        if ent is None or exi is None:
            print(f"  {nm}: no {'entrance' if ent is None else 'exit'} target")
            continue

        # --- both pockets as waypoints, in the order the substrate meets them
        ca = s.ca(ch)
        pl = [h for (c, rn, h) in pocket_ligands(s) if c == ch]
        sp = sd = None
        for (lbl, res) in (("proximal", PBP), ("distal", DBP)):
            q, _ = free_point(clear_fn, centroid(ca, res))
            if q is None:
                continue
            if lbl == "proximal":
                sp = grid.free_seed(q)[0]
            else:
                sd = grid.free_seed(q)[0]
        if (sp is None or sd is None) and pl:
            q, _ = free_point(clear_fn, min(
                (coords(h).mean(0) for h in pl),
                key=lambda c: float(np.linalg.norm(c - centroid(ca, DBP)))))
            if q is not None:
                notes.append("waypoint the ligand itself")
                if sp is None:
                    sp = grid.free_seed(q)[0]
                if sd is None:
                    sd = grid.free_seed(q)[0]
        if sp is None or sd is None:
            print(f"  {nm}: no free voxel at a pocket"); continue

        why = "proximal then distal"
        A, ra = leg(grid, clear_fn, sp, ent)            # proximal -> cleft
        tgt = np.zeros(grid.shape, bool)         # the distal seed alone
        tgt[sd] = True
        M, rm = leg(grid, clear_fn, sp, tgt)     # proximal -> distal
        B, rb = leg(grid, clear_fn, sd, exi)            # distal -> funnel
        if A is None or M is None or B is None:
            miss = "cleft" if A is None else ("switch loop" if M is None
                                              else "funnel")
            print(f"  {nm}: no route on the {miss} leg"); continue

        # A runs proximal->cleft, M proximal->distal, B distal->funnel
        pts = np.vstack([A[::-1], M[1:], B[1:]])
        pts = T.densify(T.refine_path(T.smooth_path(pts, n=1), clear_fn),
                        spacing=0.15)
        rad = clear_fn(pts)
        arc = np.concatenate([[0.0], np.cumsum(
            np.linalg.norm(np.diff(pts, axis=0), axis=1))])
        total = float(arc[-1])
        out = os.path.join(CXDIR, f"whole_{pid}_{ch}.pdb")
        T.write_trace(out, pts, rad)
        rows.append([pid, ch, nm, fmt(total), fmt(min(ra, rm, rb)),
                     fmt(float(rad.min())), why, "; ".join(notes) or "clean",
                     os.path.basename(out)])
        print(f"  {nm:16} {pid} {ch}: {total:5.1f} A end to end, "
              f"neck {min(ra, rm, rb):.2f} A  [{'; '.join(notes) or 'clean'}]"
              f"  ({time.time() - t0:.0f}s)")

        # --- where each ligand sits along the whole tunnel
        seen = {}
        for (c, rn, h) in pocket_ligands(s):
            if c != ch:
                continue
            seen[rn] = seen.get(rn, 0) + 1
            tag = f"{rn}{seen[rn]}" if sum(
                1 for (c2, r2, _) in pocket_ligands(s)
                if c2 == ch and r2 == rn) > 1 else rn
            X = coords(h)
            d = np.linalg.norm(X[:, None, :] - pts[None, :, :], axis=2)
            j = d.argmin(1)
            off = d[np.arange(len(X)), j]
            ligrows.append([pid, ch, nm, tag, len(X),
                            fmt(float(arc[j].mean())), fmt(float(arc[j].min())),
                            fmt(float(arc[j].max())), fmt(float(off.mean())),
                            fmt(float(off.min())),
                            fmt(float(np.interp(arc[j].mean(), arc, rad))),
                            fmt(total)])

    write_csv(os.path.join(TABLES, "full_tunnels.csv"),
              ["pdb", "chain", "ligand", "length_A", "leg_bottleneck_A",
               "trace_min_radius_A", "waypoint", "route_notes", "trace_file"],
              rows)
    write_csv(os.path.join(TABLES, "full_tunnel_ligands.csv"),
              ["pdb", "chain", "ligand", "resname", "heavy_atoms",
               "along_tunnel_A", "along_min_A", "along_max_A",
               "mean_offset_A", "closest_offset_A", "radius_there_A",
               "tunnel_length_A"], ligrows)
    print(f"\nwrote results/tables/full_tunnels.csv ({len(rows)} tunnels)")
    print("wrote results/tables/full_tunnel_ligands.csv "
          f"({len(ligrows)} ligands)")


if __name__ == "__main__":
    main()
