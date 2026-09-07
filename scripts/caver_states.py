#!/usr/bin/env python3
"""CAVER tunnels for all three protomers of one trimer, by functional state.

A single protomer's tunnels answer "where can something get in or out of this
one", which is not the mechanism. The mechanism is the three states side by
side: the access protomer taking substrate in from the cleft, the binding
protomer holding it in the pocket, the extrusion protomer passing it up the
funnel. This runs CAVER on each chain of one structure, seeded in that chain's
own pocket, and labels every tunnel by where it comes out - measured on the
trimer's pseudo-3-fold axis, not by eye.

MexB_DDM_3_20260730 is the useful one: its three chains call Access, Binding
and Extrusion, so one structure carries the whole cycle.

Writes results/tables/caver_states.csv and, via caver_scene, a combined scene.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_caver as rc
import tunnels as T
from caver_scene import cluster_spheres
from caver_tunnels import PROBE, SHELL_D, SHELL_R, reuse
from per_structure_tunnels import rows_of
from published_pockets import PDBDIR, load_channel, pocket_ligands
from mexb_common import (DBP, PBP, STRUCT_DIR, TABLES, WORK_DIR, Structure,
                         centroid, coords, fmt, write_csv)

DEFAULT = "MexB_DDM_3_20260730"
OUT = os.path.join(os.path.dirname(TABLES), "structures")
# the poster's protomer colours
SCOL = {"Access": "#0A9DA0", "Binding": "#CA0FC1", "Extrusion": "#0F9C1B"}
UP = 8.0            # Angstrom along the axis before an exit counts as "up"
N_CLUSTERS = 6      # how many clusters to characterise per protomer


def states_of(pid):
    """{chain: Access|Binding|Extrusion} from the protomer table."""
    out = {}
    for r in rows_of(os.path.join(TABLES, "protomer_pockets.csv")):
        if r["pdb"] == pid:
            out[r["chain"]] = r["state_call"]
    return out


def seed_of(s, ch):
    """That chain's own pocket.

    The ligand nearest the distal pocket if the chain carries one - by
    distance to the DBP centroid, which is measured in the structure's own
    frame and so cannot be thrown by a peripheral detergent the way a depth
    on the reference channel would be - otherwise the transferred DBP/PBP
    midpoint, as an empty protomer has nothing to seed on.
    """
    ca = s.ca(ch)
    dbp = centroid(ca, DBP)
    best, bd, nm = None, np.inf, ""
    for (lch, rn, heavy) in [(c, rn, h) for (c, rn, h) in pocket_ligands(s)]:
        if lch != ch:
            continue
        c = coords(heavy).mean(0)
        d = float(np.linalg.norm(c - dbp))
        if d < bd:
            best, bd, nm = c, d, rn
    if best is not None:
        return best, f"{nm} {bd:.0f} A from the distal pocket"
    return 0.5 * (dbp + centroid(ca, PBP)), "DBP/PBP midpoint (no ligand)"


def free_point(cf, p, reach=11.0, step=0.8, floor=1.4):
    """Nudge a seed off the protein.

    The transferred DBP/PBP midpoint lands inside the protein in a closed
    protomer - clearance -0.45 A in the access chain here - and CAVER then
    reports no tunnels at all rather than complaining. So take the roomiest
    point within `reach`, mildly preferring the nearest, and refuse to seed
    if nothing there is even water-wide.
    """
    g = np.arange(-reach, reach + 1e-9, step)
    off = np.array([[x, y, z] for x in g for y in g for z in g])
    off = off[np.linalg.norm(off, axis=1) <= reach]
    r = cf(p + off)
    k = int(np.argmax(r - 0.02 * np.linalg.norm(off, axis=1)))
    if r[k] < floor:
        return None, 0.0, 0.0
    return p + off[k], float(r[k]), float(np.linalg.norm(off[k]))


def exit_call(dh, dr):
    """Where a tunnel comes out, from its rise and its outward step."""
    if dh > UP:
        return "up - funnel and docking domain"
    if dh < -UP:
        return "down - towards the membrane"
    return "out - periplasmic cleft"


def pick(rows, ch):
    """The tunnels worth drawing for one protomer.

    The best-ranked route out of a pocket usually leaves through the nearest
    opening, which for a porter-domain pocket is the periplasmic cleft it came
    in by - that is a real route, not an artefact, but on its own it hides the
    one that matters. So keep the best of each exit direction rather than the
    top of the ranking.
    """
    mine = [r for r in rows if r[1] == ch]
    out = []
    for kind in ("up", "out", "down"):
        got = [r for r in mine if r[9].startswith(kind)]
        if got:
            out.append(min(got, key=lambda r: int(r[3])))
    return out


def build_scene(pid, s, rows, state):
    """One PDB + BILD + ChimeraX script showing all three states at once."""
    os.makedirs(OUT, exist_ok=True)
    chains = sorted(s.chains)
    pdb = os.path.join(OUT, f"{pid}_states.pdb")
    keep = []
    for a in s.protein_atoms:
        if not a.is_hydrogen:
            keep.append(a)
    with open(pdb, "w") as fh:
        fh.write(f"REMARK  {pid}: the trimer, chains "
                 + ", ".join(f"{c}={state.get(c, '?')}" for c in chains)
                 + "\n")
        fh.writelines(pdb_lines(keep))
        fh.writelines(pdb_lines([a for a in s.het_atoms
                                 if not a.is_hydrogen
                                 and a.resname.strip() not in ("HOH", "WAT")],
                                het=True))
        fh.write("END\n")

    bild = os.path.join(OUT, f"{pid}_states.bild")
    drawn = []
    with open(bild, "w") as fh:
        for ch in chains:
            st = state.get(ch, "")
            base = SCOL.get(st, "#7a8891")
            for i, r in enumerate(pick(rows, ch)):
                k = int(r[3])
                col = shade(base, 0.0 if i == 0 else 0.35 + 0.2 * i)
                fh.write(f".comment {ch} [{st}] cluster {k}: {r[9]}\n")
                fh.write(".color %.3f %.3f %.3f\n"
                         % tuple(int(col[j:j + 2], 16) / 255
                                 for j in (1, 3, 5)))
                for (x, y, z, rad) in cluster_spheres(
                        os.path.join(WORK_DIR, "caver_runs",
                                     f"{pid}_st{ch}"),
                        k):
                    fh.write(f".sphere {x:.3f} {y:.3f} {z:.3f} {rad:.2f}\n")
                drawn.append((ch, st, k, r[9], col))

    cxc = os.path.join(OUT, f"{pid}_states.cxc")
    with open(cxc, "w") as fh:
        fh.write(f"# ChimeraX: open this file, keeping it beside\n"
                 f"#   {os.path.basename(pdb)} and "
                 f"{os.path.basename(bild)}\n")
        fh.write(f"open {os.path.basename(pdb)}\n")
        fh.write("hide atoms\nshow cartoon\n")
        for ch in chains:
            fh.write(f"color #1/{ch} {SCOL.get(state.get(ch, ''), '#888888')}"
                     f"   # {state.get(ch, '?')}\n")
        fh.write("transparency #1 70 target c\n")
        fh.write(f"open {os.path.basename(bild)}\n")
        for (ch, st, k, call, col) in drawn:
            fh.write(f"# chain {ch} [{st}] cluster {k}: {call}\n")
        fh.write("select ligand\nstyle sel ball\ncolor sel byhetero\n")
        fh.write("~select\nset bgColor white\nlighting soft\n"
                 "graphics silhouettes true\nview\n")
    print(f"\n  scene: {os.path.relpath(pdb)}\n         "
          f"{os.path.relpath(bild)}\n         {os.path.relpath(cxc)}")
    for (ch, st, k, call, col) in drawn:
        print(f"      {ch} [{st:9}] cluster {k}: {call}")
    return drawn


def shade(hexcol, f):
    """Blend towards white by f, to separate several tunnels of one colour."""
    v = [int(hexcol[i:i + 2], 16) for i in (1, 3, 5)]
    v = [int(round(c + (255 - c) * f)) for c in v]
    return "#%02X%02X%02X" % tuple(v)


def pdb_lines(atoms, het=False):
    """PDB records. The altLoc column at 17 is easy to leave out and then
    every field after it is one place left, which no viewer will tell you -
    it just fails to recognise the residues and draws no cartoon."""
    tag = "HETATM" if het else "ATOM  "
    for i, a in enumerate(atoms, start=1):
        nm = a.name if len(a.name) >= 4 else f" {a.name:<3.3s}"
        yield (f"{tag}{i % 100000:5d} {nm:<4.4s} {a.resname:>3.3s} "
               f"{a.chain}{a.resseq:4d}    "
               f"{a.xyz[0]:8.3f}{a.xyz[1]:8.3f}{a.xyz[2]:8.3f}"
               f"  1.00{a.bfac:6.2f}          {a.element:>2.2s}\n")


def render_states(pid, s, drawn, state, axis_c, axis):
    """A side-on view with the pseudo-3-fold vertical, so up means up."""
    try:
        import pymol
        from pymol import cmd
    except ImportError:
        print("  (no PyMOL - skipping the rendered view)")
        return
    pdb = os.path.join(OUT, f"{pid}_states.pdb")
    pymol.finish_launching(["pymol", "-qc"])
    cmd.reinitialize()
    cmd.load(pdb, "trimer")
    cmd.hide("everything")
    cmd.show("cartoon", "polymer")
    for ch in sorted(s.chains):
        col = SCOL.get(state.get(ch, ""), "#888888")
        cmd.set_color(f"st{ch}", [int(col[j:j + 2], 16) / 255
                                  for j in (1, 3, 5)])
        cmd.color(f"st{ch}", f"polymer and chain {ch}")
    cmd.set("cartoon_transparency", 0.80, "polymer")

    for (ch, st, k, call, col) in drawn:
        sel = f"tun{ch}{k}"
        for (x, y, z, rad) in cluster_spheres(
                os.path.join(WORK_DIR, "caver_runs", f"{pid}_st{ch}"), k):
            cmd.pseudoatom(sel, pos=[float(x), float(y), float(z)],
                           vdw=float(rad))
        cmd.set_color(f"c{ch}{k}", [int(col[j:j + 2], 16) / 255
                                    for j in (1, 3, 5)])
        cmd.color(f"c{ch}{k}", sel)
        cmd.show("spheres", sel)
        # the route up the funnel is the one the mechanism is about; the
        # others are drawn faint so they read as context
        cmd.set("sphere_transparency", 0.0 if call.startswith("up") else 0.55,
                sel)
    lig = "hetatm and not resn HOH"
    if cmd.count_atoms(lig):
        cmd.show("sticks", lig)
        cmd.color("yellow", f"{lig} and elem C")
    cmd.hide("everything", "hydro")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("sphere_quality", 3)

    up = np.asarray(axis, float)
    ce = np.array([v for v in s.ca(sorted(s.chains)[1]).values()]).mean(0)
    back = ce - axis_c
    back = back - float(np.dot(back, up)) * up
    back /= float(np.linalg.norm(back))
    right = np.cross(up, back)
    # PyMOL's view matrix is world->camera in column-major order, so the
    # camera axes go in as columns; passing them as rows looks down the
    # 3-fold instead of across it
    cmd.set_view(list(np.array([right, up, back]).T.flatten())
                 + [0.0, 0.0, -260.0] + list(map(float, axis_c))
                 + [120.0, 400.0, -20.0])
    cmd.zoom("polymer", 3, complete=1)
    png = os.path.join(OUT, f"{pid}_states.png")
    cmd.png(png, width=2600, height=2200, dpi=300, ray=1)
    print(f"         {os.path.relpath(png)}")


def main():
    pid = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    only = sys.argv[2:] or None
    chan = load_channel()
    path = os.path.join(STRUCT_DIR, f"{pid}.pdb")
    if not os.path.exists(path):
        path = os.path.join(PDBDIR, f"{pid}.pdb")
    s = Structure(path)
    state = states_of(pid)
    axis_c, axis = T.trimer_axis(s)
    cf = T.ExactClearance([a for a in s.protein_atoms if not a.is_hydrogen])
    jar = rc.ensure_caver()
    if not jar:
        print("  CAVER not available")
        return

    print(f"=== CAVER by protomer state: {pid} ===")
    tell = ", ".join(f"{c}:{state.get(c, '?')}" for c in sorted(s.chains))
    print(f"    chains {tell}")
    rows = []
    for ch in sorted(s.chains):
        if only and ch not in only:
            continue
        seed, how = seed_of(s, ch)
        moved, clr, dist = free_point(cf, seed)
        if moved is None:
            print(f"  chain {ch} [{state.get(ch, '?')}]: nothing water-wide "
                  f"within reach of the {how} - skipped")
            continue
        if dist > 0.5:
            how += f", moved {dist:.1f} A to open space ({clr:.1f} A clear)"
        seed = moved
        cen = {c: coords([a for a in s.protein_atoms
                          if a.chain == c and a.name.strip() == "CA"]).mean(0)
               for c in s.chains}
        near = sorted(cen, key=lambda c: float(np.linalg.norm(
            cen[c] - cen[ch])))[:3]
        atoms = [a for a in s.protein_atoms
                 if not a.is_hydrogen and a.chain in set(near)]
        # its own tag: chain E already has a cached run from the tunnel
        # cross-check, seeded straight on the ligand centroid, and mixing
        # that in would leave one protomer of a three-state figure seeded by
        # a different rule from the other two
        tag = f"{pid}_st{ch}"
        cached = os.path.join(WORK_DIR, "caver_runs", tag, "out")
        if os.path.exists(os.path.join(cached, "analysis",
                                       "tunnel_profiles.csv")):
            res = reuse(cached)
            print(f"  chain {ch} [{state.get(ch, '?')}]: cached run, "
                  f"seeded on the {how}")
        else:
            print(f"  chain {ch} [{state.get(ch, '?')}]: running CAVER, "
                  f"seeded on the {how} ...")
            res = rc.run_one(jar, tag, atoms, seed, PROBE, SHELL_R, SHELL_D)
        if "error" in res:
            print(f"    CAVER failed - {res['error'][:120]}")
            continue

        d = os.path.join(res["outdir"], "data", "clusters_timeless")
        if not os.path.isdir(d):
            print("    no clusters written")
            continue
        hs = float(np.dot(seed - axis_c, axis))
        rs = float(np.linalg.norm((seed - axis_c)
                                  - hs * axis))
        for k in range(1, N_CLUSTERS + 1):
            pts = []
            for f in sorted(os.listdir(d)):
                if not f.startswith(f"tun_cl_{k:03d}_"):
                    continue
                for ln in open(os.path.join(d, f)):
                    if ln.startswith(("ATOM", "HETATM")):
                        pts.append([float(ln[30:38]), float(ln[38:46]),
                                    float(ln[46:54]),
                                    float(ln[54:].split()[0])])
            if not pts:
                continue
            P = np.array(pts)[:, :3]
            rad = np.array(pts)[:, 3]
            arc = float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum())
            h = float(np.dot(P[-1] - axis_c, axis))
            r = float(np.linalg.norm((P[-1] - axis_c) - h * axis))
            call = exit_call(h - hs, r - rs)
            rows.append([pid, ch, state.get(ch, ""), k, len(P), fmt(arc),
                         fmt(rad.min()), fmt(h - hs), fmt(r - rs), call, how])
            print(f"    cluster {k}: {arc:5.1f} A long, bottleneck "
                  f"{rad.min():.2f} A, rises {h - hs:+6.1f} A, "
                  f"{r - rs:+6.1f} A outward  ->  {call}")

    write_csv(os.path.join(TABLES, "caver_states.csv"),
              ["pdb", "chain", "state", "cluster", "spheres", "length_A",
               "bottleneck_A", "rise_along_axis_A", "outward_step_A",
               "exit", "seed"], rows)
    print("\nwrote results/tables/caver_states.csv")
    if rows and len({r[1] for r in rows}) > 1:
        drawn = build_scene(pid, s, rows, state)
        render_states(pid, s, drawn, state, axis_c, axis)


if __name__ == "__main__":
    main()
