#!/usr/bin/env python3
"""Build a viewable scene of the CAVER tunnels: PDB + ChimeraX + PyMOL files.

CAVER writes each tunnel cluster as a string of spheres (data/clusters_timeless
/tun_cl_XXX_1.pdb), with the sphere radius in the last numeric column and the
coordinates in the frame of the structure CAVER was given. This packs that back
together into one file a viewer can open: the protein CAVER actually saw, the
ligands of that protomer, and the top clusters as sphere pseudo-atoms carrying
their radius in the B-factor column, one chain per cluster.

Rendering the spheres at their stored radius is what makes the tunnel look like
a tunnel; both command files do that, so the shape on screen is the measured
one rather than a decoration.

Usage:  python3 scripts/caver_scene.py [run_tag ...]   (default: our two models)
Writes results/structures/<tag>_tunnels.pdb, .cxc and .pml
"""
from __future__ import annotations

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import STRUCT_DIR, TABLES, WORK_DIR
from published_pockets import PDBDIR

N_CLUSTERS = 3
OUT = os.path.join(os.path.dirname(TABLES), "structures")
# poster colours, cluster 1 first
COL = [("#CA0FC1", "magenta-ish"), ("#0A9DA0", "teal"), ("#C68B3C", "amber"),
       ("#0F9C1B", "green"), ("#2F54FF", "blue")]
DEFAULT = ["Amp_MexB_20260826_E", "MexB_DDM_3_20260730_E"]


def cluster_spheres(run, k):
    """[(x, y, z, radius)] of cluster k, or [] if CAVER did not write it."""
    d = os.path.join(run, "out", "data", "clusters_timeless")
    hit = [f for f in sorted(os.listdir(d))
           if re.match(rf"tun_cl_0*{k}_\d+\.pdb$", f)] if os.path.isdir(d) \
        else []
    pts = []
    for f in hit:
        for ln in open(os.path.join(d, f)):
            if not ln.startswith(("ATOM", "HETATM")):
                continue
            try:
                x, y, z = (float(ln[30:38]), float(ln[38:46]),
                           float(ln[46:54]))
            except ValueError:
                continue
            tail = ln[54:].split()
            r = float(tail[0]) if tail else 1.0
            pts.append((x, y, z, r))
    return pts


def source_pdb(tag):
    """The structure the run came from, for its ligands."""
    pid = re.sub(r"_(protein|withlig)?_?[A-Z]$", "", tag)
    for base in (STRUCT_DIR, PDBDIR):
        p = os.path.join(base, f"{pid}.pdb")
        if os.path.exists(p):
            return p, pid
    return None, pid


def build(tag):
    run = os.path.join(WORK_DIR, "caver_runs", tag)
    prot = os.path.join(run, "pdbs", "1.pdb")
    if not os.path.exists(prot):
        print(f"  {tag}: no CAVER input pdb - skipped")
        return
    lines = [ln for ln in open(prot) if ln.startswith("ATOM")]
    chains = {ln[21] for ln in lines}

    src, pid = source_pdb(tag)
    hets = []
    if src:
        for ln in open(src):
            if (ln.startswith("HETATM") and ln[21] in chains
                    and ln[17:20].strip() not in ("HOH", "WAT")):
                hets.append(ln)

    got = []
    body = []
    serial = len(lines) + len(hets)
    for k in range(1, N_CLUSTERS + 1):
        pts = cluster_spheres(run, k)
        if not pts:
            continue
        got.append((k, len(pts), max(p[3] for p in pts),
                    min(p[3] for p in pts)))
        ch = "TUVWX"[len(got) - 1]
        for i, (x, y, z, r) in enumerate(pts, start=1):
            serial += 1
            body.append(f"HETATM{serial % 100000:5d}  O   TUN {ch}"
                        f"{k:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{r:6.2f}"
                        f"           O\n")
    if not got:
        print(f"  {tag}: no tunnel clusters found - skipped")
        return

    os.makedirs(OUT, exist_ok=True)
    pdb = os.path.join(OUT, f"{tag}_tunnels.pdb")
    with open(pdb, "w") as fh:
        fh.write(f"REMARK  {pid}: CAVER 3.0.3 tunnel clusters "
                 f"{', '.join(str(g[0]) for g in got)}\n")
        fh.write("REMARK  TUN pseudo-atoms carry the local tunnel radius in "
                 "the B-factor column;\n")
        fh.write("REMARK  render them as spheres sized by B-factor to see "
                 "the measured shape.\n")
        fh.writelines(lines)
        fh.writelines(hets)
        fh.writelines(body)
        fh.write("END\n")

    tunchains = "".join("TUVWX"[i] for i in range(len(got)))

    # BILD: explicit spheres at explicit radii. ChimeraX draws these exactly
    # as written, with no dependence on how a viewer maps a B-factor column
    # onto an atom radius - which is the part that silently goes wrong.
    bild = os.path.join(OUT, f"{tag}_tunnels.bild")
    with open(bild, "w") as fh:
        for i, (k, n_sp, rmax, rmin) in enumerate(got):
            c = COL[i][0]
            fh.write(f".comment cluster {k}: {n_sp} spheres, "
                     f"r {rmin:.2f}-{rmax:.2f} A\n")
            fh.write(".color %.3f %.3f %.3f\n"
                     % tuple(int(c[j:j + 2], 16) / 255 for j in (1, 3, 5)))
            for (x, y, z, r) in cluster_spheres(run, k):
                fh.write(f".sphere {x:.3f} {y:.3f} {z:.3f} {r:.2f}\n")

    cxc = os.path.join(OUT, f"{tag}_tunnels.cxc")
    with open(cxc, "w") as fh:
        fh.write(f"# ChimeraX: open this file (File > Open, or\n"
                 f"#   'open {os.path.basename(cxc)}' in the command line).\n"
                 f"# Keep it beside {os.path.basename(pdb)} and "
                 f"{os.path.basename(bild)}.\n")
        fh.write(f"open {os.path.basename(pdb)}\n")
        fh.write("hide atoms\nshow cartoon\n")
        fh.write("color #1 #D6DEE2\ncolor #1/E #6699CC\n")
        fh.write("transparency #1 60 target c\n")
        fh.write("transparency #1/E 35 target c\n")
        fh.write(f"open {os.path.basename(bild)}\n")
        for i, (k, n_sp, rmax, rmin) in enumerate(got):
            fh.write(f"# model #2 sub-model: cluster {k}, {COL[i][1]}, "
                     f"r {rmin:.2f}-{rmax:.2f} A\n")
        fh.write("select ligand\nstyle sel ball\ncolor sel byhetero\n")
        fh.write("~select\nset bgColor white\nlighting soft\ngraphics "
                 "silhouettes true\nview #2\n")
        fh.write("# the TUN atoms in the PDB carry the same radii in their "
                 "B-factor column,\n# if you would rather drive them "
                 "yourself: show #1/" + tunchains + " atoms\n")

    pml = os.path.join(OUT, f"{tag}_tunnels.pml")
    with open(pml, "w") as fh:
        fh.write(f"load {os.path.basename(pdb)}, scene\n")
        fh.write("hide everything\nshow cartoon, polymer\n")
        fh.write("color grey80, polymer\ncolor skyblue, polymer and chain E\n")
        fh.write("select tun, resn TUN\n")
        fh.write("alter tun, vdw=b\nrebuild\nshow spheres, tun\n")
        for i, (k, n, rmax, rmin) in enumerate(got):
            fh.write(f"set_color cav{k}, [{int(COL[i][0][1:3], 16)/255:.3f},"
                     f"{int(COL[i][0][3:5], 16)/255:.3f},"
                     f"{int(COL[i][0][5:7], 16)/255:.3f}]\n")
            fh.write(f"color cav{k}, resn TUN and chain {'TUVWX'[i]}\n")
        fh.write("show sticks, not polymer and not resn TUN\n")
        fh.write("util.cbay('not polymer and not resn TUN')\n")
        fh.write("bg_color white\nset ray_opaque_background, 1\norient\n")

    print(f"  {tag}: {len(lines)} protein atoms, {len(hets)} ligand atoms, "
          f"{len(got)} clusters")
    for i, (k, n, rmax, rmin) in enumerate(got):
        print(f"      cluster {k} -> chain {'TUVWX'[i]}, {n} spheres, "
              f"radius {rmin:.2f}-{rmax:.2f} A, {COL[i][1]}")
    print(f"      {os.path.relpath(pdb)}\n      {os.path.relpath(bild)}\n"
          f"      {os.path.relpath(cxc)}\n      {os.path.relpath(pml)}")


def render(tag):
    """A poster-ready view of the same scene, if PyMOL is installed."""
    try:
        import pymol
        from pymol import cmd
    except ImportError:
        print("  (no PyMOL - skipping the rendered view)")
        return
    pdb = os.path.join(OUT, f"{tag}_tunnels.pdb")
    if not os.path.exists(pdb):
        return
    pymol.finish_launching(["pymol", "-qc"])
    cmd.reinitialize()
    cmd.load(pdb, "scene")
    cmd.hide("everything")
    # only the protomer the tunnels belong to, and see-through, or the
    # tunnels are buried behind its own porter domain
    cmd.show("cartoon", "polymer and chain E")
    cmd.color("skyblue", "polymer and chain E")
    cmd.set("cartoon_transparency", 0.70, "polymer and chain E")
    for i, ch in enumerate("TUVWX"):
        sel = f"resn TUN and chain {ch}"
        if cmd.count_atoms(sel) == 0:
            continue
        cmd.set_color(f"cav{i}", [int(COL[i][0][j:j + 2], 16) / 255
                                  for j in (1, 3, 5)])
        cmd.alter(sel, "vdw=b")
        cmd.color(f"cav{i}", sel)
        cmd.show("spheres", sel)
    cmd.rebuild()
    lig = "not polymer and not resn TUN"
    if cmd.count_atoms(lig):
        cmd.show("sticks", lig)
        cmd.color("yellow", f"{lig} and elem C")
        cmd.set("stick_radius", 0.28)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("specular", 0.2)
    cmd.set("sphere_quality", 3)
    cmd.hide("everything", "hydro")
    cmd.orient("resn TUN")
    cmd.zoom("resn TUN or (%s)" % lig, 12, complete=1)
    cmd.clip("slab", 400)
    png = os.path.join(OUT, f"{tag}_tunnels.png")
    cmd.png(png, width=2400, height=2000, dpi=300, ray=1)
    print(f"      {os.path.relpath(png)}")


def main():
    tags = sys.argv[1:] or DEFAULT
    print("=== CAVER tunnel scenes ===")
    for t in tags:
        build(t)
        render(t)


if __name__ == "__main__":
    main()
