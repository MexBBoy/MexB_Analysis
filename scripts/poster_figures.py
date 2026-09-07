#!/usr/bin/env python3
"""Poster-scale panels that extend the Combio poster's story.

Sized and typed for A0 at ~1 m viewing distance, in the poster's own palette
(dark teal #104862 headers, mint ground, protomers cyan/magenta/green for
access/binding/extrusion) so they sit beside the existing panels rather than
looking imported.

Each panel adds something the poster does not already show:
  P1  the pocket-polarity renderings, made quantitative
  P2  the detergent-occupancy problem the introduction raises, measured
  P3  the exit route, and where passage is actually hard - a dimension the
      poster does not cover at all
"""
from __future__ import annotations

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import CXDIR, FIGURES, TABLES

# ---- poster palette, sampled from Combio_Poster_20260828.pdf
TEAL = "#104862"          # header bars and titles
MINT = "#C8FEED"          # poster ground
INK = "#16202a"
INK2 = "#4a5a66"
GRID = "#dde5e8"
# protomer identity, darkened just enough to clear 3:1 on white while
# staying recognisably the poster's cyan / magenta / green
ACCESS, BINDING, EXTRUSION = "#0A9DA0", "#CA0FC1", "#0F9C1B"
# the poster's own hydrophobic/hydrophilic ramp
APOLAR, POLAR = "#C68B3C", "#0E9AA0"
WARN = "#B26A00"

OUT = os.path.join(FIGURES, "poster")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 17,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.edgecolor": "#9fb0b8", "axes.linewidth": 1.6,
    "axes.labelcolor": INK, "axes.labelsize": 18,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 16, "ytick.labelsize": 16,
    "xtick.major.width": 1.4, "ytick.major.width": 1.4,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 1.2,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 16,
    "lines.linewidth": 3.4,
})

# each substrate in the colour the poster gives it, sampled from the
# conserved-residues legend and table header of Combio_Poster_20260828.pdf,
# then darkened at constant hue to clear 4.5:1 on white. 3W9I is not on the
# poster; it takes the pink the legend uses for DDM generally.
LIGCOL = {
    "Amp_MexB_20260826": "#078A08",    # poster #1EFF21 ampicillin green
    "MexB_DDM_3_20260730": "#CF13CF",  # poster #FF29FF DDM #1 magenta
    "2V50": "#986598",                 # poster #FFB0FF DDM #2 pink
    "3W9I": "#C24A8B",                 # poster #FF6DBC legend DDM pink
    "6IIA": "#2F54FF",                 # poster #3E61FF LMNG blue
    "21FO": "#A1685E",                 # poster #FFAB9C CYMAL-7 salmon
    "3W9J": "#767676",                 # poster #A3A3A3 EPI grey
    "21FP": "#000000",                 # poster black chloramphenicol
}


def tint(hexcol, f):
    """Blend a colour towards white (f>0) or towards black (f<0)."""
    r, g, b = (int(hexcol[i:i + 2], 16) for i in (1, 3, 5))
    if f >= 0:
        return "#%02X%02X%02X" % tuple(round(c + (255 - c) * f)
                                       for c in (r, g, b))
    return "#%02X%02X%02X" % tuple(round(c * (1 + f)) for c in (r, g, b))


STATE_COLOR = {"Access": ACCESS, "Binding": BINDING, "Extrusion": EXTRUSION}
SHORT = {"Amp_MexB_20260826": "Ampicillin", "MexB_DDM_3_20260730": "DDM"}


def R(name):
    p = os.path.join(TABLES, name)
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh))


def title(fig, text, sub=None, y=0.995):
    fig.text(0.005, y, text, ha="left", va="top", fontsize=26,
             fontweight="bold", color=TEAL)
    if sub:
        fig.text(0.005, y - 0.062, sub, ha="left", va="top", fontsize=16,
                 color=INK2)


def callout(fig, x, y, big, small, color=TEAL, size=44):
    fig.text(x, y, big, ha="left", va="top", fontsize=size,
             fontweight="bold", color=color)
    fig.text(x, y - 0.075, small, ha="left", va="top", fontsize=15,
             color=INK2, linespacing=1.35)



def statcell(fig, rect, big, caption, color, big_size=58):
    """Big number over its caption, inside a reserved rectangle."""
    ax = fig.add_axes(rect)
    ax.axis("off")
    ax.text(0, 1.0, big, ha="left", va="top", fontsize=big_size,
            fontweight="bold", color=color, transform=ax.transAxes)
    ax.text(0, 0.30, caption, ha="left", va="top", fontsize=15.5,
            color=INK2, transform=ax.transAxes, linespacing=1.5)
    return ax


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=300,
                    bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)
    print(f"  wrote results/figures/poster/{name}.png / .svg")


def tunnel_profile():
    """(depth, radius) along the reference entry channel, in Angstrom.

    Straight from the tunnel trace tunnels.py writes: the B-factor column is
    the local radius and depth is arc length back from the periplasmic mouth,
    the same coordinate the ligands are placed on.
    """
    p = os.path.join(CXDIR,
                     "Amp_MexB_20260826_protein_E_ZZ72000_t1_tunnel.pdb")
    if not os.path.exists(p):
        return None
    P, B = [], []
    for ln in open(p):
        if ln.startswith(("ATOM", "HETATM")):
            P.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            B.append(float(ln[60:66]))
    if len(P) < 10:
        return None
    P = np.asarray(P, float)
    arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0),
                                                          axis=1))])
    depth = arc[-1] - arc
    o = np.argsort(depth)
    return depth[o], np.asarray(B, float)[o]


def state_of(struct, chain, states):
    for r in states:
        if r["structure"] == struct and r["chain"] == chain:
            return r["state_call"]
    return ""


# ------------------------------------------------------------------- P1
def panel_pockets():
    comp = R("pocket_composition.csv")
    states = R("states.csv")
    if not comp:
        return
    fig = plt.figure(figsize=(10.4, 6.3))
    gs = fig.add_gridspec(1, 2, left=0.095, right=0.985, top=0.575,
                          bottom=0.26, wspace=0.34)
    specs = [("pct_apolar", "apolar side-chain\natoms (%)", "{:.0f}%"),
             ("mean_KD", "mean Kyte–Doolittle\nhydropathy", "{:+.2f}")]
    for k, (key, ylab, fs) in enumerate(specs):
        ax = fig.add_subplot(gs[0, k])
        for i, pocket in enumerate(("PBP", "DBP")):
            vals, cols = [], []
            for r in comp:
                if r["pocket"] != pocket:
                    continue
                vals.append(float(r[key]))
                cols.append(STATE_COLOR.get(
                    state_of(r["structure"], r["chain"], states), INK2))
            x = np.full(len(vals), i) + np.linspace(-.16, .16, len(vals))
            ax.scatter(x, vals, s=200, c=cols, zorder=3,
                       edgecolor="white", linewidth=2.2)
            m = float(np.mean(vals))
            ax.plot([i - .30, i + .30], [m, m], color=INK, linewidth=4,
                    solid_capstyle="round", zorder=4)
            ax.annotate(fs.format(m), (i, m), textcoords="offset points",
                        xytext=(0, 16), ha="center", fontsize=24,
                        fontweight="bold",
                        color=APOLAR if pocket == "DBP" else POLAR)
        if key == "mean_KD":
            ax.axhline(0, color="#9fb0b8", linewidth=1.6,
                       linestyle=(0, (5, 4)))
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Proximal\n(PBP)", "Distal\n(DBP)"],
                           fontsize=18)
        ax.set_xlim(-.55, 1.55)
        ax.set_ylabel(ylab, labelpad=10)
        ax.margins(y=.34)
        ax.grid(axis="x", visible=False)
        ax.set_axisbelow(True)
    title(fig, "The distal pocket is the hydrophobic one — measured",
          "Atom composition of every pocket-lining side chain, scored "
          "independently in all six protomers of both maps.")
    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=13,
                          markerfacecolor=c, markeredgecolor="white",
                          markeredgewidth=2, label=s)
               for s, c in STATE_COLOR.items()]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(0.55, 0.715), fontsize=16,
               handletextpad=0.35, columnspacing=1.6)
    fig.text(0.095, 0.075,
             "Distal: 8 aromatic residues.   Proximal: none.",
             fontsize=16, color=INK)
    fig.text(0.095, 0.028,
             "The ordering is identical in every protomer of both "
             "reconstructions.", fontsize=15, color=INK2)
    save(fig, "P1_pocket_chemistry")


# ------------------------------------------------------------------- P2
def panel_occlusion():
    vol = [r for r in R("pocket_volumes.csv")
           if r["volume_ligands_stripped_A3"]]
    states = R("states.csv")
    if not vol:
        return
    fig = plt.figure(figsize=(10.2, 7.0))
    title(fig, "Detergent fills the pocket it obscures",
          "The problem raised in the introduction, measured: free volume "
          "within 16 Å of the substrate site.")
    statcell(fig, [0.095, 0.635, 0.42, 0.175], "98%",
             "of the binding-protomer site is taken by\n"
             "the three DDM molecules    2018 → 36 Å³", BINDING)
    # Ampicillin displaces a fixed 387 A^3 at every sphere radius tested,
    # but its PERCENTAGE runs 27%->12% as the sphere grows, because the
    # denominator grows and the ligand does not. Quote the volume.
    statcell(fig, [0.575, 0.635, 0.42, 0.175], "387 Å³",
             "displaced by ampicillin in the\n"
             "equivalent protomer", TEAL, big_size=52)

    ax = fig.add_axes([0.095, 0.165, 0.88, 0.365])
    labs, free, occ, cols = [], [], [], []
    for r in vol:
        labs.append(f"{SHORT[r['structure']][:3]} {r['chain']}")
        free.append(float(r["volume_ligands_stripped_A3"]))
        occ.append(float(r["volume_with_ligands_A3"]))
        cols.append(STATE_COLOR.get(
            state_of(r["structure"], r["chain"], states), INK2))
    x = np.arange(len(vol)); w = 0.38
    ax.bar(x - w / 2 - .012, free, w, color="#cbd6db",
           edgecolor="white", linewidth=1.6, label="ligand removed")
    ax.bar(x + w / 2 + .012, occ, w, color=cols,
           edgecolor="white", linewidth=1.6, label="ligand in place")
    for i, (f, o) in enumerate(zip(free, occ)):
        pct = 100 * (f - o) / f if f else 0
        if pct > 1:
            ax.annotate(f"−{f - o:.0f} Å³", (i, max(f, o)),
                        textcoords="offset points", xytext=(0, 13),
                        ha="center", fontsize=20, fontweight="bold",
                        color=WARN)
    ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=16)
    ax.set_ylabel("free volume (Å³)", fontsize=17)
    ax.set_ylim(0, max(free) * 1.22)
    ax.grid(axis="x", visible=False); ax.set_axisbelow(True)
    # direct labels on the first pair instead of a legend box: nothing to
    # collide with, and it reads without a key at poster distance
    # key goes in the reserved band between the stat row and the chart,
    # where neither the bars nor the % callouts can reach it
    fig.legend(*ax.get_legend_handles_labels(), loc="upper center", ncol=2,
               bbox_to_anchor=(0.53, 0.588), fontsize=15, handlelength=1.5,
               handleheight=1.1, columnspacing=2.0)
    fig.text(0.095, 0.075,
             "The DDM figure is a percentage because three molecules "
             "saturate the pocket: 97.8–98.8% however the site is measured. "
             "Ampicillin\ndisplaces a constant 387 Å³, so it is quoted as a "
             "volume. All three DDM are present in both half maps; LMT2003 "
             "is the weakest.",
             fontsize=14.5, color=INK2, va="top", linespacing=1.5)
    save(fig, "P2_detergent_occlusion")


# ------------------------------------------------------------------- P3
def read_trace(path):
    d, b = [], []
    for line in open(path):
        if line.startswith("HETATM"):
            d.append((float(line[30:38]), float(line[38:46]),
                      float(line[46:54])))
            b.append(float(line[60:66]))
    P = np.array(d)
    s = np.concatenate([[0], np.cumsum(
        np.linalg.norm(np.diff(P, axis=0), axis=1))])
    return s, np.array(b)


def panel_exit_route():
    cd = [r for r in R("caverdock_profile.csv") if r["bound"] == "lb"]
    if not cd:
        return
    arc = np.array([float(r["position_along_tunnel_A"]) for r in cd])
    E = np.array([float(r["energy_kcal_mol"]) for r in cd])
    rad = np.array([float(r["tunnel_radius_A"]) for r in cd])
    i, j = int(np.argmax(E)), int(np.argmin(rad))

    fig = plt.figure(figsize=(10.6, 7.0))
    title(fig, "Where the exit route is hard is not where it is narrow",
          "Ampicillin pulled along the chain E route out of the binding "
          "pocket to bulk solvent.")

    # left column: two stacked plots. right column: the number and its text.
    ax = fig.add_axes([0.095, 0.455, 0.50, 0.245])
    ax2 = fig.add_axes([0.095, 0.185, 0.50, 0.245], sharex=ax)

    ax.plot(arc, E, color=BINDING, solid_capstyle="round")
    ax.scatter([arc[i]], [E[i]], s=210, color=BINDING, zorder=5,
               edgecolor="white", linewidth=2.4)
    ax.set_ylabel("binding energy\n(kcal/mol)", fontsize=16)
    ax.margins(y=.34); ax.tick_params(labelbottom=False)
    ax.annotate("energy barrier", (arc[i], E[i]),
                textcoords="offset points", xytext=(12, 2), fontsize=16,
                fontweight="bold", color=BINDING)

    ax2.plot(arc, rad, color=TEAL, solid_capstyle="round")
    ax2.scatter([arc[j]], [rad[j]], s=210, color=TEAL, zorder=5,
                edgecolor="white", linewidth=2.4)
    ax2.set_ylabel("tunnel radius\n(Å)", fontsize=16)
    ax2.set_xlabel("distance along the exit route (Å)", fontsize=16)
    ax2.margins(y=.36)
    ax2.annotate("narrowest point", (arc[j], rad[j]),
                 textcoords="offset points", xytext=(-14, 10), ha="right",
                 fontsize=16, fontweight="bold", color=TEAL)
    for a_ in (ax, ax2):
        a_.axvline(arc[i], color=BINDING, linewidth=2, linestyle=(0, (4, 4)),
                   alpha=.6)
        a_.axvline(arc[j], color=TEAL, linewidth=2, linestyle=(0, (4, 4)),
                   alpha=.6)
        a_.grid(axis="x", visible=False); a_.set_axisbelow(True)

    statcell(fig, [0.645, 0.52, 0.34, 0.19], f"{arc[j] - arc[i]:.0f} Å",
             "between the energy barrier and the\ngeometric constriction — "
             "radius alone\nnames the wrong rate-limiting residues",
             WARN, big_size=54)
    tx = fig.add_axes([0.645, 0.175, 0.34, 0.245]); tx.axis("off")
    tx.text(0, 1.0,
            f"bottleneck   {rad[j]:.2f} Å\n"
            f"barrier         +{E[i] - E[0]:.1f} kcal/mol\n\n"
            "Exit is through the PC1/PC2 periplasmic\ncleft. The bottleneck "
            "is confirmed\nindependently by CAVER 3.0.3 — the\ntwo agree "
            "to 0.01 Å.",
            ha="left", va="top", fontsize=15.5, color=INK2,
            transform=tx.transAxes, linespacing=1.55)

    fig.text(0.095, 0.075,
             "Docking is a lower bound at low exhaustiveness: the separation "
             "is robust, the barrier value is provisional.\nThe ampicillin "
             "pose and the constriction residues are themselves weakly "
             "resolved in the density.",
             fontsize=14.5, color=WARN, va="top", linespacing=1.5)
    save(fig, "P3_exit_route")



# ------------------------------------------------------------------- P4
def panel_ligand_size():
    rows = R("published_pockets.csv")
    if not rows:
        return
    NAME = {"21FP": "chloramphenicol", "Amp_MexB_20260826": "ampicillin",
            "2V50": "DDM (2V50)", "3W9I": "DDM (3W9I)", "21FO": "CYMAL-7",
            "3W9J": "EPI", "6IIA": "LMNG",
            "MexB_DDM_3_20260730": "DDM \u00d73", "6T7S": "apo"}
    OURS = {"Amp_MexB_20260826", "MexB_DDM_3_20260730"}

    def num(r, k):
        try:
            return float(r[k])
        except (KeyError, TypeError, ValueError):
            return float("nan")

    bound = [r for r in rows if int(r["ligand_heavy_atoms"]) > 0
             and np.isfinite(num(r, "depth_from_entrance_A"))]
    if not bound:
        return
    D = np.array([num(r, "depth_from_entrance_A") for r in bound])
    V = np.array([num(r, "volume_r16_A3") for r in bound])
    L = np.array([int(r["ligand_heavy_atoms"]) for r in bound], float)
    rp = float(np.corrcoef(D, V)[0, 1])

    fig = plt.figure(figsize=(10.6, 6.6))
    title(fig, "The pocket does not enlarge, wherever the ligand sits",
          "Every published substrate- or detergent-bound MexB structure, "
          "measured in one common frame.")
    ax = fig.add_axes([0.095, 0.275, 0.60, 0.40])

    for r in bound:
        x = num(r, "depth_from_entrance_A"); y = num(r, "volume_r16_A3")
        mine = r["pdb"] in OURS
        col = LIGCOL.get(r["pdb"], TEAL)
        # marker area tracks ligand size, so both variables stay visible
        ax.scatter([x], [y], s=90 + 3.2 * int(r["ligand_heavy_atoms"]),
                   color=col, zorder=4, edgecolor="white", linewidth=2.4,
                   marker="D" if mine else "o")
        OFF = {"21FP": (17, -6, "left"), "Amp_MexB_20260826": (14, -12, "left"),
               "2V50": (17, -6, "left"), "3W9I": (17, -6, "left"),
               "21FO": (0, -30, "center"), "3W9J": (0, -30, "center"),
               "6IIA": (-16, -6, "right"),
               "MexB_DDM_3_20260730": (0, 24, "center")}
        dx, dy, ha = OFF.get(r["pdb"], (0, -28, "center"))
        ax.annotate(NAME.get(r["pdb"], r["pdb"]), (x, y),
                    textcoords="offset points", xytext=(dx, dy),
                    ha=ha, fontsize=13, color=col)

    xs = np.linspace(D.min() - 4, D.max() + 4, 10)
    ax.plot(xs, np.polyval(np.polyfit(D, V, 1), xs), color=TEAL,
            linewidth=2, linestyle=(0, (5, 4)), alpha=.55, zorder=1)
    ax.set_xlabel("depth into the porter domain from the periplasmic "
                  "entrance (\u00c5)", labelpad=26)
    ax.set_ylabel("ligand-free pocket volume (\u00c5\u00b3)")
    ax.set_xlim(25, 80); ax.margins(y=.30)
    ax.set_xticks([30, 40, 50, 60, 70])
    ax.grid(axis="x", visible=False); ax.set_axisbelow(True)
    ax.annotate("\u2190 towards the entrance", (0.0, 0.0),
                xycoords="axes fraction", xytext=(2, -42),
                textcoords="offset points", ha="left", fontsize=13,
                color=INK2, fontstyle="italic")
    ax.annotate("deeper into the porter domain \u2192", (1.0, 0.0),
                xycoords="axes fraction", xytext=(-2, -42),
                textcoords="offset points", ha="right", fontsize=13,
                color=INK2, fontstyle="italic")

    tx = fig.add_axes([0.735, 0.245, 0.25, 0.44]); tx.axis("off")
    tx.text(0, 1.0, f"r = {rp:+.2f}", ha="left", va="top", fontsize=46,
            fontweight="bold", color=TEAL, transform=tx.transAxes)
    tx.text(0, 0.79,
            "between pocket volume and\nhow deep the ligand sits,\n"
            "across eight structures.\n\n"
            "Ligands sit 33 to 63 \u00c5 in,\nover a 5.3\u00d7 range of\n"
            "size, with no systematic\nchange in pocket volume\n"
            "(r = +0.05 against ligand\nsize).",
            ha="left", va="top", fontsize=14, color=INK2,
            transform=tx.transAxes, linespacing=1.55)

    fig.text(0.095, 0.085,
             "Depth is arc length back from the periplasmic mouth along the "
             "widest ligand-free entry channel of the reference protomer, "
             "to the ligand centroid. The proximal and\ndistal pockets both "
             "lie at the far end of that channel, so depth separates pocket "
             "from entry cleft but not one pocket from the other.\nMarker "
             "area tracks ligand size. Each "
             "protomer superposed on 39 pocket-lining C\u03b1 of one "
             "reference, so the measuring sphere sits\nidentically in every "
             "structure. Engineered MexB chimeras excluded.",
             fontsize=14, color=INK2, va="top", linespacing=1.5)
    save(fig, "P4_ligand_size_vs_pocket")



# ------------------------------------------------------------------- P5
def panel_protomer_states():
    """Pocket volume in all 39 protomers, grouped by conformational state."""
    rows = R("protomer_pockets.csv")
    if not rows:
        return
    OURS = ("Amp_MexB_20260826", "MexB_DDM_3_20260730")
    STATES = ("Access", "Binding", "Extrusion")

    def num(r, k):
        try:
            return float(r[k])
        except (KeyError, TypeError, ValueError):
            return float("nan")

    rows = [r for r in rows if r["state_call"] in STATES
            and np.isfinite(num(r, "free_volume_r16_A3"))]
    if not rows:
        return

    fig = plt.figure(figsize=(10.6, 6.8))
    nprot = len(rows)
    nstruct = len({r["pdb"] for r in rows})
    title(fig, "Every protomer, not just the ligand-bound one",
          f"Free volume at the substrate site in all {nprot} protomers of "
          f"{nstruct} structures, measured in one common frame.")
    ax = fig.add_axes([0.095, 0.245, 0.58, 0.44])

    rng = np.random.default_rng(0)
    stats = {}
    for i, st in enumerate(STATES):
        grp = [r for r in rows if r["state_call"] == st]
        v = np.array([num(r, "free_volume_r16_A3") for r in grp])
        stats[st] = v
        pub = [(j, r) for j, r in enumerate(grp) if r["pdb"] not in OURS]
        mine = [(j, r) for j, r in enumerate(grp) if r["pdb"] in OURS]
        jit = rng.uniform(-.17, .17, len(grp))
        col = STATE_COLOR[st]
        for j, r in pub:
            ax.scatter([i + jit[j]], [v[j]], s=110, color=col, alpha=.55,
                       zorder=3, edgecolor="white", linewidth=1.6)
        # fixed offsets: the six labels sit close to each other and to the
        # mean bars, so automatic placement collides
        OFF = {"MexB_DDM_3_20260730|D": (16, -17), "Amp_MexB_20260826|D": (16, -18),
               "MexB_DDM_3_20260730|E": (10, 15), "Amp_MexB_20260826|E": (16, -14),
               "MexB_DDM_3_20260730|F": (10, 15), "Amp_MexB_20260826|F": (14, -18)}
        for j, r in mine:
            ax.scatter([i + jit[j]], [v[j]], s=230, color=col, zorder=5,
                       marker="D", edgecolor=INK, linewidth=1.8)
            dx, dy = OFF.get(f"{r['pdb']}|{r['chain']}", (16, -5))
            ax.annotate(SHORT.get(r["pdb"], r["pdb"]) + " " + r["chain"],
                        (i + jit[j], v[j]), textcoords="offset points",
                        xytext=(dx, dy), ha="left", fontsize=12.5, color=INK,
                        zorder=7)
        m = float(v.mean())
        ax.plot([i - .33, i + .33], [m, m], color=INK, linewidth=4,
                solid_capstyle="round", zorder=6)
        ax.annotate(f"{m:.0f}", (i - .33, m), textcoords="offset points",
                    xytext=(-6, -6), ha="right", fontsize=20,
                    fontweight="bold", color=col)

    ax.set_xticks(range(len(STATES)))
    ax.set_xticklabels([f"{st}\n(n = {len(stats[st])})" for st in STATES],
                       fontsize=18)
    ax.set_xlim(-.62, len(STATES) - .38)
    ax.set_ylabel("free volume at the substrate site (\u00c5\u00b3)",
                  labelpad=10)
    ax.margins(y=.16)
    ax.grid(axis="x", visible=False); ax.set_axisbelow(True)

    acc, bind, ext = (stats[s] for s in STATES)
    # largest departure of one of ours from the published spread for its state
    worst = 0.0
    for r in rows:
        if r["pdb"] not in OURS:
            continue
        same = [num(x, "free_volume_r16_A3") for x in rows
                if x["state_call"] == r["state_call"] and x["pdb"] not in OURS]
        if len(same) > 2:
            z = abs(num(r, "free_volume_r16_A3") - np.mean(same)) / \
                np.std(same, ddof=1)
            worst = max(worst, float(z))
    tx = fig.add_axes([0.715, 0.245, 0.27, 0.44]); tx.axis("off")
    tx.text(0, 1.0, f"{bind.mean() / acc.mean():.1f}\u00d7",
            ha="left", va="top", fontsize=46, fontweight="bold",
            color=BINDING, transform=tx.transAxes)
    tx.text(0, 0.80,
            "more room at the site in a\nBinding protomer than in an\n"
            "Access one. The ordering is\nBinding > Extrusion > Access\n"
            "in every structure.\n\n"
            "Our six protomers all fall\ninside the published spread\n"
            f"for their state (largest\ndeparture {worst:.1f} SD).",
            ha="left", va="top", fontsize=14, color=INK2,
            transform=tx.transAxes, linespacing=1.55)

    fig.text(0.095, 0.105,
             "One sphere fixed in the frame of the reference Binding "
             "protomer; every protomer superposed on its 39 pocket-lining "
             "C\u03b1, so the sphere sits at the\nsame anatomical position "
             "throughout. In an Access or Extrusion protomer that is not "
             "the protomer's own pocket as it would be defined in "
             "isolation - it\nis how open the substrate site is at the "
             "same place. Every MexB-containing entry in the PDB except "
             "22XK and 22XM, whose protomers fail the\nnumbering check as "
             "engineered chimeras at ~40% identity; four entries carry two "
             "trimers in the asymmetric unit.",
             fontsize=13.5, color=INK2, va="top", linespacing=1.5)
    save(fig, "P5_protomer_states")


# ------------------------------------------------------------------- P6
def panel_path_occupancy():
    """Where every bound ligand sits along the transport path."""
    rows = R("ligand_environment.csv")
    if not rows:
        return
    OURS = ("Amp_MexB_20260826", "MexB_DDM_3_20260730")
    NAME = {"21FP": "Chloramphenicol", "Amp_MexB_20260826": "Ampicillin",
            "2V50": "DDM", "3W9I": "DDM", "21FO": "CYMAL-7",
            "3W9J": "EPI", "6IIA": "LMNG",
            "MexB_DDM_3_20260730": "DDM \u00d73"}
    SITECOL = {"DBP": APOLAR, "PBP": POLAR, "both": "#7a8891",
               "neither": "#b9c3c8"}
    ONE = {"PHE": "F", "TYR": "Y", "TRP": "W", "ILE": "I", "VAL": "V",
           "ALA": "A", "ASN": "N", "GLU": "E", "SER": "S", "ARG": "R",
           "GLN": "Q", "LYS": "K", "MET": "M", "LEU": "L", "THR": "T",
           "GLY": "G", "PRO": "P", "ASP": "D", "HIS": "H", "CYS": "C"}

    prot = {}
    for r in rows:
        prot.setdefault((r["pdb"], r["chain"]), []).append(r)
    for v in prot.values():
        v.sort(key=lambda r: float(r["depth_from_entrance_A"]))
    order = sorted(prot, key=lambda k: (
        -(float(prot[k][-1]["depth_from_entrance_A"])
          - float(prot[k][0]["depth_from_entrance_A"])),
        -float(prot[k][-1]["depth_from_entrance_A"])))

    fig = plt.figure(figsize=(10.6, 7.0))
    title(fig, "Three ligands, three stations of one pathway",
          "Every bound ligand in every MexB structure, placed on the same "
          "entry channel.")
    ax = fig.add_axes([0.245, 0.235, 0.50, 0.47])

    for i, k in enumerate(order):
        y = len(order) - 1 - i
        grp = prot[k]
        d = [float(r["depth_from_entrance_A"]) for r in grp]
        mine = k[0] in OURS
        if len(grp) > 1:
            ax.plot([min(d), max(d)], [y, y], color=BINDING, linewidth=5,
                    alpha=.30, solid_capstyle="round", zorder=2)
        for r, x in zip(grp, d):
            ax.scatter([x], [y], s=80 + 2.6 * int(r["heavy_atoms"]),
                       color=SITECOL.get(r["site"], "#b9c3c8"), zorder=4,
                       edgecolor=INK if mine else "white",
                       linewidth=2.0 if mine else 1.6)
        lab = f"{NAME.get(k[0], k[0])}  {k[1]}"
        ax.annotate(lab, (0, y), xycoords=("axes fraction", "data"),
                    xytext=(-12, -5), textcoords="offset points", ha="right",
                    fontsize=14, color=BINDING if mine else INK2,
                    fontweight="bold" if mine else "normal")

    ax.set_yticks([]); ax.set_ylim(-0.8, len(order) - 0.2)
    ax.set_xlim(25, 70)
    ax.set_xlabel("Depth into the porter domain (\u00c5 from the "
                  "periplasmic entrance)", labelpad=12)
    ax.xaxis.label.set_size(16)
    ax.grid(axis="y", visible=False); ax.set_axisbelow(True)
    for sp in ("left",):
        ax.spines[sp].set_visible(False)

    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=12,
                          markerfacecolor=SITECOL[s], markeredgecolor="white",
                          markeredgewidth=1.6,
                          label={"DBP": "distal pocket",
                                 "PBP": "proximal pocket",
                                 "both": "spans both"}[s])
               for s in ("DBP", "PBP", "both")]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(0.50, 0.775), fontsize=15,
               handletextpad=0.35, columnspacing=1.8)

    multi = [r for r in rows if int(r["ligands_in_protomer"]) > 1]
    span = 0.0
    if multi:
        dd = [float(r["depth_from_entrance_A"]) for r in multi]
        span = max(dd) - min(dd)
    tx = fig.add_axes([0.775, 0.235, 0.215, 0.47]); tx.axis("off")
    tx.text(0, 1.0, f"{span:.0f} \u00c5", ha="left", va="top", fontsize=44,
            fontweight="bold", color=BINDING, transform=tx.transAxes)
    tx.text(0, 0.845,
            "of the transport path is\noccupied at once in the\n"
            "DDM \u00d73 structure. Every\nother MexB structure,\n"
            "published or ours, holds\none ligand at one point.\n\n"
            "The three engage ten\naromatic residues between\n"
            "them, but largely different\nones: only F615, F617 and\n"
            "F628 are shared by any\ntwo, and none by all three.",
            ha="left", va="top", fontsize=13.5, color=INK2,
            transform=tx.transAxes, linespacing=1.5)

    fig.text(0.055, 0.10,
             "Marker area tracks ligand size; pocket assignment is by which "
             "lining residues each ligand actually contacts at 4.5 \u00c5, "
             "not by depth - the proximal\nand distal pockets both lie at "
             "the far end of this channel, so depth alone does not separate "
             "them. Ligands are scored one at a time, which is what makes\n"
             "the multi-ligand protomer comparable with the single-ligand "
             "ones. Pocket labels follow this project's residue lists; other "
             "papers partition\nthe same contacts differently. Analyses "
             "after Lawrence et al., Nat Commun 2025;16:10601.",
             fontsize=13.5, color=INK2, va="top", linespacing=1.5)
    save(fig, "P6_path_occupancy")


# ------------------------------------------------------------------- P7
def panel_rotamers():
    """chi1 of the pocket aromatics in every protomer."""
    rows = R("aromatic_rotamers.csv")
    if not rows:
        return
    OURS = ("Amp_MexB_20260826", "MexB_DDM_3_20260730")
    rows = [r for r in rows if r["chi1_deg"] and r["state_call"] in STATE_COLOR]
    if not rows:
        return
    res = sorted({int(r["resseq"]) for r in rows})
    nm = {int(r["resseq"]): r["resname"] for r in rows}

    # deviation of the multi-ligand protomer from the published Binding
    # spread, per residue, so the callout cannot go stale
    def circmean(a):
        a = np.radians(np.asarray(a, float))
        return float(np.degrees(np.arctan2(np.sin(a).mean(),
                                           np.cos(a).mean())))

    devs = []
    for rid in res:
        peers = [float(r["chi1_deg"]) for r in rows
                 if int(r["resseq"]) == rid and r["state_call"] == "Binding"
                 and r["pdb"] not in OURS]
        mine = [float(r["chi1_deg"]) for r in rows
                if int(r["resseq"]) == rid and r["chain"] == "E"
                and r["pdb"] == "MexB_DDM_3_20260730"]
        if len(peers) < 4 or not mine:
            continue
        devs.append((abs((mine[0] - circmean(peers) + 180) % 360 - 180),
                     f"{nm[rid]}{rid}".replace("PHE", "F").replace("TYR", "Y")))
    devs.sort()
    worst = devs[-1] if devs else (0.0, "none")
    rest = devs[-2][0] if len(devs) > 1 else 0.0

    fig = plt.figure(figsize=(10.6, 7.2))
    title(fig, "The pocket lining does not rearrange either",
          f"\u03c71 of every pocket aromatic, in all "
          f"{len({(r['pdb'], r['chain']) for r in rows})} protomers of "
          f"{len({r['pdb'] for r in rows})} structures.")
    ax = fig.add_axes([0.135, 0.225, 0.60, 0.50])

    for i, rid in enumerate(res):
        y = len(res) - 1 - i
        if i % 2 == 0:
            ax.axhspan(y - .5, y + .5, color="#f4f7f8", zorder=0)
        grp = [r for r in rows if int(r["resseq"]) == rid]
        for r in grp:
            x = float(r["chi1_deg"])
            mine = r["pdb"] in OURS
            col = STATE_COLOR[r["state_call"]]
            if mine:
                ax.scatter([x], [y], s=130, color=col, marker="D", zorder=5,
                           edgecolor=INK, linewidth=1.5)
            else:
                ax.scatter([x], [y], s=70, color=col, alpha=.6, zorder=3,
                           edgecolor="white", linewidth=1.1)
    ax.set_yticks(range(len(res)))
    ax.set_yticklabels([f"{nm[r]}{r}" for r in reversed(res)], fontsize=15)
    ax.set_ylim(-.5, len(res) - .5)
    ax.set_xlim(-185, 185)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_xlabel("\u03c71 (degrees)", labelpad=10)
    ax.grid(axis="y", visible=False); ax.set_axisbelow(True)

    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=11,
                          markerfacecolor=c, markeredgecolor="white",
                          markeredgewidth=1.4, label=s_)
               for s_, c in STATE_COLOR.items()]
    handles.append(plt.Line2D([], [], marker="D", linestyle="", markersize=11,
                              markerfacecolor="#cfd8dc", markeredgecolor=INK,
                              markeredgewidth=1.5, label="this work"))
    fig.legend(handles=handles, loc="upper center", ncol=4,
               bbox_to_anchor=(0.43, 0.815), fontsize=15,
               handletextpad=0.35, columnspacing=1.5)

    tx = fig.add_axes([0.765, 0.225, 0.225, 0.50]); tx.axis("off")
    tx.text(0, 1.0, f"{len(devs) - 1} of {len(devs)}", ha="left", va="top",
            fontsize=40, fontweight="bold", color=BINDING,
            transform=tx.transAxes)
    tx.text(0, 0.855,
            "pocket aromatics in the\nDDM \u00d73 protomer sit within\n"
            f"{rest:.0f}\u00b0 of the published\nBinding mean. Three "
            "ligands\nat once do not rotate the\nlining.\n\n"
            f"The exception, {worst[1]}, is\n{worst[0]:.0f}\u00b0 out - "
            "but the ampicillin\nprotomer does the same\nthing, so it "
            "tracks our data,\nnot the ligand count, and\nneeds a density "
            "check.",
            ha="left", va="top", fontsize=13.5, color=INK2,
            transform=tx.transAxes, linespacing=1.45)

    fig.text(0.055, 0.095,
             "\u03c71 = N-CA-CB-CG, straight from the deposited "
             "coordinates; points near -180 and +180 are the same rotamer, "
             "split by the wrap-around. Residue 626 is a\nmethionine in "
             "MexB, so the F626 of Lawrence et al. has no counterpart here. "
             f"{worst[1]} is otherwise strictly state-coupled - t in every "
             "published Access\nprotomer, g+ in every Binding and Extrusion "
             "one. Pocket volume (P4, P5) and lining rotamers together: the "
             "site neither resizes nor rearranges for a\ndifferent or a "
             "larger ligand. What changes is the conformational state, and "
             "which stations along the path are occupied (P6).",
             fontsize=13.5, color=INK2, va="top", linespacing=1.5)
    save(fig, "P7_aromatic_rotamers")


# ------------------------------------------------------------------- P8
def panel_conservation():
    """Is the pocket lining conserved scaffold or specificity tuning?"""
    rows = R("lining_conservation.csv")
    if not rows:
        return
    key = [k for k in rows[0] if k.startswith("homologues_")
           and k != "homologues_aromatic"][0]
    NAMES = ["MexD", "MexF", "MexY", "AcrB", "AcrF", "MdtF", "AcrD"]
    AMINO = {"MexY", "AcrD"}          # aminoglycoside-preferring
    AROM = set("FYWH")

    fig = plt.figure(figsize=(10.6, 6.8))
    title(fig, "The aromatic core is family scaffold; the outer rim is not",
          "Pocket-lining residues of MexB against seven RND transporters, "
          "aligned to MexB.")

    # --- left: aromatic vs non-aromatic identity
    ax = fig.add_axes([0.095, 0.255, 0.30, 0.45])
    groups = [("aromatic", [r for r in rows if r["mexb_aromatic"] == "yes"],
               APOLAR),
              ("everything\nelse", [r for r in rows
                                    if r["mexb_aromatic"] == "no"], POLAR)]
    for i, (lab, grp, col) in enumerate(groups):
        v = np.array([float(r["percent_identical"]) for r in grp])
        x = np.full(len(v), i) + np.linspace(-.19, .19, len(v))
        ax.scatter(x, v, s=110, color=col, zorder=3, edgecolor="white",
                   linewidth=1.6, alpha=.85)
        m = float(v.mean())
        ax.plot([i - .32, i + .32], [m, m], color=INK, linewidth=4,
                solid_capstyle="round", zorder=4)
        ax.annotate(f"{m:.0f}%", (i, m), textcoords="offset points",
                    xytext=(0, 14), ha="center", fontsize=23,
                    fontweight="bold", color=col)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([g[0] for g in groups], fontsize=17)
    ax.set_xlim(-.6, 1.6); ax.margins(y=.22)
    ax.set_ylabel("identity to MexB across the\nseven transporters (%)",
                  labelpad=8, fontsize=16)
    ax.grid(axis="x", visible=False); ax.set_axisbelow(True)

    # --- right: aromatic scaffold retained, per transporter
    ax2 = fig.add_axes([0.505, 0.255, 0.44, 0.45])
    arom_rows = [r for r in rows if r["mexb_aromatic"] == "yes"]
    n = len(arom_rows)
    keep = []
    for k, nm in enumerate(NAMES):
        kept = sum(1 for r in arom_rows if r[key][k] in AROM)
        keep.append((nm, kept))
    keep.sort(key=lambda x: -x[1])
    ys = np.arange(len(keep))[::-1]
    for (nm, kept), y in zip(keep, ys):
        col = WARN if nm in AMINO else TEAL
        ax2.barh([y], [kept], color=col, height=.62, zorder=3,
                 alpha=.9 if nm in AMINO else .75)
        ax2.annotate(f"{kept}/{n}", (kept, y), textcoords="offset points",
                     xytext=(8, -6), fontsize=15, color=col,
                     fontweight="bold")
    ax2.set_yticks(ys)
    ax2.set_yticklabels([nm + (" *" if nm in AMINO else "")
                         for nm, _ in keep], fontsize=15)
    ax2.set_xlim(0, n + 1.6)
    ax2.set_xlabel(f"MexB pocket aromatics still aromatic "
                   f"(of {n})", fontsize=16)
    ax2.grid(axis="y", visible=False); ax2.set_axisbelow(True)

    look = {int(r["resseq"]): float(r["percent_identical"]) for r in rows}
    st = {k: np.mean([look[i] for i in v if i in look]) for k, v in
          (("out", [136, 573, 617, 628, 664, 666, 327]),
           ("mid", [615, 617]), ("deep", [178, 610, 615, 628]))}
    fig.text(0.095, 0.135,
             "*  MexY and AcrD prefer aminoglycosides - polar, cationic "
             "substrates MexB exports poorly - and they are the two that "
             "lose the aromatic core.\nAcross the three stations the DDM "
             "\u00d73 ligands occupy, conservation falls as they approach "
             f"the entrance: deepest {st['deep']:.0f}%, middle "
             f"{st['mid']:.0f}%, outermost {st['out']:.0f}%. The deep\n"
             "station is built from family-invariant residues (Y327, "
             "F628 identical in all seven; F610, F615 in six); the outer "
             "one is not (F573 identical in none).",
             fontsize=13.5, color=INK2, va="top", linespacing=1.5)
    save(fig, "P8_lining_conservation")


# ------------------------------------------------------------------- P9
def panel_pocket_physchem():
    """Lipophilicity and electrostatics of the site, by state."""
    pk = R("protomer_pockets.csv")
    el = R("pocket_electrostatics.csv")
    if not pk or not el:
        return
    STATES = ("Access", "Binding", "Extrusion")

    def col(rows, key, st):
        out = []
        for r in rows:
            if r["state_call"] != st:
                continue
            try:
                out.append(float(r[key]))
            except (KeyError, TypeError, ValueError):
                pass
        return np.array(out)

    fig = plt.figure(figsize=(10.6, 6.6))
    title(fig, "The site changes chemistry around the cycle, not just size",
          "Lipophilicity and electrostatic potential of the same substrate "
          "site, measured per protomer.")
    specs = [(pk, "lipophilic_index_pct",
              "lipophilic index of the\nsite surface (% C or S)", "{:.0f}%"),
             (el, "mean_potential_kT_e",
              "mean electrostatic\npotential (kT/e)", "{:+.1f}")]
    gs = fig.add_gridspec(1, 2, left=0.105, right=0.985, top=0.60,
                          bottom=0.235, wspace=0.36)
    for k, (rows, key, ylab, f) in enumerate(specs):
        ax = fig.add_subplot(gs[0, k])
        for i, st in enumerate(STATES):
            v = col(rows, key, st)
            if not len(v):
                continue
            x = np.full(len(v), i) + np.linspace(-.18, .18, len(v))
            ax.scatter(x, v, s=110, color=STATE_COLOR[st], zorder=3,
                       edgecolor="white", linewidth=1.5, alpha=.8)
            m = float(v.mean())
            ax.plot([i - .32, i + .32], [m, m], color=INK, linewidth=4,
                    solid_capstyle="round", zorder=4)
            ax.annotate(f.format(m), (i - .34, m),
                        textcoords="offset points", xytext=(-6, -8),
                        ha="right", fontsize=20, fontweight="bold",
                        color=STATE_COLOR[st])
        if key.startswith("mean_potential"):
            ax.axhline(0, color="#9fb0b8", linewidth=1.6,
                       linestyle=(0, (5, 4)))
        ax.set_xticks(range(len(STATES)))
        ax.set_xticklabels(STATES, fontsize=17)
        ax.set_xlim(-.95, len(STATES) - .4)
        ax.set_ylabel(ylab, labelpad=8, fontsize=16)
        ax.margins(y=.24)
        ax.grid(axis="x", visible=False); ax.set_axisbelow(True)

    fig.text(0.105, 0.135,
             "The site is electronegative in every protomer of every "
             "structure - it never once comes out positive - and most so in "
             "Access, the state that faces the\nperiplasmic entrance. It is "
             "at its most hydrophobic in Binding and its least in Extrusion, "
             "where the substrate is released. Lipophilic index over 48 "
             "protomers;\npotential over 34, one per state per structure "
             "plus all six of ours (PDB2PQR AMBER charges, APBS linearised "
             "PB, pdie 2 / sdie 78.54, 298 K). Both are\nsampled over the "
             "free points of the same 16 \u00c5 sphere used for volume. "
             "Spread within a state is wide, so read the ordering, not the "
             "individual values.",
             fontsize=13.5, color=INK2, va="top", linespacing=1.5)
    save(fig, "P9_pocket_physchem")


# ------------------------------------------------------------------ P10
def panel_mechanism():
    """AcrB's drugs and this structure's three DDM, on one measured axis.

    AcrB is 70% identical to MexB, so its pocket-lining residues map onto
    MexB's by alignment and its protomers superpose on the same reference
    used throughout. Its bound drugs therefore carry the same depth
    coordinate as the MexB ligands - nothing here is schematic.
    """
    env = R("ligand_environment.csv")
    acr = R("acrb_ligands.csv")
    ml = R("multiligand_survey.csv")
    multi = sorted(((float(r["depth_from_entrance_A"]), r) for r in env
                    if int(r["ligands_in_protomer"]) > 1),
                   key=lambda x: x[0])
    ours = next((r for r in ml if r["structure"].startswith("MexB_DDM")), None)
    if len(multi) < 3 or not acr or ours is None:
        return
    touch = float(ours["closest_approach_A"])
    span = multi[-1][0] - multi[0][0]
    SITECOL = {"DBP": APOLAR, "PBP": POLAR, "both": "#7a8891",
               "outside": "#b9c3c8"}

    # the most drugs AcrB puts in one protomer, and how far apart they are
    per = {}
    for r in acr:
        per.setdefault((r["pdb"], r["chain"]), []).append(r)
    biggest = max(per.values(), key=len)
    n_max = len(biggest)
    pair_gap = (abs(float(biggest[0]["depth_from_entrance_A"])
                    - float(biggest[-1]["depth_from_entrance_A"]))
                if n_max > 1 else 0.0)

    fig = plt.figure(figsize=(10.6, 9.0))
    title(fig, "One substrate at a time, or a chain across the path?",
          "Every AcrB drug and our three DDM, superposed on one reference "
          "and measured on the same channel.")
    callout(fig, 0.055, 0.875, f"{span:.0f} \u00c5",
            "of the path occupied at once, as one\ncontiguous chain of "
            "three molecules", BINDING, size=38)
    callout(fig, 0.535, 0.875, f"{n_max} at most",
            f"drugs in any one AcrB protomer, and\nthose two "
            f"{pair_gap:.0f} \u00c5 apart in the same pocket", ACCESS,
            size=38)

    ax = fig.add_axes([0.075, 0.245, 0.90, 0.46])
    ax.set_xlim(-1.5, 66); ax.set_ylim(-1.30, 4.15)
    ax.set_yticks([])
    ax.grid(axis="y", visible=False); ax.set_axisbelow(True)
    for sp in ("left",):
        ax.spines[sp].set_visible(False)

    # the channel itself, drawn at its measured radius rather than sketched
    prof = tunnel_profile()
    KY = 0.15                      # plot units per Angstrom of radius
    if prof is not None:
        dep, rad = prof
        for y0 in (2.30, 0.40):
            ax.fill_between(dep, y0 - KY * rad, y0 + KY * rad,
                            color="#eef4f6", zorder=0, linewidth=0)
            for sgn in (1, -1):
                ax.plot(dep, y0 + sgn * KY * rad, color="#bdccd2",
                        linewidth=1.5, zorder=1)
        # a scale bar, because the vertical axis is otherwise unitless
        ax.plot([1.2, 1.2], [-0.98 - KY * 4, -0.98 + KY * 4], color=INK2,
                linewidth=2.6, solid_capstyle="butt")
        ax.annotate("8 \u00c5 across", (1.2, -0.98),
                    textcoords="offset points", xytext=(9, -5), ha="left",
                    fontsize=12.5, color=INK2)
        ax.annotate(f"radius {rad.min():.1f}\u2013{rad.max():.1f} \u00c5; "
                    f"vertical scale is not the horizontal one",
                    (65.4, -0.98), ha="right", va="center", fontsize=12.5,
                    color=INK2, fontstyle="italic")

    # ---- top row: AcrB drugs, one point per bound copy
    rng = np.random.default_rng(1)
    for r in acr:
        x = float(r["depth_from_entrance_A"])
        y = 2.30 + rng.uniform(-.16, .16)
        ax.scatter([x], [y], s=150, color=SITECOL.get(r["site"], "#b9c3c8"),
                   zorder=4, edgecolor="white", linewidth=1.7, alpha=.9)
    # two label rows: doxorubicin at 28 and rifampicin at 32 are too close
    # to sit on one
    LAB = {"doxorubicin \u00d72": (31.0, 2.30, 0, 48, "center"),
           "minocycline": (62.8, 2.30, 0, 48, "center"),
           "rifampicin": (32.1, 2.30, 0, 76, "center"),
           "erythromycin": (40.6, 2.30, 0, 76, "center"),
           "MBX inhibitor": (55.9, 2.30, 0, 76, "center")}
    for nm, (x, y, dx, dy, ha) in LAB.items():
        ax.annotate(nm, (x, y), textcoords="offset points", xytext=(dx, dy),
                    ha=ha, fontsize=13, color=INK2)
    ax.text(-1.0, 3.98, f"AcrB \u2014 {len(acr)} bound drugs, "
            f"{len(per)} protomers, {len({r['pdb'] for r in acr})} structures",
            fontsize=17, fontweight="bold", color=ACCESS, va="center")

    # ---- bottom row: the three DDM, in contact
    xs = [d for d, _ in multi]
    ax.plot([min(xs), max(xs)], [0.40, 0.40], color=BINDING, linewidth=7,
            alpha=.30, solid_capstyle="round", zorder=2)
    for d, r in multi:
        ax.scatter([d], [0.40], s=330, color=SITECOL.get(r["site"], "#7a8891"),
                   zorder=5, marker="D", edgecolor=BINDING, linewidth=2.6)
        ax.annotate(f"{d:.0f} \u00c5", (d, 0.40), textcoords="offset points",
                    xytext=(0, 22), ha="center", fontsize=13.5,
                    color=BINDING, fontweight="bold")
    ax.text(-1.0, 1.32, "MexB with three DDM bound (this work), one protomer",
            fontsize=17, fontweight="bold", color=BINDING, va="center")
    ax.annotate(f"in van der Waals contact, closest approach "
                f"{touch:.1f} \u00c5", (float(np.mean(xs)), 0.40),
                textcoords="offset points", xytext=(0, -46), ha="center",
                fontsize=14, color=BINDING, fontweight="bold")

    ax.set_xlabel("Depth into the porter domain (\u00c5 from the "
                  "periplasmic entrance)", labelpad=12)
    ax.xaxis.label.set_size(16)
    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=12,
                          markerfacecolor=SITECOL[k], markeredgecolor="white",
                          markeredgewidth=1.6,
                          label={"DBP": "distal pocket",
                                 "PBP": "proximal pocket",
                                 "both": "spans both"}[k])
               for k in ("DBP", "PBP", "both")]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(0.53, 0.755), fontsize=15,
               handletextpad=0.35, columnspacing=1.8)

    fig.text(0.055, 0.115,
             "AcrB drugs span the same depths as the MexB ligands, so both "
             "transporters use the whole path - but no AcrB protomer holds "
             "more than two drugs, and\nthe one that does (4DX7 chain A, two "
             "doxorubicin) has them stacked in the same pocket, not spread "
             "along the route. Access-state protomers superpose on the "
             "binding-state\nreference at 2.9-3.0 \u00c5 against 0.8-1.9 "
             "\u00c5 for binding-state ones, so their depths are the coarser "
             "numbers here. DDM is a detergent: this shows the path can be "
             "occupied at\nthree stations at once, not that substrates are "
             "carried as a chain. Pocket assignment is by contact in the "
             "common frame; depth does not separate the two pockets.",
             fontsize=13.5, color=INK2, va="top", linespacing=1.5)
    save(fig, "P10_mechanism_contrast")


# ------------------------------------------------------------------ P11
def panel_mexb_rows():
    """Every MexB ligand-bound protomer on its own row of the same channel."""
    env = R("ligand_environment.csv")
    if not env:
        return
    OURS = ("Amp_MexB_20260826", "MexB_DDM_3_20260730")
    NAME = {"21FP": "Chloramphenicol", "Amp_MexB_20260826": "Ampicillin",
            "2V50": "DDM", "3W9I": "DDM", "21FO": "CYMAL-7",
            "3W9J": "EPI", "6IIA": "LMNG",
            "MexB_DDM_3_20260730": "DDM \u00d73"}
    SITECOL = {"DBP": APOLAR, "PBP": POLAR, "both": "#7a8891",
               "neither": "#b9c3c8"}
    ONE = {"PHE": "F", "TYR": "Y", "TRP": "W", "ILE": "I", "VAL": "V",
           "ALA": "A", "ASN": "N", "GLU": "E", "SER": "S", "ARG": "R",
           "GLN": "Q", "LYS": "K", "MET": "M", "LEU": "L", "THR": "T",
           "GLY": "G", "PRO": "P", "ASP": "D", "HIS": "H", "CYS": "C"}

    prot = {}
    for r in env:
        prot.setdefault((r["pdb"], r["chain"]), []).append(r)
    for v in prot.values():
        v.sort(key=lambda r: float(r["depth_from_entrance_A"]))

    # one protomer per ligand: several structures contribute two copies of
    # the same ligand and 2V50 and 3W9I both contribute DDM. Keep the
    # best-defined site - the copy contacting the most lining residues,
    # then the larger one - so the row set is one row per chemistry.
    best = {}
    for k, v in prot.items():
        nm = NAME.get(k[0], k[0])
        score = (sum(int(r["residues_contacted"]) for r in v),
                 sum(int(r["heavy_atoms"]) for r in v))
        if nm not in best or score > best[nm][0]:
            best[nm] = (score, k)
    prot = {k: prot[k] for _, k in best.values()}
    # ours first, then by how much of the path each one covers
    order = sorted(prot, key=lambda k: (
        k[0] not in OURS,
        -(float(prot[k][-1]["depth_from_entrance_A"])
          - float(prot[k][0]["depth_from_entrance_A"])),
        -float(prot[k][-1]["depth_from_entrance_A"])))

    n = len(order)
    H = 4.6 + 0.52 * n
    fig = plt.figure(figsize=(10.6, H))
    title(fig, "One channel, every substrate-bound MexB protomer",
          "Each structure's own tunnel, traced separately and projected onto "
          "one shared axis, with its ligands placed on it.")

    # the two numbers this panel exists to make, in P10's callout style
    multi = [r for r in env if int(r["ligands_in_protomer"]) > 1]
    span = 0.0
    if multi:
        dd = [float(r["depth_from_entrance_A"]) for r in multi]
        span = max(dd) - min(dd)
    bn = [float(r["route_bottleneck_A"])
          for r in R("per_structure_tunnel_summary.csv")
          if r.get("route_bottleneck_A")]
    callout(fig, 0.055, 1.0 - 1.02 / H, f"{span:.0f} \u00c5",
            "of the path occupied at once, by the\nthree DDM of one "
            "protomer", BINDING, size=38)
    if bn:
        callout(fig, 0.535, 1.0 - 1.02 / H,
                f"{min(bn):.1f}\u2013{max(bn):.1f} \u00c5",
                "bottleneck across the seven tunnels,\neach traced from its "
                "own structure", TEAL, size=38)

    y0, htop = 1.80 / H, 3.95 / H
    ax = fig.add_axes([0.245, y0, 0.660, 1.0 - y0 - htop])
    ax.set_xlim(-1.5, 66); ax.set_ylim(-1.85, n + 0.10)
    ax.annotate("Bottleneck", (1.0, n - 0.35),
                xycoords=("axes fraction", "data"), xytext=(8, -4),
                textcoords="offset points", ha="left", fontsize=12,
                color=INK2, annotation_clip=False)
    ax.set_yticks([]); ax.grid(axis="y", visible=False)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)

    prof = tunnel_profile()
    KY = 0.10
    dep = rad = None
    if prof is not None:
        dep, rad = prof

    # each structure's own tunnel, traced separately and re-expressed on the
    # reference axis. Coverage is partial for every one but the reference
    # itself, so the reference outline stays behind each row as the frame and
    # the row's own profile is drawn only where its tunnel actually goes.
    own = {}
    for r in R("per_structure_tunnels.csv"):
        own.setdefault((r["pdb"], r["chain"]), []).append(
            (float(r["depth_from_entrance_A"]), float(r["radius_min_A"])))
    cover = {(r["pdb"], r["chain"]): r
             for r in R("per_structure_tunnel_summary.csv")}

    for i, k in enumerate(order):
        y = n - 1 - i
        grp = prot[k]
        mine = k[0] in OURS
        lc = LIGCOL.get(k[0], TEAL)
        if dep is not None:                      # reference channel, faint
            for sgn in (1, -1):
                ax.plot(dep, y + sgn * KY * rad, color="#d3dde1",
                        linewidth=1.1, zorder=0)
        pts = sorted(own.get(k, []))
        if pts:                                  # this structure's own tunnel
            od = np.array([p[0] for p in pts])
            orr = np.array([p[1] for p in pts])
            # break the fill wherever the projection has a gap
            cut = np.where(np.diff(od) > 2.5)[0] + 1
            for seg_d, seg_r in zip(np.split(od, cut), np.split(orr, cut)):
                if len(seg_d) < 2:
                    continue
                ax.fill_between(seg_d, y - KY * seg_r, y + KY * seg_r,
                                color=tint(lc, 0.85), zorder=1, linewidth=0)
                for sgn in (1, -1):
                    ax.plot(seg_d, y + sgn * KY * seg_r, color=tint(lc, 0.30),
                            linewidth=1.6, zorder=2)
            c = cover.get(k)
            if c and c.get("route_bottleneck_A"):
                ax.annotate(f"{float(c['route_bottleneck_A']):.1f} \u00c5",
                            (1.0, y), xycoords=("axes fraction", "data"),
                            xytext=(8, -4), textcoords="offset points",
                            ha="left", fontsize=12, color=lc,
                            annotation_clip=False)
        d = [float(r["depth_from_entrance_A"]) for r in grp]
        if len(grp) > 1:
            ax.plot([min(d), max(d)], [y, y], color=lc, linewidth=6,
                    alpha=.35, solid_capstyle="round", zorder=2)
        for r, x in zip(grp, d):
            ax.scatter([x], [y], s=70 + 2.2 * int(r["heavy_atoms"]),
                       color=SITECOL.get(r["site"], "#b9c3c8"), zorder=4,
                       edgecolor=lc, linewidth=2.2)
        lab = NAME.get(k[0], k[0])
        ax.annotate(lab, (0, y), xycoords=("axes fraction", "data"),
                    xytext=(-12, -5), textcoords="offset points", ha="right",
                    fontsize=13.5, color=lc,
                    fontweight="bold" if mine else "normal")

    ax.set_xlabel("Depth into the porter domain (\u00c5 from the "
                  "periplasmic entrance)", labelpad=10)
    ax.xaxis.label.set_size(16)

    # a residue scale above the axis. Only side chains that actually lie
    # near the centreline can be placed on it, and the deep end is crowded
    # enough that the labels alternate between two heights.
    res = R("channel_residues.csv")
    look = {int(r["resseq"]): r for r in res}
    ROWS = [[617, 79, 628, 178], [676, 615], [620]]
    for lvl, pick in enumerate(ROWS):
        for rid in pick:
            r = look.get(rid)
            if r is None or float(r["offset_from_channel_A"]) > 8.0:
                continue
            x = float(r["depth_from_entrance_A"])
            site = r["site"].split(";")[0]
            col = WARN if site == "switch" else SITECOL.get(site, "#7a8891")
            ax.annotate("", (x, 1.0), xycoords=("data", "axes fraction"),
                        xytext=(x, 1.0 + 0.012 + 0.030 * lvl),
                        textcoords=("data", "axes fraction"),
                        arrowprops=dict(arrowstyle="-", color=col, lw=1.6),
                        annotation_clip=False)
            ax.annotate(f"{ONE.get(r['resname'], r['resname'])}{rid}",
                        (x, 1.0), xycoords=("data", "axes fraction"),
                        xytext=(0, 11 + 22 * lvl), textcoords="offset points",
                        ha="center", fontsize=11.5, fontweight="bold",
                        color=col, annotation_clip=False)
    ax.annotate("Pocket-lining residues, coloured by which pocket they "
                "belong to", (0.5, 1.0), xycoords="axes fraction",
                xytext=(0, 68), textcoords="offset points", ha="center",
                fontsize=14, color=INK2, annotation_clip=False)

    if rad is not None:
        ax.plot([1.2, 1.2], [-1.42 - KY * 4, -1.42 + KY * 4], color=INK2,
                linewidth=2.4, solid_capstyle="butt")
        ax.annotate("8 \u00c5 across", (1.2, -1.42),
                    textcoords="offset points", xytext=(9, -5), ha="left",
                    fontsize=12, color=INK2)

    ax.annotate("Entry cleft", (16.5, -1.42), ha="center", va="center",
                fontsize=12.5, color=INK2, fontstyle="italic")
    ax.annotate("Porter pocket", (52.0, -1.42), ha="center", va="center",
                fontsize=12.5, color=INK2, fontstyle="italic")

    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=11,
                          markerfacecolor=SITECOL[q], markeredgecolor="white",
                          markeredgewidth=1.5,
                          label={"DBP": "distal pocket",
                                 "PBP": "proximal pocket",
                                 "both": "spans both"}[q])
               for q in ("DBP", "PBP", "both")]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(0.62, 1.0 - 2.15 / H), fontsize=14,
               handletextpad=0.35, columnspacing=1.6)

    fig.text(0.045, 0.085,
             "Each row is that structure's own tunnel, traced from its own "
             "coordinates and re-expressed on the reference axis.\n"
             "The faint grey outline behind every row is the reference "
             "channel, the widest ligand-free route out of ampicillin\n"
             "chain E. A tunnel is drawn only where it follows that route "
             "within 6 \u00c5; coverage runs 22\u2013100%, so a gap is "
             "where a\nstructure's widest route goes elsewhere, not where "
             "it is blocked. Tube half-width is the local radius, on a "
             "vertical\nscale that is not the horizontal one. The "
             "bottleneck at the right is the narrowest point on the route "
             "itself, taken\nwith the terminal 3 \u00c5 at each end "
             "trimmed off: the trace is seeded beside the ligand, so its "
             "deep cap measures the\nclearance of the pocket the ligand "
             "sits in rather than any constriction the route passes "
             "through. Both numbers\nare tabulated "
             "(per_structure_tunnel_summary.csv); they differ by up to "
             "1.0 \u00c5.\n\n"
             "Marker area tracks ligand size and its fill gives which "
             "lining set the ligand contacts at 4.5 \u00c5. One row per "
             "ligand:\nwhere a structure or a pair of structures gave more "
             "than one copy, the copy contacting the most lining residues "
             "is\nkept. Only the DDM \u00d73 protomer carries more than "
             "one ligand at once. The two pockets are not two stretches of "
             "this\naxis \u2014 both have residues at 26\u201335 and "
             "again at 62\u201363 \u00c5, since depth is arc length along "
             "a winding path \u2014 so they are\nmarked by residue rather "
             "than shaded as bands, with the switch loop in orange. No "
             "lining residue and no modelled\nligand lies shallower than "
             "26 \u00c5; that stretch is the open entry cleft.",
             fontsize=13, color=INK2, va="top", linespacing=1.5)
    save(fig, "P11_mexb_rows")


# ------------------------------------------------------------------ P12
def panel_regional_rmsd():
    """Localised backbone RMSD by region, after Lawrence et al. Fig. 2B."""
    rows = R("regional_rmsd.csv")
    if not rows:
        return
    ORDER = ["TM1", "PN1", "PN2", "DN", "TM2", "Ialpha", "TM3-6",
             "loop496-515", "TM6b", "PC1", "PC2", "DC", "junction859-875",
             "TM7-12"]
    STATES = ("Access", "Binding", "Extrusion")

    vs_acrb = {}
    for r in rows:
        if r["comparison"] == "MexB (DDM) vs AcrB":
            vs_acrb.setdefault(r["state"], {})[r["region"]] = float(r["rmsd_A"])
    ours = {}
    for r in rows:
        if r["comparison"].startswith("ampicillin vs DDM"):
            ours.setdefault(r["comparison"][-1], {})[r["region"]] = \
                float(r["rmsd_A"])
    if not vs_acrb:
        return
    regions = [g for g in ORDER if any(g in v for v in vs_acrb.values())]

    fig = plt.figure(figsize=(11.4, 7.6))
    title(fig, "Where MexB differs from AcrB, region by region",
          "Backbone C\u03b1 RMSD per region after superposing whole "
          "protomers, for each state of the functional rotation.")
    ax = fig.add_axes([0.075, 0.365, 0.63, 0.375])
    x = np.arange(len(regions))
    w = 0.26
    for k, st in enumerate(STATES):
        v = [vs_acrb.get(st, {}).get(g, np.nan) for g in regions]
        ax.bar(x + (k - 1) * w, v, width=w, color=STATE_COLOR[st],
               alpha=.85, label=st, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(regions, rotation=38, ha="right", fontsize=13)
    ax.set_ylabel("C\u03b1 RMSD to AcrB (\u00c5)", labelpad=8, fontsize=16)
    ax.grid(axis="x", visible=False); ax.set_axisbelow(True)
    ax.legend(loc="upper left", fontsize=14, ncol=3)

    glob = {st: float(next(r["global_rmsd_A"] for r in rows
                           if r["comparison"] == "MexB (DDM) vs AcrB"
                           and r["state"] == st))
            for st in STATES if st in vs_acrb}
    worst = max(((v, g, st) for st, d in vs_acrb.items()
                 for g, v in d.items()), key=lambda t: t[0])

    tx = fig.add_axes([0.765, 0.335, 0.22, 0.405]); tx.axis("off")
    tx.text(0, 1.0, f"{worst[0]:.1f} \u00c5", ha="left", va="top",
            fontsize=42, fontweight="bold", color=STATE_COLOR[worst[2]],
            transform=tx.transAxes)
    body = (f"at {worst[1]}, the largest\nregional difference.\n\n"
            + "Global RMSD to AcrB:\n"
            + "\n".join(f"  {st}  {glob[st]:.2f} \u00c5"
                         for st in STATES if st in glob)
            + "\n\nOur two models, mean\nover regions:\n"
            + "  " + ",  ".join(
                f"{c} {np.nanmean(list(d.values())):.2f}"
                for c, d in sorted(ours.items())) + " \u00c5")
    tx.text(0, 0.80, body, ha="left", va="top", fontsize=13.5, color=INK2,
            transform=tx.transAxes, linespacing=1.5)

    fig.text(0.055, 0.105,
             "Each MexB protomer of the DDM \u00d73 model is superposed on "
             "the AcrB 4DX5 protomer in the same state over all ~1030 shared "
             "C\u03b1, then RMSD is\ntaken per region on that one fit, so a "
             "region's value is how far it sits from where the global "
             "superposition puts it. MexB residues map onto AcrB by the same "
             "global\nBLOSUM62 alignment used for the conservation panel "
             "(1045 of 1046 map). AcrB states are assigned by ranking 4DX5's "
             "own three protomers on their PN1\u2013PN2 and\nPC1\u2013PC2"
             " separations, not by this project's absolute cutoffs, which are "
             "calibrated on MexB and put two of AcrB's chains in the same "
             "state.",
             fontsize=13, color=INK2, va="top", linespacing=1.5)
    save(fig, "P12_regional_rmsd")


# ------------------------------------------------------------------ P13
def panel_tm_overlay():
    """TM cartoon overlays with AcrB, after Lawrence et al. Fig. 2C."""
    import matplotlib.image as mpimg
    src = os.path.join(FIGURES, "tm_overlay")
    STATES = ("Access", "Binding", "Extrusion")
    REPS = ("R1", "R2")
    have = {(st, rp): os.path.join(src, f"{st}_{rp}.png")
            for st in STATES for rp in REPS}
    have = {k: v for k, v in have.items() if os.path.exists(v)}
    if len(have) < 6:
        return
    fit = {r["state"]: r["tm_fit_rmsd_A"] for r in R("tm_overlay.csv")}
    MEXB, ACRB, SWING = "#2E5FE8", "#E59BD8", "#D11149"
    hd = [r for r in R("helix_displacement.csv")
          if "vs Access" in r["comparison"] and r["helix"]]
    top = (max(hd, key=lambda r: float(r["centroid_shift_A"]))
           if hd else None)

    fig = plt.figure(figsize=(10.6, 8.4))
    title(fig, "MexB and AcrB transmembrane domains, state by state",
          "The two pseudo-symmetric repeats superposed and drawn as "
          "cartoons, from one camera, for each protomer state.")
    fig.text(0.055, 0.876, "MexB", fontsize=20, fontweight="bold",
             color=MEXB, va="top")
    fig.text(0.175, 0.876, "AcrB (4DX5)", fontsize=20, fontweight="bold",
             color=ACRB, va="top")
    if top:
        fig.text(0.395, 0.876,
                 f"{top['helix']} \u2014 moves {float(top['centroid_shift_A']):.1f} "
                 f"\u00c5 from access to {top['state'].lower()}",
                 fontsize=18, fontweight="bold", color=SWING, va="top")

    def crop(a):
        """Trim the white margin PyMOL leaves around the cartoon."""
        ink = (a[..., :3] < 0.97).any(-1) if a.ndim == 3 else a < 0.97
        ys, xs = np.where(ink)
        if not len(ys):
            return a
        m = 8
        return a[max(ys.min() - m, 0):ys.max() + m,
                 max(xs.min() - m, 0):xs.max() + m]

    left, right, top, bot = 0.075, 0.985, 0.805, 0.185
    wcell = (right - left) / 3
    hcell = (top - bot) / 2
    for j, st in enumerate(STATES):
        for i, rp in enumerate(REPS):
            ax = fig.add_axes([left + j * wcell, top - (i + 1) * hcell,
                               wcell * 0.97, hcell * 0.94])
            ax.imshow(crop(mpimg.imread(have[(st, rp)])))
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            if i == 0:
                ax.set_title(st, fontsize=19, fontweight="bold",
                             color=STATE_COLOR[st], pad=8)
            if j == 0:
                ax.set_ylabel(rp, fontsize=19, fontweight="bold",
                              color=INK, labelpad=10, rotation=0,
                              va="center")
    for j, st in enumerate(STATES):
        if st in fit:
            fig.text(left + j * wcell + wcell * 0.48, bot - 0.018,
                     f"TM fit {float(fit[st]):.2f} \u00c5", ha="center",
                     va="top", fontsize=14, color=INK2)

    fig.text(0.055, 0.108,
             "R1 is TM1\u20136 with the I\u03b1 helix, R2 is TM7\u201312. "
             "Each MexB protomer of the DDM \u00d73 model is superposed on "
             "the AcrB 4DX5\nprotomer in the same state over their ~390 "
             "shared transmembrane C\u03b1, then that rigid pair is moved "
             "onto the AcrB access\nprotomer, so all six panels share one "
             "camera and differences between them are real. AcrB states are "
             "assigned by ranking\n4DX5's own protomers on their "
             "PN1\u2013PN2 and PC1\u2013PC2 separations. Rendered in "
             "PyMOL as cylindrical helices, with helix\nboundaries from its "
             "secondary-structure assignment rather than hard-coded. The "
             "highlighted helix is the one measured to\nmove furthest "
             "between MexB protomer states, not one picked by eye; every "
             "helix's displacement and axis rotation is in\n"
             "helix_displacement.csv.",
             fontsize=13, color=INK2, va="top", linespacing=1.5)
    save(fig, "P13_tm_overlay")


# ------------------------------------------------------------------ P14
def panel_caver():
    """Our widest-path bottlenecks against CAVER 3.0.3."""
    cav = R("caver_tunnels.csv")
    own = {r["ligand"]: r for r in R("own_axis_summary.csv")}
    if not cav:
        return
    NAME2PDB = {"Ampicillin": "Amp_MexB_20260826",
                "DDM x3": "MexB_DDM_3_20260730", "DDM \u00d73":
                "MexB_DDM_3_20260730", "Chloramphenicol": "21FP",
                "DDM": "2V50", "EPI": "3W9J", "LMNG": "6IIA",
                "CYMAL-7": "21FO"}
    rows = []
    for r in cav:
        nm = r["ligand"]
        # CAVER reports the narrowest point along its route, so compare it
        # with ours measured the same way - the route bottleneck, not the
        # clearance at the seed where the trace starts beside the ligand.
        o = float(r["our_route_bottleneck_A"])
        c = float(r["caver_bottleneck_A"])
        clen = float(r["caver_length_A"])
        olen = (float(own[nm]["tunnel_length_A"])
                if nm in own else float("nan"))
        rows.append((nm, o, c, clen, olen, c - o))
    rows.sort(key=lambda x: abs(x[5]))

    fig = plt.figure(figsize=(11.0, 6.8))
    title(fig, "An independent tool agrees where it measures the same route",
          "CAVER 3.0.3 against this project's widest-path search, on the "
          "same protomers, seeded on the same points.")

    ax = fig.add_axes([0.155, 0.28, 0.42, 0.44])
    for i, (nm, o, c, clen, olen, d) in enumerate(rows):
        y = len(rows) - 1 - i
        col = LIGCOL.get(NAME2PDB.get(nm, ""), TEAL)
        ax.plot([o, c], [y, y], color=col, linewidth=3, alpha=.45, zorder=2,
                solid_capstyle="round")
        ax.scatter([o], [y], s=150, facecolor="white", edgecolor=col,
                   linewidth=2.4, zorder=4)
        ax.scatter([c], [y], s=150, color=col, zorder=4,
                   edgecolor="white", linewidth=1.6)
        ax.annotate(nm, (0, y), xycoords=("axes fraction", "data"),
                    xytext=(-10, -5), textcoords="offset points",
                    ha="right", fontsize=13.5, color=col)
        ax.annotate(f"{d:+.2f}", (1, y), xycoords=("axes fraction", "data"),
                    xytext=(8, -5), textcoords="offset points", ha="left",
                    fontsize=12.5, color=INK2, annotation_clip=False)
    ax.set_yticks([]); ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlim(1.0, 3.05)
    ax.set_xlabel("Route bottleneck radius (\u00c5)", labelpad=10,
                  fontsize=15)
    ax.grid(axis="y", visible=False); ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.annotate("\u0394", (1, len(rows) - 0.35),
                xycoords=("axes fraction", "data"), xytext=(8, -5),
                textcoords="offset points", ha="left", fontsize=13,
                color=INK2, annotation_clip=False)
    h = [plt.Line2D([], [], marker="o", linestyle="", markersize=11,
                    markerfacecolor="white", markeredgecolor=INK2,
                    markeredgewidth=2, label="this work"),
         plt.Line2D([], [], marker="o", linestyle="", markersize=11,
                    markerfacecolor=INK2, markeredgecolor="white",
                    markeredgewidth=1.5, label="CAVER")]
    ax.legend(handles=h, loc="upper left", fontsize=13.5, ncol=1)

    # how much of our route CAVER's best-ranked tunnel actually spans
    ax2 = fig.add_axes([0.735, 0.28, 0.235, 0.44])
    rat = np.array([r[3] / r[4] for r in rows])
    dif = np.array([abs(r[5]) for r in rows])
    for (nm, o, c, clen, olen, d), x in zip(rows, rat):
        col = LIGCOL.get(NAME2PDB.get(nm, ""), TEAL)
        ax2.scatter([100 * x], [abs(d)], s=150, color=col, zorder=4,
                    edgecolor="white", linewidth=1.6)
    r_p = float(np.corrcoef(rat, dif)[0, 1])
    keep = rat > 0.2                       # without the 9% CYMAL-7 point
    r_k = float(np.corrcoef(rat[keep], dif[keep])[0, 1])
    ax2.set_xlabel("CAVER tunnel as % of our route", labelpad=10,
                   fontsize=14)
    ax2.set_ylabel("|difference| (\u00c5)", labelpad=8, fontsize=14)
    ax2.margins(0.18)
    ax2.set_axisbelow(True)
    ax2.annotate(f"r = {r_p:+.2f}, n = {len(rows)}", (0.96, 0.94),
                 xycoords="axes fraction", ha="right", va="top",
                 fontsize=13.5, color=INK2)

    fig.text(0.055, 0.135,
             "Both tools were given the same protomer with ligands stripped "
             "and the same seed, and both numbers are route bottlenecks: "
             "the narrowest point\nthe tunnel passes through, with the "
             "terminal 3 \u00c5 trimmed so the seed cavity beside the "
             "ligand is not counted as a constriction. They agree to\n"
             "0.01 and 0.03 \u00c5 on ampicillin and DDM \u00d73 and "
             "diverge by up to 1.00 \u00c5 elsewhere \u2014 but CAVER's "
             "best-ranked cluster is always far shorter than our route\n"
             "(14\u201334 \u00c5 against 52\u2013154 \u00c5), so the "
             "two are not always measuring the same passage. The right "
             f"panel is the pattern that suggests: the more of\nour route "
             f"CAVER's tunnel spans, the closer the bottlenecks "
             f"(r = {r_p:+.2f}, n = 7; {r_k:+.2f} without the CYMAL-7 point "
             "at 9% coverage, which carries most of it). With\nseven "
             "points that is suggestive, not established. CAVER ranks "
             "clustered candidate tunnels by bottleneck; ours takes the "
             "single widest path from the seed to\nbulk solvent. "
             "Per-tunnel profiles are in caver_tunnel_profiles.csv for a "
             "route-by-route comparison.",
             fontsize=13, color=INK2, va="top", linespacing=1.5)
    save(fig, "P14_caver_crosscheck")


def main():
    print("=== poster panels ===")
    panel_pockets()
    panel_occlusion()
    panel_exit_route()
    panel_ligand_size()
    panel_protomer_states()
    panel_path_occupancy()
    panel_rotamers()
    panel_conservation()
    panel_pocket_physchem()
    panel_mechanism()
    panel_mexb_rows()
    panel_regional_rmsd()
    panel_tm_overlay()
    panel_caver()
    print(f"\n  A0 portrait: each panel is ~250 mm wide as rendered; "
          f"SVG scales losslessly.")


if __name__ == "__main__":
    main()
