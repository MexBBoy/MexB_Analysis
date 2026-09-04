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
    # each substrate in the colour the poster already gives it, sampled from
    # the conserved-residues legend and table header of
    # Combio_Poster_20260828.pdf, then darkened at constant hue to clear
    # 4.5:1 on white (the poster's own colours are set on dark panels and
    # run 1.4-2.5:1 here). 3W9I is not on the poster; it takes the pink the
    # poster legend uses for DDM generally.
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
    title(fig, "Every protomer, not just the ligand-bound one",
          "Free volume at the substrate site in all three protomers of each "
          "structure, measured in one common frame.")
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
        OFF = {"MexB_DDM_3_20260730|D": (16, -17), "Amp_MexB_20260826|D": (16, -5),
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
    tx = fig.add_axes([0.715, 0.245, 0.27, 0.44]); tx.axis("off")
    tx.text(0, 1.0, f"{bind.mean() / acc.mean():.1f}\u00d7",
            ha="left", va="top", fontsize=46, fontweight="bold",
            color=BINDING, transform=tx.transAxes)
    tx.text(0, 0.80,
            "more room at the site in a\nBinding protomer than in an\n"
            "Access one. The ordering is\nBinding > Extrusion > Access\n"
            "in every structure.\n\n"
            "Our six protomers all fall\ninside the published spread\n"
            "for their state (largest\ndeparture 2.0 SD).",
            ha="left", va="top", fontsize=14, color=INK2,
            transform=tx.transAxes, linespacing=1.55)

    fig.text(0.095, 0.105,
             "One sphere fixed in the frame of the reference Binding "
             "protomer; every protomer superposed on its 39 pocket-lining "
             "C\u03b1, so the sphere sits at the\nsame anatomical position "
             "throughout. In an Access or Extrusion protomer that is not "
             "the protomer's own pocket as it would be defined in "
             "isolation - it\nis how open the substrate site is at the "
             "same place. 39 protomers from 9 structures; four of them "
             "carry two trimers in the asymmetric unit. 22XK and 22XM "
             "excluded -\nall six of their protomers fail the numbering "
             "check, being engineered chimeras at ~40% identity.",
             fontsize=13.5, color=INK2, va="top", linespacing=1.5)
    save(fig, "P5_protomer_states")


# ------------------------------------------------------------------- P6
def panel_path_occupancy():
    """Where every bound ligand sits along the transport path."""
    rows = R("ligand_environment.csv")
    if not rows:
        return
    OURS = ("Amp_MexB_20260826", "MexB_DDM_3_20260730")
    NAME = {"21FP": "chloramphenicol", "Amp_MexB_20260826": "ampicillin",
            "2V50": "DDM", "3W9I": "DDM", "21FO": "CYMAL-7",
            "3W9J": "EPI", "6IIA": "LMNG",
            "MexB_DDM_3_20260730": "DDM \u00d73"}
    SITECOL = {"DBP": APOLAR, "PBP": POLAR, "both": "#7a8891",
               "neither": "#b9c3c8"}

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
    ax.set_xlabel("depth into the porter domain (\u00c5 from the "
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
             "ones. Analyses after Lawrence et al., Nat Commun 2025;16:10601.",
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

    fig = plt.figure(figsize=(10.6, 7.2))
    title(fig, "The pocket lining does not rearrange either",
          "\u03c71 of every pocket aromatic, in all 39 protomers of all "
          "nine structures.")
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
    tx.text(0, 1.0, "9 of 10", ha="left", va="top", fontsize=40,
            fontweight="bold", color=BINDING, transform=tx.transAxes)
    tx.text(0, 0.855,
            "pocket aromatics in the\nDDM \u00d73 protomer sit within\n"
            "19\u00b0 of the published\nBinding mean. Three ligands\n"
            "at once do not rotate the\nlining.\n\n"
            "The exception, F664, is\n146\u00b0 out - but the "
            "ampicillin\nprotomer does the same\nthing, so it tracks our "
            "data,\nnot the ligand count, and\nneeds a density check.",
            ha="left", va="top", fontsize=13.5, color=INK2,
            transform=tx.transAxes, linespacing=1.45)

    fig.text(0.055, 0.095,
             "\u03c71 = N-CA-CB-CG, straight from the deposited "
             "coordinates; points near -180 and +180 are the same rotamer, "
             "split by the wrap-around. Residue 626 is a\nmethionine in "
             "MexB, so the F626 of Lawrence et al. has no counterpart here. "
             "F664 is otherwise strictly state-coupled - t in every "
             "published Access\nprotomer, g+ in every Binding and Extrusion "
             "one. Pocket volume (P4, P5) and lining rotamers together: the "
             "site neither resizes nor rearranges for a\ndifferent or a "
             "larger ligand. What changes is the conformational state, and "
             "which stations along the path are occupied (P6).",
             fontsize=13.5, color=INK2, va="top", linespacing=1.5)
    save(fig, "P7_aromatic_rotamers")


def main():
    print("=== poster panels ===")
    panel_pockets()
    panel_occlusion()
    panel_exit_route()
    panel_ligand_size()
    panel_protomer_states()
    panel_path_occupancy()
    panel_rotamers()
    print(f"\n  A0 portrait: each panel is ~250 mm wide as rendered; "
          f"SVG scales losslessly.")


if __name__ == "__main__":
    main()
