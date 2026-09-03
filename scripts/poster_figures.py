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
            "2V50": "DDM", "3W9I": "DDM", "21FO": "CYMAL-7",
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
    apo = [r for r in rows if int(r["ligand_heavy_atoms"]) == 0]
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

    # apo as a reference line: it has no ligand, so it has no depth
    if apo:
        av = num(apo[0], "volume_r16_A3")
        ax.axhline(av, color="#8c9aa1", linewidth=2, linestyle=(0, (5, 4)),
                   zorder=1)
        ax.annotate("apo (4.5 \u00c5)", (25.6, av), textcoords="offset points",
                    xytext=(0, 8), ha="left", fontsize=13, color="#8c9aa1")

    # noise floor: two independent structures with the identical ligand
    same = [r for r in bound if int(r["ligand_heavy_atoms"]) == 35]
    if len(same) == 2:
        vs = [num(r, "volume_r16_A3") for r in same]
        xs_ = [num(r, "depth_from_entrance_A") for r in same]
        ax.plot(xs_, vs, color=WARN, linewidth=4, zorder=2,
                solid_capstyle="round", alpha=.85)
        ax.annotate(f"identical ligand,\n{abs(vs[0]-vs[1]):.0f} \u00c5\u00b3 apart",
                    (float(np.mean(xs_)), float(np.mean(vs))),
                    xytext=(27.0, 1715), textcoords="data",
                    ha="left", va="center", fontsize=14, color=WARN,
                    fontweight="bold", linespacing=1.3,
                    arrowprops=dict(arrowstyle="-", color=WARN, lw=1.6,
                                    alpha=.7,
                                    connectionstyle="arc3,rad=-0.18"))

    for r in bound:
        x = num(r, "depth_from_entrance_A"); y = num(r, "volume_r16_A3")
        lo, hi = num(r, "shallowest_atom_depth_A"), num(r, "deepest_atom_depth_A")
        mine = r["pdb"] in OURS
        col = BINDING if mine else TEAL
        # the bar is the span of the ligand itself along the channel
        if np.isfinite(lo) and np.isfinite(hi):
            ax.plot([lo, hi], [y, y], color=col, linewidth=2.6, alpha=.42,
                    solid_capstyle="round", zorder=3)
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
    ax.set_xlabel("depth into the pocket from the periplasmic entrance "
                  "(\u00c5)", labelpad=26)
    ax.set_ylabel("ligand-free pocket volume (\u00c5\u00b3)")
    ax.set_xlim(25, 80); ax.margins(y=.30)
    ax.set_xticks([30, 40, 50, 60, 70])
    ax.grid(axis="x", visible=False); ax.set_axisbelow(True)
    ax.annotate("\u2190 towards the entrance", (0.0, 0.0),
                xycoords="axes fraction", xytext=(2, -42),
                textcoords="offset points", ha="left", fontsize=13,
                color=INK2, fontstyle="italic")
    ax.annotate("deeper into the distal pocket \u2192", (1.0, 0.0),
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
             "widest ligand-free entry channel of the reference protomer; "
             "bars span the\nshallowest to deepest atom of each ligand, "
             "marker area tracks ligand size. Each protomer superposed on 39 "
             "pocket-lining C\u03b1 of one reference,\nso the measuring "
             "sphere sits identically in every structure. Engineered MexB "
             "chimeras excluded.",
             fontsize=14, color=INK2, va="top", linespacing=1.5)
    save(fig, "P4_ligand_size_vs_pocket")



def main():
    print("=== poster panels ===")
    panel_pockets()
    panel_occlusion()
    panel_exit_route()
    panel_ligand_size()
    print(f"\n  A0 portrait: each panel is ~250 mm wide as rendered; "
          f"SVG scales losslessly.")


if __name__ == "__main__":
    main()
