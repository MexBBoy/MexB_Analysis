#!/usr/bin/env python3
"""Generate the figures for the MexB analysis into results/figures/.

Palette: validated categorical slots 1-3 (blue / orange / aqua) for series,
single-hue sequential blue for magnitude, neutral grey for reference values.
Every series is direct-labelled as well as coloured, so identity is never
carried by colour alone.
"""
from __future__ import annotations

import csv
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import CXDIR, DBP, FIGURES, PBP, REFERENCE_STATES, \
    SWITCH_LOOP, TABLES

S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"     # blue, orange, aqua
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"
SEQ = LinearSegmentedColormap.from_list(
    "seqblue", ["#eef4fc", "#b9d2f2", "#7db0e6", "#4a8fdb", "#2a78d6",
                "#14508f"])

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.titlecolor": INK, "axes.titleweight": "bold",
    "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 8,
    "lines.linewidth": 2.0,
})

SHORT = {"Amp_MexB_20260826": "Ampicillin", "MexB_DDM_3_20260730": "DDM"}


def rows(name):
    p = os.path.join(TABLES, name)
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh))


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIGURES, f"{name}.{ext}"), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote results/figures/{name}.png / .svg")


def finish(ax):
    ax.set_axisbelow(True)
    ax.grid(axis="x", visible=False)


# ------------------------------------------------------------------ fig 1

def fig_states():
    st = rows("states.csv")
    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    refoff = {"Access": (-16, 14, "right"), "Binding": (16, 6, "left"),
              "Extrusion": (0, -34, "center")}
    for name, (pn, pc) in REFERENCE_STATES.items():
        dx, dy, ha = refoff[name]
        ax.scatter([pn], [pc], s=200, marker="s", facecolor="none",
                   edgecolor=MUTED, linewidth=1.6, zorder=2)
        ax.annotate(f"{name}\n(published)", (pn, pc),
                    textcoords="offset points", xytext=(dx, dy),
                    ha=ha, fontsize=8, color=MUTED)
    for i, (struct, col) in enumerate(
            [("Amp_MexB_20260826", S1), ("MexB_DDM_3_20260730", S2)]):
        xs, ys = [], []
        for r in st:
            if r["structure"] != struct:
                continue
            x, y = float(r["PN1_PN2_sep"]), float(r["PC1_PC2_sep"])
            xs.append(x)
            ys.append(y)
            dx, ha = ((10, "left") if struct == "Amp_MexB_20260826"
                      else (-10, "right"))
            ax.annotate(f"{r['chain']}", (x, y),
                        textcoords="offset points", xytext=(dx, -4),
                        ha=ha, fontsize=8.5, color=col, weight="bold")
        ax.scatter(xs, ys, s=90, color=col, zorder=3,
                   edgecolor=SURFACE, linewidth=1.5,
                   label=SHORT[struct])
    ax.set_xlabel("PN1–PN2 centroid separation (Å)")
    ax.set_ylabel("PC1–PC2 centroid separation (Å)")
    ax.margins(0.14)
    ax.set_title("Conformational state of each protomer", pad=14)
    ax.legend(loc="upper right")
    ax.text(0.5, -0.15,
            "Open squares are the published reference states. Chain D of the "
            "ampicillin model sits between\nAccess and Binding — its two "
            "diagnostics disagree, so it is reported as a conflict.",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
            color=MUTED)
    ax.set_axisbelow(True)
    save(fig, "fig1_state_diagnostics")


# ------------------------------------------------------------------ fig 2

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


def fig_tunnel_profiles():
    tr = rows("tunnels.csv")

    def meta(struct, mode, ch, seed):
        for r in tr:
            if (r["structure"] == struct and r["mode"] == mode
                    and r["chain"] == ch and r["seed"] == seed
                    and r["tunnel_rank"] == "1"):
                return r
        return {}

    panels = [("Amp_MexB_20260826", "E", "ZZ72000", "Ampicillin, chain E",
               "ampicillin"),
              ("MexB_DDM_3_20260730", "E", "LMT2001", "DDM, chain E",
               "3 DDM")]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.6), sharey=True)
    for ax, (struct, ch, tag, title, ligname) in zip(axes, panels):
        notes = []
        for mode, col, lab in (("protein", S1, "ligands stripped"),
                               ("withlig", S2, f"{ligname} in place")):
            f = os.path.join(
                CXDIR, f"{struct}_{mode}_{ch}_{tag}_t1_tunnel.pdb")
            if not os.path.exists(f):
                continue
            s, b = read_trace(f)
            ax.plot(s, b, color=col, label=lab, solid_capstyle="round")
            j = int(np.argmin(b))
            ax.scatter([s[j]], [b[j]], s=45, color=col, zorder=4,
                       edgecolor=SURFACE, linewidth=1.4)
            ax.annotate(f"{b[j]:.2f} Å", (s[j], b[j]),
                        textcoords="offset points", xytext=(0, -16),
                        ha="center", fontsize=8, color=col, weight="bold")
            m = meta(struct, mode, ch, tag)
            chan = (m.get("channel_call") or "?").split(" ")[0]
            notes.append((col, f"{lab}: {chan}, {s[-1]:.0f} Å long"))
        ax.axhline(1.4, color=MUTED, linewidth=1.2, linestyle=(0, (4, 3)))
        ax.text(0.98, 1.45, "water probe 1.4 Å", transform=
                ax.get_yaxis_transform(), ha="right", va="bottom",
                fontsize=7.5, color=MUTED)
        ax.set_ylim(top=5.15)
        y0 = 0.985
        for col, txt in notes:
            ax.text(0.03, y0, txt, transform=ax.transAxes, fontsize=7.5,
                    color=col, va="top", weight="bold")
            y0 -= 0.062
        ax.set_title(title)
        ax.set_xlabel("arc length along that tunnel's own centreline (Å)")
        finish(ax)
    axes[0].set_ylabel("local tunnel radius (Å)")
    handles = [Line2D([], [], color=S1, linewidth=2.0,
                      label="ligands stripped"),
               Line2D([], [], color=S2, linewidth=2.0,
                      label="ligand in place")]
    fig.legend(handles=handles, loc="upper center", ncol=2,
               bbox_to_anchor=(0.5, 1.005), frameon=False)
    fig.suptitle("Tunnel radius profile, widest route from the substrate "
                 "site to bulk solvent", y=1.10, fontsize=11,
                 fontweight="bold", color=INK)
    fig.text(0.5, -0.10,
             "x is cumulative distance along each curve's own path, starting "
             "at the seed in the substrate site. The two curves in a panel "
             "are DIFFERENT routes —\nwith the ligand in place the widest "
             "way out changes channel — so equal x is not the same location "
             "and the curves should not be compared point by point.\n"
             "Computed on the trimer. Radii are provisional — see REPORT.md.",
             ha="center", fontsize=7.5, color=MUTED)
    save(fig, "fig2_tunnel_profiles")


# ------------------------------------------------------------------ fig 3

def fig_occlusion():
    vr = [r for r in rows("pocket_volumes.csv")
          if r["volume_ligands_stripped_A3"]]
    labels = [f"{SHORT[r['structure']]}\nchain {r['chain']}" for r in vr]
    free = [float(r["volume_ligands_stripped_A3"]) for r in vr]
    occ = [float(r["volume_with_ligands_A3"]) for r in vr]
    x = np.arange(len(vr))
    w = 0.38
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    ax.bar(x - w / 2 - 0.01, free, w, color=S1, label="ligands stripped")
    ax.bar(x + w / 2 + 0.01, occ, w, color=S2, label="with ligands")
    for i, (f, o) in enumerate(zip(free, occ)):
        ax.text(i - w / 2 - 0.01, f + 30, f"{f:.0f}", ha="center",
                fontsize=7.5, color=INK2)
        ax.text(i + w / 2 + 0.01, o + 30, f"{o:.0f}", ha="center",
                fontsize=7.5, color=INK2)
        pct = 100 * (f - o) / f if f else 0
        if pct > 1:
            ax.text(i, max(f, o) + 155, f"−{pct:.0f}%", ha="center",
                    fontsize=8.5, color=S2, weight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("substrate-site volume (Å³)")
    ax.set_title("Occlusion of the substrate site by the bound ligand")
    ax.legend(loc="upper left")
    ax.set_ylim(0, max(free) * 1.22)
    finish(ax)
    fig.text(0.5, -0.04,
             "Grid volume within 16 Å of the site, connected to it. "
             "Internally comparable across these two structures only — "
             "not an fpocket or CASTp volume.",
             ha="center", fontsize=7.5, color=MUTED)
    save(fig, "fig3_site_occlusion")


# ------------------------------------------------------------------ fig 4

def fig_contact_matrix():
    cm = rows("contact_matrix.csv")
    ligcols = [c for c in cm[0].keys()
               if c not in ("resseq", "resname", "pocket", "switch_loop",
                            "n_ligands")]
    keep = [r for r in cm if int(r["n_ligands"]) >= 1]
    keep.sort(key=lambda r: int(r["resseq"]))
    M = np.full((len(keep), len(ligcols)), np.nan)
    for i, r in enumerate(keep):
        for j, c in enumerate(ligcols):
            if r[c]:
                M[i, j] = float(r[c])
    fig, ax = plt.subplots(figsize=(6.4, max(5.0, 0.19 * len(keep))))
    im = ax.imshow(M, cmap=SEQ.reversed(), aspect="auto", vmin=2.4,
                   vmax=4.5)
    ax.set_xticks(range(len(ligcols)))
    short = [c.split(":")[1] if ":" in c else c for c in ligcols]
    struct = [SHORT.get(c.split(":")[0], "") for c in ligcols]
    ax.set_xticklabels([f"{s}\n{t}" for s, t in zip(short, struct)],
                       fontsize=8)
    ax.set_yticks(range(len(keep)))
    ylab = []
    for r in keep:
        mark = ""
        if int(r["resseq"]) in DBP:
            mark = "  ● DBP"
        elif int(r["resseq"]) in PBP:
            mark = "  ○ PBP"
        if int(r["resseq"]) in SWITCH_LOOP:
            mark += " * switch"
        ylab.append(f"{r['resname']}{r['resseq']}{mark}")
    ax.set_yticklabels(ylab, fontsize=7)
    for i in range(len(keep)):
        for j in range(len(ligcols)):
            if not np.isnan(M[i, j]):
                ax.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center",
                        fontsize=6.5,
                        color="white" if M[i, j] < 3.3 else INK)
    ax.set_xticks(np.arange(-.5, len(ligcols), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(keep), 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=1.6)
    ax.grid(which="major", visible=False)
    ax.tick_params(which="minor", length=0)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.035)
    cb.set_label("minimum heavy-atom distance (Å)", fontsize=8)
    cb.outline.set_visible(False)
    ax.set_title("Cross-ligand contact matrix\n(blank = no contact within "
                 "4.5 Å)", fontsize=10)
    save(fig, "fig4_contact_matrix")


# ------------------------------------------------------------------ fig 5

def fig_pockets():
    pc = rows("pocket_composition.csv")
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.0))
    for ax, (key, lab, ref) in zip(axes, [
            ("pct_apolar", "apolar side-chain atoms (%)", None),
            ("mean_KD", "mean Kyte–Doolittle hydropathy", 0.0)]):
        for k, (pocket, col) in enumerate([("PBP", S3), ("DBP", S1)]):
            vals = [float(r[key]) for r in pc if r["pocket"] == pocket]
            xs = np.random.default_rng(0).normal(k, 0.055, len(vals))
            ax.scatter(xs, vals, s=70, color=col, edgecolor=SURFACE,
                       linewidth=1.4, zorder=3)
            m = float(np.mean(vals))
            ax.plot([k - 0.22, k + 0.22], [m, m], color=col, linewidth=2.5,
                    solid_capstyle="round", zorder=4)
            ax.annotate(f"{m:+.2f}" if key == "mean_KD" else f"{m:.1f}%",
                        (k + 0.26, m), fontsize=8.5, color=col,
                        weight="bold", va="center")
        if ref is not None:
            ax.axhline(ref, color=MUTED, linewidth=1.0,
                       linestyle=(0, (4, 3)))
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Proximal (PBP)\n0 aromatics",
                            "Distal (DBP)\n8 aromatics"])
        ax.set_xlim(-0.5, 1.55)
        ax.set_ylabel(lab)
        finish(ax)
    fig.suptitle("Pocket hydrophobicity — all six protomers, both "
                 "structures", y=1.02, fontsize=11, fontweight="bold",
                 color=INK)
    fig.text(0.5, -0.07,
             "The distal pocket is the hydrophobic one, as in AcrB. This "
             "settles known issue 4 against the earlier surface-rendering "
             "reading.\nEach point is one protomer; the bar is the mean.",
             ha="center", fontsize=7.5, color=MUTED)
    save(fig, "fig5_pocket_hydrophobicity")


# ------------------------------------------------------------------ fig 6

def fig_switch_gate():
    sg = [r for r in rows("switch_gate.csv")
          if r["gate_widest_clearance_A"]]
    labels = [f"{SHORT[r['structure']]}\nchain {r['chain']}" for r in sg]
    wide = [float(r["gate_widest_clearance_A"]) for r in sg]
    forced = [float(r["bottleneck_forced_through_gate_A"]) for r in sg]
    x = np.arange(len(sg))
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    for i, (w, f) in enumerate(zip(wide, forced)):
        ax.plot([i, i], [f, w], color=GRID, linewidth=6,
                solid_capstyle="round", zorder=1)
    ax.scatter(x, wide, s=95, color=S1, zorder=3, edgecolor=SURFACE,
               linewidth=1.5, label="widest clearance inside the gate")
    ax.scatter(x, forced, s=95, color=S2, zorder=3, edgecolor=SURFACE,
               linewidth=1.5,
               label="bottleneck of a path forced through the gate")
    for i, (w, f) in enumerate(zip(wide, forced)):
        ax.annotate(f"{w:.2f}", (i, w), xytext=(0, 10),
                    textcoords="offset points", ha="center", fontsize=8,
                    color=S1, weight="bold")
        ax.annotate(f"{f:.2f}", (i, f), xytext=(0, -16),
                    textcoords="offset points", ha="center", fontsize=8,
                    color=S2, weight="bold")
    ax.axhline(2.01, color=MUTED, linewidth=1.2, linestyle=(0, (4, 3)))
    ax.text(len(sg) - 0.5, 2.06, "2.01 Å expected by PROTOCOL §6",
            ha="right", va="bottom", fontsize=7.5, color=MUTED)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("clearance (Å)")
    ax.set_ylim(0, max(wide) * 1.25)
    ax.set_title("The F615 switch-loop gate is a local widening, not a "
                 "through-route")
    ax.legend(loc="upper left")
    finish(ax)
    fig.text(0.5, -0.05,
             "The gate is roomy inside (blue) but any path forced through it "
             "pinches far below 2.01 Å (orange) in every protomer.\n"
             "This is why the two tunnel validation checks fail. The gate is "
             "widest in the two Binding protomers — the same ordering in "
             "both reconstructions.",
             ha="center", fontsize=7.5, color=MUTED)
    save(fig, "fig6_switch_gate")


# ------------------------------------------------------------------ fig 7

def fig_per_residue():
    f = os.path.join(TABLES,
                     "per_residue_Amp_MexB_20260826_EvsF_TM_trimmed.csv")
    g = os.path.join(TABLES,
                     "per_residue_MexB_DDM_3_20260730_EvsF_TM_trimmed.csv")
    if not (os.path.exists(f) and os.path.exists(g)):
        return
    fig, ax = plt.subplots(figsize=(11.0, 4.0))
    bands = [(40, 328, "porter"), (181, 277, "docking DN"),
             (571, 858, "porter"), (718, 813, "docking DC")]
    for lo, hi, lab in bands:
        ax.axvspan(lo, hi, color=GRID, alpha=0.55, zorder=0, linewidth=0)
    for path, col, lab in ((f, S1, "Ampicillin"), (g, S2, "DDM")):
        with open(path) as fh:
            rs = list(csv.DictReader(fh))
        x = [int(r["resseq"]) for r in rs]
        y = [float(r["sliding_rms_A"]) for r in rs]
        ax.plot(x, y, color=col, linewidth=1.5, label=lab)
    ax.axvspan(613, 622, color=S3, alpha=0.28, zorder=1, linewidth=0)
    ax.annotate("switch loop\n613–622", (617, ax.get_ylim()[1] * 0.93),
                ha="center", fontsize=7.5, color="#0f7d57")
    ax.set_xlabel("residue number (MexB / UniProt P52002 numbering)")
    ax.set_ylabel("sliding-window Cα deviation (Å)")
    ax.set_title("Binding (chain E) versus Extrusion (chain F), fitted on "
                 "the trimmed TM domain")
    ax.legend(loc="upper right")
    ax.set_xlim(1, 1035)
    finish(ax)
    fig.text(0.5, -0.06,
             "Shaded bands mark the porter and docking domains. The two "
             "reconstructions trace each other closely, so the E→F "
             "difference replicates.",
             ha="center", fontsize=7.5, color=MUTED)
    save(fig, "fig7_per_residue_EvsF")


def main():
    os.makedirs(FIGURES, exist_ok=True)
    print("=== figures ===")
    fig_states()
    fig_tunnel_profiles()
    fig_occlusion()
    fig_contact_matrix()
    fig_pockets()
    fig_switch_gate()
    fig_per_residue()


if __name__ == "__main__":
    main()
