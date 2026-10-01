#!/usr/bin/env python3
"""Figure 1 for the MtrE binder-design write-up.

Panel a is measured, not drawn: the channel wall comes from the coordinates of
targets/mtre_trimer_framed.pdb, so the depths and the open lumen are real.
Panels b-d are mock-ups of the three binder concepts on that same scale, to
show what each campaign is trying to achieve mechanistically.

Usage:
    python3 scripts/figure_targets.py --out figures/fig1_targets.png
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent.parent

WALL = "#5a5f66"
WALL_FILL = "#e3e5e8"
MEMBRANE = "#efe7d8"
MEMBRANE_EDGE = "#d6c9ae"
BINDER = "#1f6fb4"
BINDER_FILL = "#bcd8ef"
BLOCK = "#c0392b"
FLOW = "#7a8089"
INK = "#1a1d21"

MEMBRANE_Z = (35.0, 60.0)
LUMEN_FLOOR = 34.0


def wall_profile(step: float = 3.0) -> list[dict]:
    """Outer and lumen radius of the trimer as a function of depth."""
    model = PDBParser(QUIET=True).get_structure("m", ROOT / "targets" / "mtre_trimer_framed.pdb")[0]
    atoms = [a for ch in model for r in ch for a in r if a.element != "H"]
    xyz = np.array([a.coord for a in atoms])
    vdw = np.array([{"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}.get(a.element, 1.7)
                    for a in atoms])
    rows = []
    for z in np.arange(-62, 69, step):
        sel = np.abs(xyz[:, 2] - z) < step
        if sel.sum() < 8:
            continue
        rad = np.hypot(xyz[sel, 0], xyz[sel, 1])
        rows.append({"z": float(z),
                     "pore": float((rad - vdw[sel]).min()),
                     "outer": float(np.percentile(rad, 99))})
    return rows


def draw_channel(ax, rows, zlim, lw=1.1, label_sites=False):
    """The protein cross-section: two mirrored bands between lumen and outer wall."""
    rows = [r for r in rows if zlim[0] - 4 <= r["z"] <= zlim[1] + 4]
    for sign in (-1, 1):
        outer = [(sign * r["outer"], r["z"]) for r in rows]
        inner = [(sign * r["pore"], r["z"]) for r in reversed(rows)]
        ax.add_patch(Polygon(outer + inner, closed=True, facecolor=WALL_FILL,
                             edgecolor=WALL, linewidth=lw, zorder=2))

    ax.axhspan(MEMBRANE_Z[0], MEMBRANE_Z[1], color=MEMBRANE, zorder=0)
    ax.axhline(MEMBRANE_Z[0], color=MEMBRANE_EDGE, lw=0.8, zorder=1)
    ax.axhline(MEMBRANE_Z[1], color=MEMBRANE_EDGE, lw=0.8, zorder=1)


def efflux_arrow(ax, z0, z1, blocked_at=None):
    """Drug efflux runs bottom to top. A blocked arrow stops short of the
    binder, with the cross on the arrow rather than on the binder."""
    top = z1 if blocked_at is None else blocked_at - 7
    ax.add_patch(FancyArrowPatch((0, z0), (0, top), arrowstyle="-|>", mutation_scale=12,
                                 color=FLOW, lw=2.0, zorder=4, alpha=0.9))
    if blocked_at is not None:
        for sign in (-1, 1):
            ax.plot([-6, 6], [top - 4 * sign, top + 4 * sign], color=BLOCK, lw=2.6,
                    zorder=9, solid_capstyle="round")


def panel_a(ax, rows, sites):
    draw_channel(ax, rows, (-62, 70))
    ax.set_xlim(-50, 50)
    ax.set_ylim(-66, 86)

    for z, txt in ((66, None), (LUMEN_FLOOR, None)):
        ax.plot([-16, 16], [z, z], color=BINDER, lw=1.4, zorder=5)

    ax.add_patch(Polygon([(-11, LUMEN_FLOOR), (11, LUMEN_FLOOR), (11, 66), (-11, 66)],
                         closed=True, facecolor=BINDER_FILL, alpha=0.55,
                         edgecolor=BINDER, linestyle=(0, (4, 3)), lw=1.2, zorder=3))

    asp = -46.0
    ax.plot([-11, 11], [asp, asp], color=BLOCK, lw=3.2, zorder=5,
            solid_capstyle="round")

    ax.annotate("reachable lumen:\n32 Å, no constriction",
                xy=(11, 52), xytext=(15, 74), fontsize=7.2, color=INK, ha="left",
                arrowprops=dict(arrowstyle="-", color=WALL, lw=0.8))
    ax.annotate("aspartate ring\nD422/D425 — the only\ngate, and periplasmic",
                xy=(-11, asp), xytext=(-49, -60), fontsize=7.2, color=INK, ha="left",
                arrowprops=dict(arrowstyle="-", color=WALL, lw=0.8))

    ax.text(-48, 78, "EXTRACELLULAR", fontsize=7, color=FLOW, weight="bold")
    ax.text(-48, 45, "OUTER\nMEMBRANE", fontsize=7, color=FLOW, weight="bold", va="center")
    ax.text(-48, 8, "PERIPLASM", fontsize=7, color=FLOW, weight="bold")

    ax.set_ylabel("depth along channel, Å", fontsize=8)
    ax.set_yticks([-60, -40, -20, 0, 20, 40, 60])
    ax.tick_params(axis="y", labelsize=7.5)
    ax.set_xticks([])
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(WALL)
    ax.set_title("a   MtrE, measured from 4MT0", fontsize=8.6, loc="left",
                 weight="bold", color=INK, pad=6)


def helix(ax, z0, z1, width=7.0, turns=7):
    t = np.linspace(0, 1, 260)
    z = z0 + (z1 - z0) * t
    x = width * np.sin(2 * np.pi * turns * t)
    ax.plot(x, z, color=BINDER, lw=2.0, zorder=7, solid_capstyle="round")


def panel_plug(ax, rows):
    draw_channel(ax, rows, (20, 76))
    efflux_arrow(ax, 21, 70, blocked_at=42)
    helix(ax, 40, 68)
    ax.add_patch(FancyBboxPatch((-13, 71), 26, 8, boxstyle="round,pad=1.2,rounding_size=3",
                                facecolor=BINDER_FILL, edgecolor=BINDER, lw=1.3, zorder=8))
    ax.text(0, 75, "cap", fontsize=7, ha="center", va="center", color=INK, zorder=9)
    ax.annotate("single helix, 60–110 aa,\nfills the lumen",
                xy=(7, 55), xytext=(17, 68), fontsize=7.2, color=INK, ha="left",
                arrowprops=dict(arrowstyle="-", color=WALL, lw=0.8))
    ax.set_title("b   A — lumen plug: a helix fills the conduit",
                 fontsize=8.2, loc="left", weight="bold", color=INK, pad=4)


def panel_macrocycle(ax, rows):
    draw_channel(ax, rows, (20, 76))
    efflux_arrow(ax, 21, 70, blocked_at=60)
    for rx, ry, lw in ((9.0, 4.0, 2.4), (5.2, 2.1, 1.4)):
        ax.add_patch(Ellipse((0, 65), 2 * rx, 2 * ry, facecolor="white",
                             edgecolor=BINDER, lw=lw, zorder=8))
    ax.annotate("8–16mer macrocycle\ncorks the 15 Å mouth",
                xy=(9, 65), xytext=(17, 70), fontsize=7.2, color=INK, ha="left",
                arrowprops=dict(arrowstyle="-", color=WALL, lw=0.8))
    ax.set_title("c   B — macrocycle: corks the mouth",
                 fontsize=8.2, loc="left", weight="bold", color=INK, pad=4)


def panel_crown(ax, rows):
    draw_channel(ax, rows, (20, 76))
    efflux_arrow(ax, 21, 80)
    for dx in (11, 17, 23):
        ax.add_patch(FancyBboxPatch((dx, 69), 4, 10,
                                    boxstyle="round,pad=0.8,rounding_size=2",
                                    facecolor=BINDER_FILL, edgecolor=BINDER, lw=1.2, zorder=8))
    ax.annotate("minibinder clamps the\nLoop1/Loop2 seam\n(F326, W333)",
                xy=(14, 74), xytext=(-48, 68), fontsize=7.2, color=INK, ha="left",
                arrowprops=dict(arrowstyle="-", color=WALL, lw=0.8))
    ax.text(4, 44, "efflux continues", fontsize=7.2, color=BLOCK, ha="left", zorder=9)
    ax.set_title("d   C — crown groove: binds, but may not block",
                 fontsize=8.2, loc="left", weight="bold", color=INK, pad=4)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="figures/fig1_targets.png")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    rows = wall_profile()
    sites = json.loads((ROOT / "targets" / "mtre_sites.json").read_text())

    fig = plt.figure(figsize=(7.4, 7.6))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.25], wspace=0.18, hspace=0.38,
                          left=0.085, right=0.985, top=0.945, bottom=0.03)

    panel_a(fig.add_subplot(gs[:, 0]), rows, sites)

    for i, fn in enumerate((panel_plug, panel_macrocycle, panel_crown)):
        ax = fig.add_subplot(gs[i, 1])
        fn(ax, rows)
        ax.set_xlim(-50, 50)
        ax.set_ylim(20, 84)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ("top", "right", "bottom", "left"):
            ax.spines[s].set_visible(False)

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=args.dpi, facecolor="white")
    print(f"wrote {out} ({args.dpi} dpi)")


if __name__ == "__main__":
    main()
