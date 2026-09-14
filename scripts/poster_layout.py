#!/usr/bin/env python3
"""A proposed block layout for the conference poster.

Not a figure for the poster - a schematic of where the panels go, drawn at A0
proportions so the block heights can be read off directly. The point is the
grid: twelve columns with one gutter, so every panel edge lands on a column
line and the vertical gutters stay straight down the poster. The current
arrangement changes column widths between rows, which is what makes it read
as crowded even where there is room.
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from poster_figures import APOLAR, INK2, OUT, POLAR, TEAL, tint

# (label, col, span, row_y, row_h, note)
BLOCKS = [
    ("Title, authors, logos", 0, 12, 0.945, 0.055, ""),
    ("Introduction", 0, 7, 0.845, 0.098, "fix the overlapping line"),
    ("Pump schematic\n+ three states", 7, 5, 0.845, 0.098, ""),
    ("Methodology  (one band, not two)", 0, 12, 0.762, 0.080,
     "compress to a single row; gains ~5% height"),
    ("Consensus structure\n2.03 A", 0, 4, 0.545, 0.212, ""),
    ("Conformations\nAccess / Binding / Extrusion", 4, 4, 0.545, 0.212, ""),
    ("Binding pocket surfaces\nproximal / distal", 8, 4, 0.545, 0.212, ""),
    ("Maps\nampicillin + DDM density", 0, 4, 0.318, 0.222, ""),
    ("Ligand comparison matrix\n\npolarity  (row 1)\ncontacts  (row 2)\n\n"
     "one header: Ampicillin | DDM | DDM | DDM", 4, 8, 0.318, 0.222,
     "merges your two 4-across strips"),
    ("Shared residue contacts\n+ depth along the porter domain", 0, 8,
     0.108, 0.205, "the untitled ??? panel, split and named"),
    ("Conclusions", 8, 4, 0.108, 0.205, ""),
    ("References", 0, 12, 0.062, 0.040, ""),
]


def main():
    W, H = 8.41, 11.89                      # A0 proportions
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.add_patch(mp.Rectangle((0, 0), 1, 1, color="#eef6f4"))

    M, G = 0.022, 0.008                     # margin and gutter
    span = (1 - 2 * M - 11 * G) / 12        # one column

    def x_of(col, ncol):
        return M + col * (span + G), ncol * span + (ncol - 1) * G

    for (lab, col, ncol, y, h, note) in BLOCKS:
        x, w = x_of(col, ncol)
        head = lab.startswith(("Title", "References"))
        col_ = TEAL if head else "white"
        ax.add_patch(mp.FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.006",
            facecolor=col_, edgecolor=tint(TEAL, 0.35), linewidth=1.6))
        ax.text(x + w / 2, y + h / 2, lab, ha="center", va="center",
                fontsize=11 if not head else 12,
                color="white" if head else INK2,
                fontweight="bold" if head else "normal", linespacing=1.5)
        if note:
            ax.text(x + w / 2, y + 0.012, note, ha="center", va="bottom",
                    fontsize=8.5, color=APOLAR, style="italic")

    for c in range(13):                     # the grid the panels snap to
        gx = M + c * (span + G) - G / 2
        ax.plot([gx, gx], [0.05, 0.955], color=POLAR, linewidth=0.4,
                alpha=.30, zorder=0)

    ax.text(M, 0.035, "12-column grid, equal gutters — every panel edge "
            "lands on a column line", fontsize=9.5, color=INK2)
    out = os.path.join(OUT, "poster_layout_proposal")
    fig.savefig(out + ".png", dpi=160)
    fig.savefig(out + ".svg")
    plt.close(fig)
    print(f"  wrote {os.path.relpath(out)}.png / .svg")


if __name__ == "__main__":
    main()
