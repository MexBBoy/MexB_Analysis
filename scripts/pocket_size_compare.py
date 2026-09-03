#!/usr/bin/env python3
"""Is the substrate-site pocket actually a different size in the two models?

The headline occlusion numbers are computed with each structure's ligands
stripped, so they should already be like-for-like. Two things could still
make them incomparable, and both are tested here:

  1. Centring. Each site sphere is placed at that chain's own DBP/PBP
     midpoint. If the porter domains differ slightly the spheres sit in
     slightly different places, which moves the volume for reasons that have
     nothing to do with pocket size. Fixed by superposing on the pocket
     lining and measuring both structures in one common frame.

  2. Method sensitivity. A grid volume depends on grid step, probe radius and
     sphere radius. Unless the Amp-vs-DDM difference is larger than the
     spread those parameters induce, "different" is not a claim worth making.

Only chain E is a valid comparison: it is Binding in both models and the two
agree to 0.86 A whole-protomer. Chain D is a different state in each model
and chain F of the ampicillin map is not supported by its density, so both
are reported for completeness but must not be read as pocket-size changes.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pockets import site_volume
from mexb_common import (DBP, PBP, SWITCH_LOOP, TABLES, Structure, apply_rt,
                         centroid, coords, density_warning, fmt, kabsch,
                         write_csv)

AMP, DDM = "Amp_MexB_20260826", "MexB_DDM_3_20260730"
LINING = sorted(set(DBP) | set(PBP) | set(SWITCH_LOOP))


def load():
    return (Structure(os.path.join("structures", f"{AMP}.pdb")),
            Structure(os.path.join("structures", f"{DDM}.pdb")))


def heavy_protein(s):
    return [a for a in s.protein_atoms if not a.is_hydrogen]


def common_frame(mob, ref, chain):
    """Superpose `mob` onto `ref` on the pocket-lining CA of `chain`.

    Returns the mobile structure's protein heavy-atom coordinates in the
    reference frame, plus the fit RMSD.
    """
    mca, rca = mob.ca(chain), ref.ca(chain)
    common = [r for r in LINING if r in mca and r in rca]
    M = np.array([mca[r] for r in common])
    T = np.array([rca[r] for r in common])
    R, t = kabsch(M, T)
    moved = apply_rt(R, t, coords(heavy_protein(mob)))
    fit = float(np.sqrt(((apply_rt(R, t, M) - T) ** 2).sum(axis=1).mean()))
    return moved, fit, len(common)


class Frozen:
    """Minimal stand-in exposing .xyz/.element for site_volume."""
    __slots__ = ("xyz", "element", "hetatm", "is_hydrogen")

    def __init__(self, xyz, element):
        self.xyz = tuple(xyz)
        self.element = element
        self.hetatm = False
        self.is_hydrogen = False


def as_atoms(xyz, template):
    return [Frozen(p, a.element) for p, a in zip(xyz, template)]


def main():
    amp, ddm = load()
    print("=== substrate-site pocket size: ampicillin vs DDM ===")
    print("    both measured with every ligand stripped\n")

    rows = []
    for chain in ("D", "E", "F"):
        note = []
        if density_warning(AMP, chain):
            note.append("Amp chain not supported by density")
        if chain == "D":
            note.append("different state in each model (Amp D ambiguous, "
                        "DDM D access)")
        valid = chain == "E"

        amp_at = heavy_protein(amp)
        moved, fit, nfit = common_frame(ddm, amp, chain)
        ddm_at = as_atoms(moved, heavy_protein(ddm))

        # one centre for both, taken from the reference structure
        ca = amp.ca(chain)
        centre = 0.5 * (centroid(ca, DBP) + centroid(ca, PBP))

        print(f"  chain {chain}"
              + (f"   [{'; '.join(note)}]" if note else "   [valid comparison]"))
        print(f"    common frame: fit on {nfit} pocket-lining CA, "
              f"RMSD {fit:.2f} A")

        diffs = []
        for step in (0.4, 0.5, 0.6):
            for probe in (1.2, 1.4, 1.6):
                for radius in (14.0, 16.0, 18.0):
                    va, _ = site_volume(amp_at, centre, radius=radius,
                                        step=step, probe=probe)
                    vd, _ = site_volume(ddm_at, centre, radius=radius,
                                        step=step, probe=probe)
                    if va is None or vd is None:
                        continue
                    diffs.append((va, vd, vd - va, step, probe, radius))
                    rows.append([chain, "yes" if valid else "no", fmt(step, 2),
                                 fmt(probe, 2), fmt(radius, 1), fmt(va, 0),
                                 fmt(vd, 0), fmt(vd - va, 0),
                                 fmt(100 * (vd - va) / va, 1)])
        d = np.array([x[2] for x in diffs])
        a = np.array([x[0] for x in diffs])
        rel = 100 * d / a
        print(f"    across {len(diffs)} parameter sets "
              f"(step 0.4-0.6, probe 1.2-1.6, radius 14-18 A):")
        print(f"      Amp {a.min():.0f}-{a.max():.0f} A^3, "
              f"DDM {np.array([x[1] for x in diffs]).min():.0f}-"
              f"{np.array([x[1] for x in diffs]).max():.0f} A^3")
        print(f"      DDM − Amp: median {np.median(d):+.0f} A^3 "
              f"({np.median(rel):+.1f}%), range {d.min():+.0f} to "
              f"{d.max():+.0f}")
        sign = "consistent" if (d > 0).all() or (d < 0).all() else \
            "SIGN FLIPS across parameters"
        print(f"      direction: {sign}\n")

    # --- occlusion robustness -------------------------------------------
    # Occlusion is a within-structure comparison, so the sphere ought to
    # cancel. It does for DDM and does not for ampicillin, and the reason
    # matters: three DDM molecules saturate the pocket, so free volume
    # collapses to a constant and the PERCENTAGE is stable. Ampicillin
    # displaces a fixed amount of a much larger pocket, so its ABSOLUTE
    # displaced volume is stable and its percentage is not.
    print("  --- occlusion robustness, chain E ---")
    orows = []
    for s_, tag in ((amp, "Amp"), (ddm, "DDM")):
        prot = heavy_protein(s_)
        lig = [a for a in s_.het_atoms if not a.is_hydrogen]
        ca = s_.ca("E")
        cen = 0.5 * (centroid(ca, DBP) + centroid(ca, PBP))
        disp, pct = [], []
        for radius in (14.0, 16.0, 18.0, 20.0):
            vf, _ = site_volume(prot, cen, radius=radius, step=0.5,
                                probe=1.4)
            vo, _ = site_volume(prot + lig, cen, radius=radius, step=0.5,
                                probe=1.4)
            disp.append(vf - vo)
            pct.append(100 * (vf - vo) / vf)
            orows.append([tag, fmt(radius, 1), fmt(vf, 0), fmt(vo, 0),
                          fmt(vf - vo, 0), fmt(100 * (vf - vo) / vf, 1)])
        print(f"    {tag}: displaced volume {min(disp):.0f}-{max(disp):.0f} "
              f"A^3, occluded {min(pct):.1f}-{max(pct):.1f}%")
    print("    -> quote ampicillin as an absolute volume, DDM as a "
          "percentage\n")
    write_csv(os.path.join(TABLES, "occlusion_robustness.csv"),
              ["structure", "sphere_radius_A", "free_A3", "with_ligand_A3",
               "displaced_A3", "occluded_pct"], orows)

    write_csv(os.path.join(TABLES, "pocket_size_compare.csv"),
              ["chain", "valid_comparison", "grid_step_A", "probe_A",
               "sphere_radius_A", "amp_volume_A3", "ddm_volume_A3",
               "difference_A3", "difference_pct"], rows)
    print("wrote results/tables/pocket_size_compare.csv")


if __name__ == "__main__":
    main()
