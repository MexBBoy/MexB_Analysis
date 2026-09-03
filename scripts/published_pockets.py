#!/usr/bin/env python3
"""Does the MexB binding pocket enlarge to accommodate bigger ligands?

Tests it across every published substrate- or detergent-bound MexB structure
plus the two in this project, spanning a 5x range of bound ligand size
(chloramphenicol, 20 heavy atoms, to three DDM molecules, 105).

Method. For each structure the ligand-bound protomer is identified, its
numbering checked against the expected pocket residue identities, and the
protomer superposed onto a single reference on the pocket-lining CA. The
ligand-free site volume is then measured in that one common frame, so the
measuring sphere sits in the same place for every structure. Volume is
reported at several sphere radii because a single radius is not meaningful
on its own (PROTOCOL known issue 10).

Correlating that volume against ligand size answers the question directly.
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pockets import site_volume
from mexb_common import (CXDIR, DBP, PBP, SWITCH_LOOP, STRUCT_DIR, TABLES,
                         Structure, apply_rt, centroid, coords, fmt, kabsch,
                         write_csv)

PDBDIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "work", "pdb")
LINING = sorted(set(DBP) | set(PBP) | set(SWITCH_LOOP))
CRYO = {"HOH", "GOL", "EDO", "SO4", "PO4", "CL", "NA", "MG", "K", "ACT",
        "MPD", "PEG", "DMS", "IOD", "CA", "ZN", "TRS", "FMT", "NO3", "BME",
        "EPE", "UNX", "PG4", "1PE", "P6G", "PGE"}
# a numbering sanity check: these must be these residues in MexB numbering
EXPECT = {151: "LYS", 178: "PHE", 610: "PHE", 615: "PHE", 617: "PHE",
          620: "ARG", 630: "MET"}
LABEL = {
    "21FP": "chloramphenicol (2.89 A)", "21FO": "CYMAL-7 (2.30 A)",
    "2V50": "DDM x1 (3.00 A)", "3W9I": "DDM x1 (2.71 A)",
    "3W9J": "pyridopyrimidine EPI (3.15 A)", "6IIA": "LMNG (2.91 A)",
    "22XK": "apo (3.60 A)", "22XM": "apo (3.55 A)", "6T7S": "apo (4.50 A)",
    "Amp_MexB_20260826": "ampicillin (2.19 A, this work)",
    "MexB_DDM_3_20260730": "DDM x3 (2.11 A, this work)",
}
RADII = (14.0, 16.0, 18.0, 20.0)
# widest ligand-free route out of the reference protomer: it leaves by the
# PC1/PC2 periplasmic cleft, i.e. the substrate entry channel (CH1). Its
# coordinates are already in the reference frame, so the common-frame
# transform puts every ligand onto the same axis.
REF_CHANNEL = os.path.join(CXDIR,
                           "Amp_MexB_20260826_protein_E_ZZ72000_t1_tunnel.pdb")


def load_channel(path=REF_CHANNEL):
    """Entry-channel centreline as (points, arc-length-from-site, total)."""
    if not os.path.exists(path):
        return None
    pts = []
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")):
            pts.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    if len(pts) < 10:
        return None
    P = np.asarray(pts, float)
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(seg)])
    return P, arc, float(arc[-1])


def channel_depth(cen, chan):
    """(depth from the periplasmic entrance, offset from the centreline).

    Arc length runs 0 at the deep site seed to `total` at the bulk exit, so
    depth into the pocket measured from the entrance is total - arc.
    """
    if chan is None:
        return None, None
    P, arc, total = chan
    d = np.linalg.norm(P - np.asarray(cen, float), axis=1)
    i = int(np.argmin(d))
    return float(total - arc[i]), float(d[i])


def numbering_ok(s, ch):
    good = bad = 0
    for r, want in EXPECT.items():
        at = s.residue_atoms(ch, r)
        if not at:
            continue
        (good := good + 1) if at[0].resname == want else (bad := bad + 1)
    return good, bad


def pocket_ligands(s):
    """[(chain, resname, n_heavy, atoms)] for ligands in a porter pocket."""
    out = []
    for (ch, rs, rn, ats) in s.ligands():
        heavy = [a for a in ats if not a.is_hydrogen]
        if rn in CRYO or len(heavy) < 6:
            continue
        cen = coords(heavy).mean(0)
        best = None
        for pch in s.chains:
            ca = s.ca(pch)
            pk = [ca[r] for r in DBP + PBP if r in ca]
            if len(pk) < 10:
                continue
            d = float(np.linalg.norm(cen - np.mean(pk, axis=0)))
            if best is None or d < best[1]:
                best = (pch, d)
        if best and best[1] < 18:
            out.append((best[0], rn, heavy))
    return out


def binding_protomer(s):
    """Chain whose PC1-PC2 separation is nearest the binding reference."""
    best, bd = None, 1e9
    for ch in s.chains:
        ca = s.ca(ch)
        c1, c2 = centroid(ca, [r for r in range(571, 667)]), None
        pc2 = [r for r in list(range(679, 718)) + list(range(814, 819))
               + list(range(821, 859))]
        c2 = centroid(ca, pc2)
        if c1 is None or c2 is None:
            continue
        d = abs(float(np.linalg.norm(c1 - c2)) - 29.1)
        if d < bd:
            best, bd = ch, d
    return best


class Frozen:
    __slots__ = ("xyz", "element", "hetatm", "is_hydrogen")

    def __init__(self, xyz, element):
        self.xyz = tuple(xyz); self.element = element
        self.hetatm = False; self.is_hydrogen = False


ENTRIES = ["2V50", "3W9I", "3W9J", "6IIA", "21FO", "21FP", "22XK", "22XM",
           "6T7S"]


def ensure_pdbs():
    """Fetch the published MexB entries if they are not already here."""
    import urllib.request
    os.makedirs(PDBDIR, exist_ok=True)
    for e in ENTRIES:
        dest = os.path.join(PDBDIR, f"{e}.pdb")
        if os.path.exists(dest) and os.path.getsize(dest) > 10000:
            continue
        url = f"https://files.rcsb.org/download/{e}.pdb"
        try:
            urllib.request.urlretrieve(url, dest)
            print(f"  fetched {e}")
        except Exception as exc:
            print(f"  could not fetch {e}: {exc}")


def main():
    ensure_pdbs()
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    ref_ch = "E"
    rca = ref.ca(ref_ch)
    ref_centre = 0.5 * (centroid(rca, DBP) + centroid(rca, PBP))
    chan = load_channel()
    if chan is None:
        print("  (entry channel trace missing - run tunnels.py first; "
              "depth columns will be blank)")
    # cross-check axis: proximal-pocket centroid -> distal-pocket centroid
    pbp_c, dbp_c = centroid(rca, PBP), centroid(rca, DBP)
    ax = dbp_c - pbp_c
    ax = ax / np.linalg.norm(ax)

    targets = [(os.path.join(STRUCT_DIR, f"{n}.pdb"), n)
               for n in ("Amp_MexB_20260826", "MexB_DDM_3_20260730")]
    targets += [(f, os.path.basename(f)[:-4])
                for f in sorted(glob.glob(os.path.join(PDBDIR, "*.pdb")))]

    print("=== does the pocket enlarge for bigger ligands? ===")
    print(f"    common frame: every protomer superposed on the "
          f"{len(LINING)} pocket-lining CA of {ref.name} chain {ref_ch}\n")
    rows = []
    for path, pid in targets:
        s = Structure(path)
        ligs = pocket_ligands(s)
        if ligs:
            ch = ligs[0][0]
            same = [l for l in ligs if l[0] == ch]
            n_heavy = sum(len(l[2]) for l in same)
            names = "+".join(sorted({l[1] for l in same}))
            n_lig = len(same)
            lig_atoms = [a for l in same for a in l[2]]
        else:
            ch = binding_protomer(s)
            n_heavy, names, n_lig, lig_atoms = 0, "apo", 0, []
        if ch is None:
            print(f"  {pid}: no usable protomer, skipped"); continue
        good, bad = numbering_ok(s, ch)
        if bad or good < 4:
            print(f"  {pid}: numbering check FAILED "
                  f"({good} match, {bad} mismatch) - skipped")
            continue

        mca = s.ca(ch)
        common = [r for r in LINING if r in mca and r in rca]
        if len(common) < 25:
            print(f"  {pid}: only {len(common)} lining CA in common - skipped")
            continue
        M = np.array([mca[r] for r in common])
        T = np.array([rca[r] for r in common])
        R, t = kabsch(M, T)
        fit = float(np.sqrt(((apply_rt(R, t, M) - T) ** 2).sum(1).mean()))
        prot = [a for a in s.protein_atoms if not a.is_hydrogen
                and a.chain == ch]
        moved = apply_rt(R, t, coords(prot))
        atoms = [Frozen(p, a.element) for p, a in zip(moved, prot)]

        # put the ligand in the same common frame and ask how deep it sits
        if lig_atoms:
            lxyz = apply_rt(R, t, coords(lig_atoms))
            lcen = lxyz.mean(0)
            depth, offset = channel_depth(lcen, chan)
            axial = float(np.dot(lcen - pbp_c, ax))
            # the ligands are elongated, so record the span they occupy
            per = [channel_depth(x, chan)[0] for x in lxyz] \
                if chan is not None else []
            deep = max(per) if per else None
            shal = min(per) if per else None
        else:
            depth = offset = axial = deep = shal = None

        vols = {}
        for rad in RADII:
            v, _ = site_volume(atoms, ref_centre, radius=rad, step=0.5,
                               probe=1.4)
            vols[rad] = v
        dtxt = "  apo " if depth is None else f"{depth:5.1f}"
        print(f"  {pid:22} {LABEL.get(pid,''):32} ch{ch} "
              f"lig {n_heavy:>3} atoms  depth {dtxt} A  fit {fit:4.2f} A   "
              + "  ".join(f"r{int(r)}={vols[r]:.0f}" for r in RADII))
        rows.append([pid, LABEL.get(pid, ""), ch, names, n_lig, n_heavy,
                     fmt(depth), fmt(shal), fmt(deep), fmt(offset), fmt(axial),
                     fmt(fit), len(common)]
                    + [fmt(vols[r], 0) for r in RADII])

    write_csv(os.path.join(TABLES, "published_pockets.csv"),
              ["pdb", "description", "chain", "ligands", "n_ligands",
               "ligand_heavy_atoms", "depth_from_entrance_A",
               "shallowest_atom_depth_A", "deepest_atom_depth_A",
               "offset_from_channel_A",
               "axial_PBP_to_DBP_A", "fit_rmsd_A", "n_lining_CA"]
              + [f"volume_r{int(r)}_A3" for r in RADII], rows)

    print("\n  --- correlation of pocket volume with ligand size ---")
    # apo entries carry no ligand and only one exists (6T7S, 4.5 A); leaving
    # it in would manufacture a correlation out of a single low-resolution
    # point, so the size trend is fitted on the ligand-bound structures only.
    bound = [r for r in rows if float(r[5]) > 0]
    lig = np.array([float(r[5]) for r in bound])
    for i, rad in enumerate(RADII):
        v = np.array([float(r[13 + i]) for r in bound])
        ok = np.isfinite(v) & np.isfinite(lig)
        if ok.sum() < 4:
            continue
        r_p = float(np.corrcoef(lig[ok], v[ok])[0, 1])
        sl = float(np.polyfit(lig[ok], v[ok], 1)[0])
        print(f"    sphere {int(rad)} A: Pearson r = {r_p:+.3f}, "
              f"slope {sl:+.2f} A^3 per ligand heavy atom "
              f"(volume spread {v[ok].min():.0f}-{v[ok].max():.0f})")
    print("\n  --- how far into the pocket each ligand sits ---")
    dep = np.array([float(r[6]) if r[6] not in ("", "NA") else np.nan
                    for r in bound])
    for i, rad in enumerate(RADII):
        v = np.array([float(r[13 + i]) for r in bound])
        ok = np.isfinite(v) & np.isfinite(dep)
        if ok.sum() < 4:
            continue
        r_p = float(np.corrcoef(dep[ok], v[ok])[0, 1])
        print(f"    sphere {int(rad)} A: volume vs depth Pearson "
              f"r = {r_p:+.3f}  (depth range "
              f"{dep[ok].min():.1f}-{dep[ok].max():.1f} A)")
    print("\nwrote results/tables/published_pockets.csv")


if __name__ == "__main__":
    main()
