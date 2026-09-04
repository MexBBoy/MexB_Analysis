#!/usr/bin/env python3
"""Per-ligand analysis across every MexB structure, after Lawrence et al.
(Nat Commun 2025;16:10601), the MdtF paper this project follows.

Two of that paper's structure-only analyses, applied to every bound ligand in
every MexB entry rather than to one ligand in one structure:

1. Where along the transport path each ligand sits (their CH1/CH2/CH3 and
   PBP/DBP assignment, which they reach by docking; here it is measured
   directly from the coordinates as depth along the reference entry channel
   plus which lining set the ligand actually contacts).

2. What kind of contacts hold it there (their Fig. 3E treatment of R6G:
   "exclusively hydrophobic ... stabilised by cation-pi (F626) and pi-pi
   (F178)"), typed per ligand: aromatic stacking, hydrophobic, hydrogen
   bonding, and which aromatic residues are engaged.

The point of doing this per ligand rather than per structure is the DDM x3
model: three ligands in one porter domain can only be compared with the
one-ligand structures if each is scored separately.

Writes results/tables/ligand_environment.csv.
"""
from __future__ import annotations

import glob
import os
import sys
from itertools import combinations

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from published_pockets import (LABEL, LINING, PDBDIR, channel_depth,
                               ensure_pdbs, load_channel, numbering_ok,
                               pocket_ligands)
from mexb_common import (DBP, PBP, SWITCH_LOOP, STRUCT_DIR, TABLES, Structure,
                         apply_rt, centroid, coords, fmt, kabsch, write_csv)

AROMATIC = {"PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
            "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
            "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
            "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"]}
CONTACT = 4.5      # heavy-atom contact cutoff, A
HBOND = 3.5        # N/O to N/O, A
STACK_CENTRE = 5.5  # ring-centroid separation for pi-pi, A
CATPI = 6.0        # cation to ring centroid, A


def ring_frame(atoms, resname):
    """(centroid, unit normal) of a side-chain aromatic ring, or None."""
    names = AROMATIC.get(resname)
    if not names:
        return None
    pts = np.array([a.xyz for a in atoms if a.name.strip() in names])
    if len(pts) < 5:
        return None
    c = pts.mean(0)
    # normal = smallest-variance direction of the ring atoms
    _, _, vt = np.linalg.svd(pts - c)
    return c, vt[2] / np.linalg.norm(vt[2])


def ligand_rings(heavy):
    """Planar 5- and 6-membered carbon/hetero rings in a ligand.

    Bonds are inferred by distance, which is all a PDB HETATM record
    supports. Only rings that are actually flat are kept, so an aliphatic
    sugar or cyclohexyl does not get counted as aromatic.
    """
    X = np.array([a.xyz for a in heavy])
    n = len(X)
    if n < 5:
        return []
    d = np.linalg.norm(X[:, None] - X[None, :], axis=-1)
    adj = [set(np.where((d[i] > 0.1) & (d[i] < 1.75))[0]) for i in range(n)]
    rings, seen = [], set()
    for i in range(n):                      # short cycles through each atom
        for a, b in combinations(sorted(adj[i]), 2):
            for size in (5, 6):
                paths = [[a]]
                for _ in range(size - 3):
                    nxt = []
                    for p in paths:
                        for k in adj[p[-1]]:
                            if k != i and k not in p:
                                nxt.append(p + [k])
                    paths = nxt
                for p in paths:
                    if b in adj[p[-1]]:
                        cyc = tuple(sorted([i, b] + p))
                        if len(cyc) == size and cyc not in seen:
                            seen.add(cyc)
                            rings.append(list(cyc))
    out = []
    for cyc in rings:
        pts = X[cyc]
        c = pts.mean(0)
        _, s, vt = np.linalg.svd(pts - c)
        if s[2] / max(s[0], 1e-9) < 0.08:   # flat -> aromatic or conjugated
            out.append((c, vt[2] / np.linalg.norm(vt[2])))
    return out


def main():
    ensure_pdbs()
    ref = Structure(os.path.join(STRUCT_DIR, "Amp_MexB_20260826.pdb"))
    rca = ref.ca("E")
    chan = load_channel()
    if chan is None:
        print("  entry channel trace missing - run tunnels.py first")
        return

    targets = [(os.path.join(STRUCT_DIR, f"{n}.pdb"), n)
               for n in ("Amp_MexB_20260826", "MexB_DDM_3_20260730")]
    targets += [(f, os.path.basename(f)[:-4])
                for f in sorted(glob.glob(os.path.join(PDBDIR, "*.pdb")))]

    print("=== every bound ligand, one at a time ===")
    print("    depth measured along the same reference entry channel; "
          "contacts at 4.5 A\n")
    rows = []
    for path, pid in targets:
        s = Structure(path)
        ligs = pocket_ligands(s)
        if not ligs:
            continue
        by_chain = {}
        for (ch, rn, heavy) in ligs:
            by_chain.setdefault(ch, []).append((rn, heavy))
        for ch, items in sorted(by_chain.items()):
            good, bad = numbering_ok(s, ch)
            if bad or good < 4:
                print(f"  {pid} {ch}: numbering check failed - skipped")
                continue
            mca = s.ca(ch)
            common = [r for r in LINING if r in mca and r in rca]
            if len(common) < 25:
                continue
            R, t = kabsch(np.array([mca[r] for r in common]),
                          np.array([rca[r] for r in common]))

            prot = [a for a in s.protein_atoms
                    if not a.is_hydrogen and a.chain == ch]
            P = coords(prot)
            tree = cKDTree(P)
            # side-chain aromatic rings of this protomer, in place
            rings_prot = []
            for rid in sorted({a.resseq for a in prot}):
                at = s.residue_atoms(ch, rid)
                if not at:
                    continue
                fr = ring_frame(at, at[0].resname)
                if fr:
                    rings_prot.append((rid, at[0].resname, fr[0], fr[1]))

            # order the ligands of this protomer by depth, deepest last
            scored = []
            for rn, heavy in items:
                L = coords(heavy)
                cen = L.mean(0)
                depth, offset = channel_depth(apply_rt(R, t, cen[None])[0],
                                              chan)
                scored.append((depth, rn, heavy, L, cen, offset))
            scored.sort(key=lambda x: x[0] if x[0] is not None else 0)

            for k, (depth, rn, heavy, L, cen, offset) in enumerate(scored):
                idx = tree.query_ball_point(L, CONTACT)
                touched = {prot[j].resseq for g in idx for j in g}
                n_dbp = len(touched & set(DBP))
                n_pbp = len(touched & set(PBP))
                n_sw = len(touched & set(SWITCH_LOOP))
                site = ("DBP" if n_dbp > n_pbp else
                        "PBP" if n_pbp > n_dbp else
                        "both" if n_dbp else "neither")

                # contact typing, after their Fig. 3E
                lig_el = [a.element.strip().upper() for a in heavy]
                pol_l = np.array([e in ("N", "O") for e in lig_el])
                apol_l = np.array([e in ("C", "S") for e in lig_el])
                nhb = napol = 0
                for j, g in enumerate(idx):
                    for m in g:
                        pe = prot[m].element.strip().upper()
                        dist = float(np.linalg.norm(L[j] - P[m]))
                        if pol_l[j] and pe in ("N", "O") and dist <= HBOND:
                            nhb += 1
                        elif apol_l[j] and pe in ("C", "S"):
                            napol += 1

                arom_res, stack, catpi = set(), [], []
                lrings = ligand_rings(heavy)
                for rid, rname, rc, rn_ in rings_prot:
                    if np.min(np.linalg.norm(L - rc, axis=1)) <= CONTACT + 1.5:
                        arom_res.add(f"{rname}{rid}")
                    for lc, ln in lrings:
                        dd = float(np.linalg.norm(lc - rc))
                        if dd <= STACK_CENTRE:
                            ang = np.degrees(np.arccos(
                                min(1.0, abs(float(np.dot(ln, rn_))))))
                            stack.append(f"{rname}{rid}({dd:.1f}A,{ang:.0f}deg)")
                    # cation-pi: a ligand N within reach of the ring face
                    for j, e in enumerate(lig_el):
                        if e == "N" and np.linalg.norm(L[j] - rc) <= CATPI:
                            catpi.append(f"{rname}{rid}")
                            break

                tag = f"{rn}#{k + 1}" if len(scored) > 1 else rn
                print(f"  {pid:22} {ch} {tag:10} {len(heavy):>3} atoms  "
                      f"depth {depth:5.1f} A  {site:7} "
                      f"apolar {napol:3d}  H-bond {nhb:2d}  "
                      f"aromatic {len(arom_res)}"
                      + (f"  stack {stack[0]}" if stack else "")
                      + (f"  cation-pi {catpi[0]}" if catpi else ""))
                rows.append([pid, LABEL.get(pid, ""), ch, rn, k + 1,
                             len(scored), len(heavy), fmt(depth), fmt(offset),
                             site, n_dbp, n_pbp, n_sw, len(touched),
                             napol, nhb,
                             fmt(napol / max(nhb + napol, 1) * 100, 1),
                             len(arom_res), ";".join(sorted(arom_res)),
                             ";".join(stack), ";".join(sorted(set(catpi)))])

    write_csv(os.path.join(TABLES, "ligand_environment.csv"),
              ["pdb", "description", "chain", "ligand", "ligand_index",
               "ligands_in_protomer", "heavy_atoms", "depth_from_entrance_A",
               "offset_from_channel_A", "site", "DBP_residues_contacted",
               "PBP_residues_contacted", "switch_loop_contacted",
               "residues_contacted", "apolar_contacts", "hbond_contacts",
               "percent_apolar", "aromatic_residues_engaged",
               "aromatic_residues", "pi_stacking", "cation_pi"], rows)

    print("\n  --- what the multi-ligand structure adds ---")
    multi = [r for r in rows if r[5] > 1]
    if multi:
        d = [float(r[7]) for r in multi]
        allres = set()
        for r in multi:
            allres |= set(r[18].split(";")) - {""}
        print(f"    {multi[0][0]} chain {multi[0][2]}: {len(multi)} ligands "
              f"spanning {min(d):.1f}-{max(d):.1f} A of the channel "
              f"({max(d) - min(d):.1f} A of path occupied at once)")
        for r in multi:
            print(f"      {r[3]}#{r[4]}  depth {float(r[7]):5.1f} A  "
                  f"{r[9]:7}  {r[13]:2d} residues  "
                  f"{float(r[16]):.0f}% apolar  aromatics: {r[18] or 'none'}")
        singles = [r for r in rows if r[5] == 1]
        sd = [float(r[7]) for r in singles]
        print(f"    single-ligand structures each occupy one point: "
              f"{min(sd):.1f}-{max(sd):.1f} A across all {len(singles)} "
              f"of them, but no one structure spans more than 0 A")
    print("\nwrote results/tables/ligand_environment.csv")


if __name__ == "__main__":
    main()
