#!/usr/bin/env python3
"""MexB analysis driver: stages 1-5.

Subcommands
  inventory  stage 1  ingest, sequence validation, ligand inventory
  state      stage 2  centroid separations, pseudo-contacts, state call
  relay      stage 2  proton relay distance matrix
  cleft      stage 2  PC1-PC2 / PN1-PN2 pseudo-contact counts only
  rmsd       stage 3  inter-protomer RMSD, per-region, rigid-body, .defattr
  cross      stage 4  cross-structure comparison
  contacts   stage 5  ligand environment and cross-ligand contact matrix
  all        every stage above
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.spatial import cKDTree, distance as spdist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import (  # noqa: E402
    AA3to1, AROMATIC, CXDIR, DBP, DOMAINS, KD, PBP, REFERENCE_STATES,
    REGIONS, RELAY, RELAY_ATOMS, SUBDOMAINS, SWITCH_LOOP, TABLES,
    DETERGENTS, Structure, apply_rt, centroid, coords,
    density_warning, fmt,
    load_reference_sequence, load_structures, rmsd, rotation_angle_axis,
    superpose_on, write_csv,
)

FLAGS: list[str] = []


def flag(msg):
    FLAGS.append(msg)
    print("  [FLAG] " + msg)


# ------------------------------------------------------------------ stage 1

def stage_inventory(structs):
    print("\n=== Stage 1: ingest and validate ===")
    ref = load_reference_sequence()
    rows, lrows = [], []
    for s in structs:
        print(f"\n{s.name}")
        for ch in s.chains:
            seq = s.sequence(ch)
            resids = sorted(seq)
            lo, hi = resids[0], resids[-1]
            gaps = [r for r in range(lo, hi + 1) if r not in seq]
            # sequence check at zero offset against UniProt P52002
            mism = []
            for r, aa in seq.items():
                if r < 1 or r > len(ref):
                    mism.append((r, aa, "-"))
                elif ref[r - 1] != aa:
                    mism.append((r, aa, ref[r - 1]))
            # best offset scan, to report if zero is not optimal
            best_off, best_score = 0, -1
            for off in range(-5, 6):
                sc = sum(1 for r, aa in seq.items()
                         if 1 <= r + off <= len(ref) and ref[r + off - 1] == aa)
                if sc > best_score:
                    best_off, best_score = off, sc
            print(f"  chain {ch}: {lo}-{hi}  n={len(seq)}  "
                  f"gaps={gaps if gaps else 'none'}  mismatches={len(mism)}  "
                  f"offset={best_off}")
            if mism:
                flag(f"{s.name} chain {ch}: {len(mism)} sequence mismatches "
                     f"vs P52002 (first: {mism[:3]})")
            if best_off != 0:
                flag(f"{s.name} chain {ch}: non-zero numbering offset "
                     f"{best_off} vs P52002")
            rows.append([s.name, ch, lo, hi, len(seq),
                         ";".join(map(str, gaps)) or "none",
                         len(mism), best_off])
        for (lch, lres, lname, ats) in s.ligands():
            heavy = [a for a in ats if not a.is_hydrogen]
            b = np.array([a.bfac for a in heavy])
            group = "group" if np.allclose(b, b[0]) else "per-atom"
            print(f"  ligand {lname} {lch}{lres}: {len(heavy)} heavy atoms, "
                  f"B {b.min():.2f}/{np.median(b):.2f}/{b.max():.2f} ({group})")
            if group == "group":
                flag(f"{s.name} {lname} {lch}{lres}: B-factors are "
                     f"group-refined (single value {b[0]:.2f}) - known issue 5")
            lrows.append([s.name, lname, lch, lres, len(heavy), len(ats),
                          fmt(b.min()), fmt(float(np.median(b))),
                          fmt(b.max()), group,
                          "detergent" if lname in DETERGENTS else "substrate"])

    write_csv(os.path.join(TABLES, "inventory.csv"),
              ["structure", "chain", "first_res", "last_res", "n_res",
               "gaps", "seq_mismatches_P52002", "numbering_offset"], rows)
    write_csv(os.path.join(TABLES, "ligand_inventory.csv"),
              ["structure", "ligand", "chain", "resseq", "n_heavy", "n_atoms",
               "b_min", "b_median", "b_max", "b_mode", "class"], lrows)
    return rows, lrows


# ------------------------------------------------------------------ stage 2

def cleft_metrics(s, ch, cutoff=10.0):
    ca = s.ca(ch)
    out = {}
    for a, b in (("PN1", "PN2"), ("PC1", "PC2")):
        A = np.array([ca[r] for r in SUBDOMAINS[a] if r in ca])
        B = np.array([ca[r] for r in SUBDOMAINS[b] if r in ca])
        d = spdist.cdist(A, B)
        n = int((d < cutoff).sum())
        out[f"{a}-{b}"] = {
            "sep": float(np.linalg.norm(A.mean(0) - B.mean(0))),
            "contacts": n,
            "contacts_norm": n / (len(A) + len(B)),
            "nA": len(A), "nB": len(B),
        }
    return out


def assign_state(pn, pc):
    scores = {k: abs(pn - a) + abs(pc - b)
              for k, (a, b) in REFERENCE_STATES.items()}
    best = min(scores, key=scores.get)
    pn_best = min(REFERENCE_STATES, key=lambda k: abs(pn - REFERENCE_STATES[k][0]))
    pc_best = min(REFERENCE_STATES, key=lambda k: abs(pc - REFERENCE_STATES[k][1]))
    return best, scores[best], pn_best, pc_best


def stage_state(structs):
    print("\n=== Stage 2: state assignment ===")
    rows = []
    for s in structs:
        for ch in s.chains:
            m = cleft_metrics(s, ch)
            pn = m["PN1-PN2"]["sep"]
            pc = m["PC1-PC2"]["sep"]
            call, dev, pn_best, pc_best = assign_state(pn, pc)
            agree = "yes" if pn_best == pc_best else "no"
            print(f"  {s.name} {ch}: PN1-PN2 {pn:6.2f}  PC1-PC2 {pc:6.2f}  "
                  f"-> {call} (L1 dev {dev:.2f}); PN says {pn_best}, "
                  f"PC says {pc_best}, agree={agree}")
            print(f"      pseudo-contacts PN1-PN2 {m['PN1-PN2']['contacts']} "
                  f"({m['PN1-PN2']['contacts_norm']:.3f}/res), "
                  f"PC1-PC2 {m['PC1-PC2']['contacts']} "
                  f"({m['PC1-PC2']['contacts_norm']:.3f}/res)")
            dw = density_warning(s.name, ch)
            if dw:
                flag(f"{s.name} chain {ch}: NOT SUPPORTED BY DENSITY - {dw}. "
                     f"State assignment for this protomer is not "
                     f"determinable; exclude it from downstream analysis "
                     f"(PROTOCOL known issues 2/3)")
            if agree == "no":
                flag(f"{s.name} chain {ch}: cleft diagnostics conflict - "
                     f"PN1-PN2 nearest {pn_best}, PC1-PC2 nearest {pc_best}; "
                     f"reporting conflict rather than forcing a call")
            rows.append([s.name, ch, fmt(pn), fmt(pc),
                         m["PN1-PN2"]["contacts"],
                         fmt(m["PN1-PN2"]["contacts_norm"], 3),
                         m["PC1-PC2"]["contacts"],
                         fmt(m["PC1-PC2"]["contacts_norm"], 3),
                         call, fmt(dev), pn_best, pc_best, agree,
                         "no - not supported by density" if dw else "yes"])
    write_csv(os.path.join(TABLES, "states.csv"),
              ["structure", "chain", "PN1_PN2_sep", "PC1_PC2_sep",
               "PN1_PN2_contacts10A", "PN1_PN2_contacts_per_res",
               "PC1_PC2_contacts10A", "PC1_PC2_contacts_per_res",
               "state_call", "L1_deviation", "PN_nearest", "PC_nearest",
               "diagnostics_agree", "density_supported"], rows)
    return rows


def relay_atoms_for(s, ch, resid):
    want = RELAY[resid]
    ats = s.residue_atoms(ch, resid)
    if not ats:
        return None, []
    obs = ats[0].resname
    names = RELAY_ATOMS.get(obs, ())
    sel = [a for a in ats if a.name in names]
    return obs, sel


def stage_relay(structs):
    print("\n=== Stage 2: proton relay ===")
    rows = []
    keys = sorted(RELAY)
    for s in structs:
        for ch in s.chains:
            vals = {}
            for i, r1 in enumerate(keys):
                n1, a1 = relay_atoms_for(s, ch, r1)
                for r2 in keys[i + 1:]:
                    n2, a2 = relay_atoms_for(s, ch, r2)
                    if not a1 or not a2:
                        d = None
                        note = "missing side-chain atoms"
                    else:
                        d = float(spdist.cdist(coords(a1), coords(a2)).min())
                        note = ""
                    lab = f"{n1 or '?'}{r1}-{n2 or '?'}{r2}"
                    vals[lab] = d
                    rows.append([s.name, ch, f"{n1 or '?'}{r1}",
                                 f"{n2 or '?'}{r2}", fmt(d), note,
                                 "no - not supported by density"
                                 if density_warning(s.name, ch) else "yes"])
            txt = "  ".join(f"{k} {fmt(v)}" for k, v in vals.items())
            print(f"  {s.name} {ch}: {txt}")
            dw = density_warning(s.name, ch)
            if dw:
                flag(f"{s.name} chain {ch}: relay distances NOT "
                     f"DETERMINABLE - {dw}. This is the source of the "
                     f"known chain F relay inconsistency; exclude these "
                     f"values (PROTOCOL known issue 2)")
            # known issue 1
            for lab in ("ASP407-ARG971", "ASP408-ARG971"):
                d = vals.get(lab)
                if d is not None and d > 10.0:
                    flag(f"{s.name} chain {ch}: {lab} = {d:.2f} A "
                         f"(known issue 1, R971 rotamer) - still unresolved")
    write_csv(os.path.join(TABLES, "proton_relay.csv"),
              ["structure", "chain", "res1", "res2", "min_dist_A", "note",
               "density_supported"], rows)
    return rows


# ------------------------------------------------------------------ stage 3

def per_residue_dev(mob, mch, ref, rch, R, t):
    mca, rca = mob.ca(mch), ref.ca(rch)
    common = sorted(set(mca) & set(rca))
    M = apply_rt(R, t, np.array([mca[r] for r in common]))
    T = np.array([rca[r] for r in common])
    d = np.linalg.norm(M - T, axis=1)
    return common, d


def sliding(resids, dev, win=9):
    half = win // 2
    idx = {r: i for i, r in enumerate(resids)}
    out = []
    for i, r in enumerate(resids):
        vals = [dev[idx[q]] for q in range(r - half, r + half + 1) if q in idx]
        out.append(float(np.sqrt(np.mean(np.square(vals)))))
    return out


def write_defattr(path, name, chain, resids, values):
    with open(path, "w") as fh:
        fh.write(f"attribute: {name}\nmatch mode: any\nrecipient: residues\n")
        for r, v in zip(resids, values):
            fh.write(f"\t/{chain}:{r}\t{v:.3f}\n")


def region_summary(resids, dev):
    m = dict(zip(resids, dev))
    out = {}
    for reg, rr in REGIONS.items():
        v = [m[r] for r in rr if r in m]
        out[reg] = (float(np.sqrt(np.mean(np.square(v)))), float(max(v)),
                    len(v)) if v else (None, None, 0)
    return out


def domain_rigid_body(mob, mch, ref, rch, R, t, resids):
    """Residual rigid-body transform of one domain after the global fit."""
    mca, rca = mob.ca(mch), ref.ca(rch)
    common = [r for r in resids if r in mca and r in rca]
    if len(common) < 4:
        return None
    M = apply_rt(R, t, np.array([mca[r] for r in common]))
    T = np.array([rca[r] for r in common])
    from mexb_common import kabsch
    R2, t2 = kabsch(M, T)
    ang, axis = rotation_angle_axis(R2)
    shift_vec = T.mean(0) - M.mean(0)
    return {
        "n": len(common),
        "rmsd_before": rmsd(M, T),
        "rmsd_after": rmsd(apply_rt(R2, t2, M), T),
        "angle": ang,
        "axis": axis,
        "shift": float(np.linalg.norm(shift_vec)),
        "shift_along_axis": float(np.dot(shift_vec, axis)),
    }


def stage_rmsd(structs):
    print("\n=== Stage 3: inter-protomer conformational analysis ===")
    rows, reg_rows, rb_rows = [], [], []
    for s in structs:
        chs = s.chains
        pairs = [(chs[i], chs[j]) for i in range(len(chs))
                 for j in range(i + 1, len(chs))]
        for frame, fit_res in (("TM_trimmed", DOMAINS["TM"]),
                               ("porter", DOMAINS["porter"])):
            for a, b in pairs:
                R, t, n, fitr = superpose_on(s, a, s, b, fit_res)
                resids, dev = per_residue_dev(s, a, s, b, R, t)
                whole = float(np.sqrt(np.mean(np.square(dev))))
                dom = {}
                for dname in ("porter", "docking", "TM"):
                    m = dict(zip(resids, dev))
                    v = [m[r] for r in DOMAINS[dname] if r in m]
                    dom[dname] = float(np.sqrt(np.mean(np.square(v))))
                print(f"  {s.name} {a}->{b} fit={frame} (n={n}, fit "
                      f"{fitr:.2f} A): whole {whole:.2f}  porter "
                      f"{dom['porter']:.2f}  docking {dom['docking']:.2f}  "
                      f"TM {dom['TM']:.2f}")
                rows.append([s.name, frame, a, b, n, fmt(fitr), fmt(whole),
                             fmt(dom["porter"]), fmt(dom["docking"]),
                             fmt(dom["TM"])])
                sw = sliding(resids, dev)
                write_defattr(
                    os.path.join(CXDIR,
                                 f"{s.name}_{a}vs{b}_{frame}.defattr"),
                    f"dev{a}{b}{frame}".replace("_", ""), a, resids, sw)
                write_csv(os.path.join(
                    TABLES, f"per_residue_{s.name}_{a}vs{b}_{frame}.csv"),
                    ["resseq", "deviation_A", "sliding_rms_A"],
                    [[r, fmt(d, 3), fmt(w, 3)]
                     for r, d, w in zip(resids, dev, sw)])
                for reg, (rms, mx, nn) in region_summary(resids, dev).items():
                    reg_rows.append([s.name, frame, a, b, reg, nn,
                                     fmt(rms), fmt(mx)])
                for dname in ("porter", "docking", "TM"):
                    rb = domain_rigid_body(s, a, s, b, R, t, DOMAINS[dname])
                    if rb:
                        rb_rows.append([
                            s.name, frame, a, b, dname, rb["n"],
                            fmt(rb["rmsd_before"]), fmt(rb["rmsd_after"]),
                            fmt(rb["angle"]),
                            ",".join(f"{x:.3f}" for x in rb["axis"]),
                            fmt(rb["shift"]), fmt(rb["shift_along_axis"])])
    write_csv(os.path.join(TABLES, "rmsd_protomer_pairs.csv"),
              ["structure", "fit_frame", "chain_A", "chain_B", "n_fit",
               "fit_rmsd_A", "whole_protomer_rmsd_A", "porter_rmsd_A",
               "docking_rmsd_A", "TM_rmsd_A"], rows)
    write_csv(os.path.join(TABLES, "rmsd_by_region.csv"),
              ["structure", "fit_frame", "chain_A", "chain_B", "region",
               "n_res", "rms_dev_A", "max_dev_A"], reg_rows)
    write_csv(os.path.join(TABLES, "rigid_body.csv"),
              ["structure", "fit_frame", "chain_A", "chain_B", "domain",
               "n_res", "rmsd_before_A", "rmsd_after_A", "rotation_deg",
               "axis", "centroid_shift_A", "shift_along_axis_A"], rb_rows)

    # highlight the mechanistically important quantities
    print("\n  -- key mechanical quantities (TM-frame fit) --")
    for r in reg_rows:
        if r[1] == "TM_trimmed" and r[4] in ("TM2", "Ialpha", "TM7-12"):
            print(f"     {r[0]} {r[2]}vs{r[3]} {r[4]:>8}: rms {r[6]} A, "
                  f"max {r[7]} A")
    for r in rb_rows:
        if r[1] == "TM_trimmed" and r[4] == "porter":
            print(f"     {r[0]} {r[2]}vs{r[3]} porter rigid-body swing: "
                  f"{r[8]} deg, shift {r[10]} A")
    return rows, reg_rows, rb_rows


# ------------------------------------------------------------------ stage 4

def stage_cross(structs):
    print("\n=== Stage 4: cross-structure comparison ===")
    rows = []
    for i in range(len(structs)):
        for j in range(i + 1, len(structs)):
            A, B = structs[i], structs[j]
            for ch in sorted(set(A.chains) & set(B.chains)):
                # Whole-protomer RMSD is an all-CA best fit (this is the
                # definition that reproduces the protocol's validation value
                # of 0.86 A for chain E; a TM-frame fit gives 1.13 A).
                allres = sorted(set(A.ca(ch)) & set(B.ca(ch)))
                _, _, n, whole = superpose_on(A, ch, B, ch, allres)
                # domain deviations under the two named frames
                R, t, _, _ = superpose_on(A, ch, B, ch, DOMAINS["TM"])
                resids, dev = per_residue_dev(A, ch, B, ch, R, t)
                whole_tm = float(np.sqrt(np.mean(np.square(dev))))
                m = dict(zip(resids, dev))
                porter_after_tm = float(np.sqrt(np.mean(np.square(
                    [m[r] for r in DOMAINS["porter"] if r in m]))))
                # porter frame
                R2, t2, n2, _ = superpose_on(A, ch, B, ch, DOMAINS["porter"])
                resids2, dev2 = per_residue_dev(A, ch, B, ch, R2, t2)
                m2 = dict(zip(resids2, dev2))
                tm_after_porter = float(np.sqrt(np.mean(np.square(
                    [m2[r] for r in DOMAINS["TM"] if r in m2]))))
                switch = float(np.sqrt(np.mean(np.square(
                    [m2[r] for r in SWITCH_LOOP if r in m2]))))
                same = "same state (<1.0 A)" if whole < 1.0 else "differs"
                print(f"  {A.name} vs {B.name} chain {ch}: whole "
                      f"{whole:.2f} (TM-frame {whole_tm:.2f})  porter|TM-fit "
                      f"{porter_after_tm:.2f}  TM|porter-fit "
                      f"{tm_after_porter:.2f}  switch-loop "
                      f"{switch:.2f}  -> {same}")
                rows.append([A.name, B.name, ch, n, fmt(whole), fmt(whole_tm),
                             fmt(porter_after_tm), fmt(tm_after_porter),
                             fmt(switch), same])
    write_csv(os.path.join(TABLES, "cross_structure.csv"),
              ["structure_A", "structure_B", "chain", "n_fit",
               "whole_protomer_rmsd_allCA_fit_A",
               "whole_protomer_rmsd_TM_fit_A", "porter_dev_after_TM_fit_A",
               "TM_dev_after_porter_fit_A", "switch_loop_rmsd_A",
               "interpretation"], rows)
    return rows


# ------------------------------------------------------------------ stage 5

def ligand_contacts(s, lig, cutoff=4.5):
    lch, lres, lname, lats = lig
    heavy = [a for a in lats if not a.is_hydrogen]
    L = coords(heavy)
    prot = [a for a in s.protein_atoms if not a.is_hydrogen]
    P = coords(prot)
    tree = cKDTree(P)
    pairs = tree.query_ball_point(L, r=cutoff)
    per_res = {}
    for li, idxs in enumerate(pairs):
        for pi in idxs:
            pa = prot[pi]
            d = float(np.linalg.norm(L[li] - P[pi]))
            key = (pa.chain, pa.resseq, pa.resname)
            e = per_res.setdefault(key, {"min": 1e9, "n": 0, "arom": 0,
                                         "minpair": None, "dbp_sc": 0})
            e["n"] += 1
            if pa.name in AROMATIC.get(pa.resname, ()):
                e["arom"] += 1
            if (pa.resseq in DBP and pa.name not in
                    ("N", "CA", "C", "O")):
                e["dbp_sc"] += 1
            if d < e["min"]:
                e["min"] = d
                e["minpair"] = (heavy[li].name, pa.name)
    return heavy, prot, per_res


def stage_contacts(structs):
    print("\n=== Stage 5: ligand environment ===")
    matrix = {}
    ligcols = []
    summary = []
    for s in structs:
        for lig in s.ligands():
            lch, lres, lname, lats = lig
            tag = f"{s.name}:{lname}{lch}{lres}"
            ligcols.append(tag)
            heavy, prot, per_res = ligand_contacts(s, lig)
            L = coords(heavy)
            lcen = L.mean(0)
            # pocket centroids from this chain's own lining residues
            ca = s.ca(lch)
            dbp_c = centroid(ca, DBP)
            pbp_c = centroid(ca, PBP)
            d_dbp = float(np.linalg.norm(lcen - dbp_c))
            d_pbp = float(np.linalg.norm(lcen - pbp_c))
            # per-atom nearest pocket, by min distance to pocket side chains
            def pocket_atoms(resids):
                out = []
                for r in resids:
                    out += [a for a in s.residue_atoms(lch, r)
                            if a.name not in ("N", "C", "O")]
                return coords(out) if out else np.zeros((0, 3))
            DA, PA = pocket_atoms(DBP), pocket_atoms(PBP)
            # Primary per-atom call: nearer pocket centroid. This is
            # consistent with the reported centroid distances. The
            # nearest-lining-atom variant is reported alongside because the
            # two disagree for elongated ligands (a DDM tail can touch a DBP
            # side chain while its bulk sits in the PBP).
            n_dbp = n_pbp = 0
            n_dbp_lin = 0
            for p in L:
                if np.linalg.norm(p - dbp_c) < np.linalg.norm(p - pbp_c):
                    n_dbp += 1
                else:
                    n_pbp += 1
                dd = np.linalg.norm(DA - p, axis=1).min() if len(DA) else 1e9
                dp = np.linalg.norm(PA - p, axis=1).min() if len(PA) else 1e9
                if dd < dp:
                    n_dbp_lin += 1
            frac = n_dbp / len(L)
            if 0.2 < frac < 0.8:
                pocket_call = "spans DBP and PBP"
            elif frac >= 0.8:
                pocket_call = "distal (DBP)"
            else:
                pocket_call = "proximal (PBP)"
            # hydrogen bonds
            hb = []
            for la in heavy:
                if la.element not in ("N", "O"):
                    continue
                for pa in prot:
                    if pa.element not in ("N", "O"):
                        continue
                    d = float(np.linalg.norm(np.array(la.xyz) -
                                             np.array(pa.xyz)))
                    if d <= 3.6:
                        hb.append((la.name, pa.chain, pa.resname, pa.resseq,
                                   pa.name, d))
            hb.sort(key=lambda x: x[-1])
            # clashes
            # Steric clashes. Protocol: any heavy-atom pair under 2.6 A is a
            # hard flag. We additionally record whether both partners are
            # N/O, because such a pair at ~2.6 A is a short hydrogen bond
            # rather than a steric overlap - the flag still stands, but the
            # distinction changes what needs fixing.
            clash = []
            for la in heavy:
                for pa in prot:
                    d = float(np.linalg.norm(np.array(la.xyz) -
                                             np.array(pa.xyz)))
                    if d < 2.6:
                        polar = (la.element in ("N", "O") and
                                 pa.element in ("N", "O"))
                        clash.append((la.name, pa.chain, pa.resname,
                                      pa.resseq, pa.name, d,
                                      "short H-bond (both N/O)" if polar
                                      else "steric overlap"))
            dbp_sc = sum(e["dbp_sc"] for e in per_res.values())
            print(f"\n  {tag}: {len(heavy)} heavy atoms, "
                  f"{len(per_res)} contact residues at 4.5 A")
            print(f"    centroid to DBP {d_dbp:.2f} A, to PBP {d_pbp:.2f} A "
                  f"-> {pocket_call} ({n_dbp}/{len(L)} atoms nearer the DBP "
                  f"centroid; {n_dbp_lin}/{len(L)} by nearest lining atom)")
            print(f"    contacts to DBP side-chain atoms: {dbp_sc}")
            top = sorted(per_res.items(), key=lambda kv: kv[1]["min"])[:8]
            for (c, r, rn), e in top:
                print(f"      {rn}{r}{c}: min {e['min']:.2f} A, "
                      f"{e['n']} contacts, {e['arom']} aromatic")
            if clash:
                for c in clash:
                    flag(f"{tag}: CLOSE CONTACT {c[0]}-{c[2]}{c[3]}{c[1]}."
                         f"{c[4]} = {c[5]:.2f} A (<2.6) [{c[6]}] - pose or "
                         f"rotamer must be revisited before figures")
            else:
                print("    no heavy-atom pair under 2.6 A")
            rows = []
            for (c, r, rn), e in sorted(per_res.items(),
                                        key=lambda kv: kv[1]["min"]):
                rows.append([s.name, lname, f"{lch}{lres}", c, r, rn,
                             fmt(e["min"]), e["n"], e["arom"], e["dbp_sc"],
                             f"{e['minpair'][0]}-{e['minpair'][1]}",
                             "DBP" if r in DBP else
                             ("PBP" if r in PBP else "")])
                matrix.setdefault((r, rn), {})[tag] = e["min"]
            write_csv(os.path.join(
                TABLES, f"contacts_{s.name}_{lname}{lch}{lres}.csv"),
                ["structure", "ligand", "lig_id", "chain", "resseq",
                 "resname", "min_dist_A", "n_contacts_4.5A",
                 "n_aromatic_contacts", "n_DBP_sidechain_contacts",
                 "closest_atom_pair", "pocket"], rows)
            write_csv(os.path.join(
                TABLES, f"hbonds_{s.name}_{lname}{lch}{lres}.csv"),
                ["lig_atom", "chain", "resname", "resseq", "prot_atom",
                 "dist_A"],
                [[a, c, rn, r, pa, fmt(d)] for a, c, rn, r, pa, d in hb])
            write_csv(os.path.join(
                TABLES, f"clashes_{s.name}_{lname}{lch}{lres}.csv"),
                ["lig_atom", "chain", "resname", "resseq", "prot_atom",
                 "dist_A", "interpretation"],
                [[a, c, rn, r, pa, fmt(d), k]
                 for a, c, rn, r, pa, d, k in clash])
            summary.append([s.name, lname, f"{lch}{lres}", len(heavy),
                            len(per_res), fmt(d_dbp), fmt(d_pbp),
                            n_dbp, n_pbp, n_dbp_lin, pocket_call, dbp_sc,
                            len(hb), len(clash),
                            sum(1 for c in clash if c[6] == "steric overlap"),
                            "detergent" if lname in DETERGENTS
                            else "substrate"])

    # cross-ligand contact matrix
    resrows = []
    for (r, rn), d in sorted(matrix.items()):
        cells = [fmt(d.get(t)) for t in ligcols]
        resrows.append([r, rn, "DBP" if r in DBP else
                        ("PBP" if r in PBP else ""),
                        "switch" if r in SWITCH_LOOP else ""] +
                       cells + [len(d)])
    write_csv(os.path.join(TABLES, "contact_matrix.csv"),
              ["resseq", "resname", "pocket", "switch_loop"] + ligcols +
              ["n_ligands"], resrows)
    write_csv(os.path.join(TABLES, "ligand_summary.csv"),
              ["structure", "ligand", "lig_id", "n_heavy",
               "n_contact_residues", "centroid_to_DBP_A",
               "centroid_to_PBP_A", "atoms_nearer_DBP_centroid",
               "atoms_nearer_PBP_centroid", "atoms_nearer_DBP_lining_atom",
               "pocket_call", "n_DBP_sidechain_contacts", "n_hbond_cands",
               "n_close_contacts_2.6A", "n_true_steric_overlaps",
               "class"], summary)
    shared = [r for r in resrows if r[-1] > 1]
    print(f"\n  contact matrix: {len(resrows)} residues over "
          f"{len(ligcols)} ligands; {len(shared)} contacted by more than one")
    return summary, resrows, ligcols


# ---------------------------------------------------------------- pocket chem

def pocket_composition(structs):
    """Known issue 4: atom-composition hydrophobicity of PBP vs DBP."""
    print("\n=== Known issue 4: pocket hydrophobicity ===")
    rows = []
    for s in structs:
        for ch in s.chains:
            for pname, resids in (("PBP", PBP), ("DBP", DBP)):
                apolar = polar = 0
                arom = 0
                kd = []
                for r in resids:
                    ats = s.residue_atoms(ch, r)
                    if not ats:
                        continue
                    rn = ats[0].resname
                    kd.append(KD.get(AA3to1.get(rn, "X"), 0.0))
                    if rn in AROMATIC and rn != "HIS":
                        arom += 1
                    for a in ats:
                        if a.name in ("N", "CA", "C", "O"):
                            continue
                        if a.element == "C":
                            apolar += 1
                        elif a.element in ("N", "O", "S"):
                            polar += 1
                tot = apolar + polar
                pct = 100.0 * apolar / tot if tot else float("nan")
                rows.append([s.name, ch, pname, len(kd), tot, apolar, polar,
                             fmt(pct, 1), fmt(float(np.mean(kd))), arom])
                print(f"  {s.name} {ch} {pname}: {pct:.1f}% apolar "
                      f"side-chain atoms, mean KD {np.mean(kd):+.2f}, "
                      f"{arom} aromatics")
    write_csv(os.path.join(TABLES, "pocket_composition.csv"),
              ["structure", "chain", "pocket", "n_res", "n_sidechain_atoms",
               "n_apolar_C", "n_polar_NOS", "pct_apolar", "mean_KD",
               "n_aromatic_res"], rows)
    return rows


# ---------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["inventory", "state", "relay", "cleft",
                                      "rmsd", "cross", "contacts",
                                      "composition", "all"])
    a = ap.parse_args()
    structs = load_structures()
    print(f"Loaded {len(structs)} structures: "
          f"{', '.join(s.name for s in structs)}")
    run = a.stage
    if run in ("inventory", "all"):
        stage_inventory(structs)
    if run in ("state", "cleft", "all"):
        stage_state(structs)
    if run in ("relay", "all"):
        stage_relay(structs)
    if run in ("rmsd", "all"):
        stage_rmsd(structs)
    if run in ("cross", "all"):
        stage_cross(structs)
    if run in ("contacts", "all"):
        stage_contacts(structs)
    if run in ("composition", "all"):
        pocket_composition(structs)
    if FLAGS:
        with open(os.path.join(TABLES, f"flags_{run}.txt"), "w") as fh:
            fh.write("\n".join(FLAGS) + "\n")
        print(f"\n{len(FLAGS)} flags raised (see "
              f"results/tables/flags_{run}.txt)")


if __name__ == "__main__":
    main()
