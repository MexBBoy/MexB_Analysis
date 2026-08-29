#!/usr/bin/env python3
"""Check the pipeline against the known values in PROTOCOL.md section 6.

Every check is re-derived from the generated tables, so this fails loudly if
a code change moves a number that is supposed to be fixed.
"""
from __future__ import annotations

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import TABLES

AMP = "Amp_MexB_20260826"
DDM = "MexB_DDM_3_20260730"


def load(name):
    p = os.path.join(TABLES, name)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return list(csv.DictReader(fh))


def close(a, b, tol=0.011):
    try:
        return abs(float(a) - float(b)) <= tol
    except (TypeError, ValueError):
        return False


CHECKS = []

# Constants demoted by map validation (2026-08-30). The ampicillin pose is
# only supported at the 26th-34th percentile of its own map, so these still
# have to reproduce - they guard against silent code drift - but they are no
# longer evidence that the pose is right. PROTOCOL known issue 6.
PROVISIONAL = "provisional - pose weakly supported by density (issue 6)"


def check(label, got, want, ok=None, note=""):
    if ok is None:
        ok = close(got, want)
    CHECKS.append((label, got, want, ok, note))


def run():
    inv = load("inventory.csv") or []
    for r in inv:
        check(f"{r['structure']} {r['chain']} mismatches vs P52002",
              r["seq_mismatches_P52002"], "0",
              r["seq_mismatches_P52002"] == "0")
        check(f"{r['structure']} {r['chain']} numbering offset",
              r["numbering_offset"], "0", r["numbering_offset"] == "0")
    want_gaps = {("D", "none"), ("E", "229;360"), ("F", "663")}
    for r in inv:
        exp = dict(want_gaps)[r["chain"]]
        check(f"{r['structure']} chain {r['chain']} gaps", r["gaps"], exp,
              r["gaps"] == exp)
    want_range = {"D": ("1", "1032"), "E": ("1", "1030"), "F": ("1", "1033")}
    for r in inv:
        lo, hi = want_range[r["chain"]]
        check(f"{r['structure']} chain {r['chain']} residue range",
              f"{r['first_res']}-{r['last_res']}", f"{lo}-{hi}",
              r["first_res"] == lo and r["last_res"] == hi)

    st = load("states.csv") or []
    exp_pn = {"D": 26.30, "E": 27.82, "F": 29.96}
    exp_pc = {"D": 28.94, "E": 28.38, "F": 24.22}
    exp_cpn = {"D": "21", "E": "1", "F": "4"}
    exp_cpc = {"D": "5", "E": "5", "F": "16"}
    for r in st:
        if r["structure"] != AMP:
            continue
        c = r["chain"]
        check(f"amp {c} PN1-PN2 separation", r["PN1_PN2_sep"], exp_pn[c])
        check(f"amp {c} PC1-PC2 separation", r["PC1_PC2_sep"], exp_pc[c])
        check(f"amp {c} PN1-PN2 pseudo-contacts",
              r["PN1_PN2_contacts10A"], exp_cpn[c],
              r["PN1_PN2_contacts10A"] == exp_cpn[c])
        check(f"amp {c} PC1-PC2 pseudo-contacts",
              r["PC1_PC2_contacts10A"], exp_cpc[c],
              r["PC1_PC2_contacts10A"] == exp_cpc[c])

    relay = load("proton_relay.csv") or []

    def rel(struct, ch, a, b):
        for r in relay:
            if (r["structure"] == struct and r["chain"] == ch
                    and {r["res1"], r["res2"]} == {a, b}):
                return r["min_dist_A"]
        return None
    check("amp D ASP407-LYS939", rel(AMP, "D", "ASP407", "LYS939"), 2.71)
    check("amp D ASP408-LYS939", rel(AMP, "D", "ASP408", "LYS939"), 2.70)
    check("amp F LYS939-THR976", rel(AMP, "F", "LYS939", "THR976"), 4.81,
          note="chain F of this map is not supported by density; value "
               "reproduces but is not determinable (issue 2)")

    con = load(f"contacts_{AMP}_ZZ7E2000.csv") or []

    def cres(n):
        for r in con:
            if r["resseq"] == str(n):
                return r
        return {}
    check("amp K151 min distance", cres(151).get("min_dist_A"), 2.81,
          note=PROVISIONAL)
    check("amp F178 min distance", cres(178).get("min_dist_A"), 3.25,
          note=PROVISIONAL)
    check("amp F178 contact count", cres(178).get("n_contacts_4.5A"), "37",
          cres(178).get("n_contacts_4.5A") == "37", note=PROVISIONAL)
    check("amp F615 min distance", cres(615).get("min_dist_A"), 3.35,
          note=PROVISIONAL)
    check("amp F610 min distance", cres(610).get("min_dist_A"), 3.76,
          note=PROVISIONAL)
    hb = load(f"hbonds_{AMP}_ZZ7E2000.csv") or []
    got = {(r["lig_atom"], r["prot_atom"]): r["dist_A"] for r in hb
           if r["resseq"] == "151"}
    check("amp K151 NZ-O1", got.get(("O1", "NZ")), 2.81, note=PROVISIONAL)
    check("amp K151 NZ-O2", got.get(("O2", "NZ")), 2.84, note=PROVISIONAL)
    ls = load("ligand_summary.csv") or []
    amp_row = next((r for r in ls if r["structure"] == AMP), {})
    check("amp: no heavy-atom pair under 2.6 A",
          amp_row.get("n_close_contacts_2.6A"), "0",
          amp_row.get("n_close_contacts_2.6A") == "0")

    cs = load("cross_structure.csv") or []
    e = next((r for r in cs if r["chain"] == "E"), {})
    check("cross-structure chain E whole protomer",
          e.get("whole_protomer_rmsd_allCA_fit_A"), 0.86)
    check("cross-structure chain E switch loop",
          e.get("switch_loop_rmsd_A"), 0.30)

    tun = load("tunnels.csv") or []
    t = next((r for r in tun if r["structure"] == AMP
              and r["mode"] == "protein" and r["chain"] == "E"
              and r.get("tunnel_rank") == "1"), {})
    check("amp chain E tunnel bottleneck (seeded on ligand)",
          t.get("bottleneck_radius_A"), 2.01)
    check("amp chain E constriction residue is F615",
          (t.get("constriction_lining_clearance_A") or "?").split(":")[0],
          "PHE615E",
          (t.get("constriction_lining_clearance_A") or "").startswith(
              "PHE615E"))

    # Independent cross-check: CAVER 3.0.3 against our own tunnel code.
    cav = load("caver.csv") or []
    for r in cav:
        if not r["caver_bottleneck_A"] or not r["our_bottleneck_A"]:
            continue
        # CAVER discards most ligand atoms, so its ligand-in-place runs are
        # not a valid comparison and are not checked here.
        if r.get("valid_comparison", "yes") != "yes":
            continue
        lab = (f"CAVER agrees with ours: {r['structure']} {r['mode']} "
               f"chain {r['chain']}")
        check(lab, r["caver_bottleneck_A"], r["our_bottleneck_A"],
              close(r["caver_bottleneck_A"], r["our_bottleneck_A"], 0.05))


    npass = sum(1 for c in CHECKS if c[3])
    nprov = sum(1 for c in CHECKS if c[4])
    print(f"PROTOCOL section 6 validation: {npass}/{len(CHECKS)} pass"
          f"{f', {nprov} provisional' if nprov else ''}\n")
    for label, got, want, ok, note in CHECKS:
        if not ok:
            print(f"  FAIL  {label}: got {got!r}, expected {want!r}")
    print()
    for label, got, want, ok, note in CHECKS:
        if ok:
            print(f"  pass  {label}: {got}"
                  + (f"   [{note}]" if note else ""))
    with open(os.path.join(TABLES, "validation.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["check", "value", "expected", "status", "note"])
        for label, got, want, ok, note in CHECKS:
            w.writerow([label, got, want, "PASS" if ok else "FAIL", note])
    return npass, len(CHECKS)


if __name__ == "__main__":
    p, n = run()
    sys.exit(0 if p == n else 1)
