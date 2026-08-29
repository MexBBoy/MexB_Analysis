#!/usr/bin/env python3
"""Assemble results/paper/mexb_paper.html - the analysis written up in the
format of the Nature Communications article the protocol was modelled on
(Lawrence et al. 2025, MdtF).

Every number in the prose is substituted from results/tables at build time,
so the manuscript cannot drift from the data. Figures are embedded as
data URIs so the page is self-contained.
"""
from __future__ import annotations

import base64
import csv
import datetime
import os
import statistics as stat
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import FIGURES, REPO, TABLES

OUT = os.path.join(REPO, "results", "paper", "mexb_paper.html")
TPL = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "paper_template.html")


def R(name):
    p = os.path.join(TABLES, name)
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh))


def img(name):
    p = os.path.join(FIGURES, name)
    if not os.path.exists(p):
        return ""
    b = base64.b64encode(open(p, "rb").read()).decode()
    return f"data:image/png;base64,{b}"


def table_html(rows, cols, headers=None, cls="tbl"):
    if not rows:
        return "<p><em>(no rows)</em></p>"
    headers = headers or cols
    h = "".join(f"<th>{x}</th>" for x in headers)
    body = ""
    for r in rows:
        body += "<tr>" + "".join(
            f"<td>{r.get(c,'')}</td>" for c in cols) + "</tr>"
    return (f'<div class="twrap"><table class="{cls}"><thead><tr>{h}</tr>'
            f'</thead><tbody>{body}</tbody></table></div>')


def main():
    states = {(r["structure"], r["chain"]): r for r in R("states.csv")}
    cross = {r["chain"]: r for r in R("cross_structure.csv")}
    ligs = R("ligand_summary.csv")
    vols = {(r["structure"], r["chain"]): r for r in R("pocket_volumes.csv")}
    cav = R("cavities.csv")
    comp = R("pocket_composition.csv")
    tun = R("tunnels.csv")
    cavr = R("caver.csv")
    gate = R("switch_gate.csv")
    relay = R("proton_relay.csv")
    inv = R("inventory.csv")
    val = R("validation.csv")
    cm = R("contact_matrix.csv")

    AMP, DDM = "Amp_MexB_20260826", "MexB_DDM_3_20260730"

    def t1(struct, mode, ch):
        for r in tun:
            if (r["structure"] == struct and r["mode"] == mode
                    and r["chain"] == ch and r["tunnel_rank"] == "1"):
                return r
        return {}

    def pocket_stat(pocket, key):
        v = [float(r[key]) for r in comp if r["pocket"] == pocket]
        return stat.mean(v)

    amp_p, amp_w = t1(AMP, "protein", "E"), t1(AMP, "withlig", "E")
    ddm_p, ddm_w = t1(DDM, "protein", "E"), t1(DDM, "withlig", "E")
    ampv, ddmv = vols[(AMP, "E")], vols[(DDM, "E")]
    ampl = next(r for r in ligs if r["ligand"] == "ZZ7")

    cav_ok = [r for r in cavr if r["difference_A"]]
    cav_prot = [r for r in cavr if r["mode"] == "protein"]
    cav_agree = [r for r in cav_prot
                 if r["difference_A"] and abs(float(r["difference_A"])) <= .05]

    npass = sum(1 for r in val if r["status"] == "PASS")

    V = {
        "DATE": datetime.datetime.now().strftime("%d %B %Y"),
        # states
        "AMP_D_PN": states[(AMP, "D")]["PN1_PN2_sep"],
        "AMP_D_PC": states[(AMP, "D")]["PC1_PC2_sep"],
        "AMP_E_PN": states[(AMP, "E")]["PN1_PN2_sep"],
        "AMP_E_PC": states[(AMP, "E")]["PC1_PC2_sep"],
        "AMP_F_PN": states[(AMP, "F")]["PN1_PN2_sep"],
        "AMP_F_PC": states[(AMP, "F")]["PC1_PC2_sep"],
        "DDM_D_PN": states[(DDM, "D")]["PN1_PN2_sep"],
        "DDM_D_PC": states[(DDM, "D")]["PC1_PC2_sep"],
        "DDM_E_PN": states[(DDM, "E")]["PN1_PN2_sep"],
        "DDM_E_PC": states[(DDM, "E")]["PC1_PC2_sep"],
        "DDM_F_PN": states[(DDM, "F")]["PN1_PN2_sep"],
        "DDM_F_PC": states[(DDM, "F")]["PC1_PC2_sep"],
        # cross structure
        "X_D": cross["D"]["whole_protomer_rmsd_allCA_fit_A"],
        "X_E": cross["E"]["whole_protomer_rmsd_allCA_fit_A"],
        "X_F": cross["F"]["whole_protomer_rmsd_allCA_fit_A"],
        "X_E_SW": cross["E"]["switch_loop_rmsd_A"],
        "X_D_SW": cross["D"]["switch_loop_rmsd_A"],
        "X_F_SW": cross["F"]["switch_loop_rmsd_A"],
        # ligand
        "AMP_NRES": ampl["n_contact_residues"],
        "AMP_DBP": ampl["centroid_to_DBP_A"],
        "AMP_PBP": ampl["centroid_to_PBP_A"],
        "AMP_HB": ampl["n_hbond_cands"],
        # pockets
        "PBP_APOLAR": f"{pocket_stat('PBP','pct_apolar'):.1f}",
        "DBP_APOLAR": f"{pocket_stat('DBP','pct_apolar'):.1f}",
        "PBP_KD": f"{pocket_stat('PBP','mean_KD'):+.2f}",
        "DBP_KD": f"{pocket_stat('DBP','mean_KD'):+.2f}",
        # volumes
        "AMP_VFREE": ampv["volume_ligands_stripped_A3"],
        "AMP_VOCC": ampv["volume_with_ligands_A3"],
        "AMP_PCT": ampv["occluded_pct"],
        "DDM_VFREE": ddmv["volume_ligands_stripped_A3"],
        "DDM_VOCC": ddmv["volume_with_ligands_A3"],
        "DDM_PCT": ddmv["occluded_pct"],
        # tunnels
        "AMP_BN": amp_p.get("bottleneck_radius_A", "—"),
        "AMP_BN_W": amp_w.get("bottleneck_radius_A", "—"),
        "DDM_BN": ddm_p.get("bottleneck_radius_A", "—"),
        "DDM_BN_W": ddm_w.get("bottleneck_radius_A", "—"),
        "AMP_CONS": amp_p.get("constriction_lining_clearance_A", ""),
        "AMP_CHAN": amp_p.get("channel_call", ""),
        "N_MATRIX": str(len(cm)),
        "N_SHARED": str(sum(1 for r in cm if int(r["n_ligands"]) > 1)),
        "NPASS": str(npass), "NCHECK": str(len(val)),
        "CAVER_AGREE": str(len(cav_agree)), "CAVER_N": str(len(cav_prot)),
        # figures
        "FIG1": img("fig1_state_diagnostics.png"),
        "FIG2": img("fig2_tunnel_profiles.png"),
        "FIG3": img("fig3_site_occlusion.png"),
        "FIG4": img("fig4_contact_matrix.png"),
        "FIG5": img("fig5_pocket_hydrophobicity.png"),
        "FIG6": img("fig6_switch_gate.png"),
        "FIG7": img("fig7_per_residue_EvsF.png"),
        # tables
        "T_STATES": table_html(
            R("states.csv"),
            ["structure", "chain", "PN1_PN2_sep", "PC1_PC2_sep",
             "PN1_PN2_contacts10A", "PC1_PC2_contacts10A", "state_call",
             "diagnostics_agree"],
            ["Structure", "Protomer", "PN1–PN2 (Å)", "PC1–PC2 (Å)",
             "PN1–PN2 contacts", "PC1–PC2 contacts", "State",
             "Diagnostics agree"]),
        "T_TUNNELS": table_html(
            [r for r in tun if r["tunnel_rank"] == "1" and not r["note"]],
            ["structure", "mode", "chain", "bottleneck_radius_A",
             "geodesic_path_length_A", "channel_call",
             "constriction_lining_clearance_A"],
            ["Structure", "Ligand", "Protomer", "Bottleneck (Å)",
             "Path length (Å)", "Channel", "Constriction lining"]),
        "T_CAVER": table_html(
            cavr, ["structure", "mode", "chain", "caver_bottleneck_A",
                   "our_bottleneck_A", "difference_A",
                   "caver_bottleneck_residues"],
            ["Structure", "Ligand", "Protomer", "CAVER (Å)",
             "This work (Å)", "Δ (Å)", "CAVER bottleneck residues"]),
        "T_VOL": table_html(
            R("pocket_volumes.csv"),
            ["structure", "chain", "volume_ligands_stripped_A3",
             "volume_with_ligands_A3", "occluded_volume_A3", "occluded_pct"],
            ["Structure", "Protomer", "Ligand-free (Å³)",
             "With ligand (Å³)", "Occluded (Å³)", "Occluded (%)"]),
        "T_GATE": table_html(
            gate, ["structure", "chain", "gate_widest_clearance_A",
                   "bottleneck_forced_through_gate_A"],
            ["Structure", "Protomer", "Widest clearance in gate (Å)",
             "Bottleneck forced through gate (Å)"]),
        "T_INV": table_html(
            inv, ["structure", "chain", "first_res", "last_res", "gaps",
                  "seq_mismatches_P52002", "numbering_offset"],
            ["Structure", "Chain", "First", "Last", "Gaps",
             "Mismatches vs P52002", "Offset"]),
        "T_RELAY": table_html(
            [r for r in relay if r["chain"] == "E"],
            ["structure", "res1", "res2", "min_dist_A"],
            ["Structure", "Residue 1", "Residue 2", "Min. distance (Å)"]),
        "T_CAV": table_html(
            [r for r in cav if r["detection"].startswith("unguided")],
            ["structure", "chain", "volume_A3", "area_A2", "max_depth_A",
             "avg_hydropathy", "DBP_overlap", "PBP_overlap", "continuity"],
            ["Structure", "Protomer", "Volume (Å³)", "Area (Å²)",
             "Max depth (Å)", "Avg hydropathy", "DBP overlap",
             "PBP overlap", "Continuity"]),
    }

    html = open(TPL, encoding="utf-8").read()
    for k, v in V.items():
        html = html.replace(f"§{k}§", str(v))
    left = [x for x in html.split("§") if x and x.isupper() and len(x) < 22]
    if left:
        print(f"  WARNING unfilled placeholders: {sorted(set(left))[:8]}")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    open(OUT, "w", encoding="utf-8").write(html)
    print(f"wrote {OUT} ({os.path.getsize(OUT)/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
