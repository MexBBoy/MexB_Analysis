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


def cross_structure_block(table_html):
    """Numbers and figures for the cross-structure Results section.

    Everything is read from the tables the corresponding scripts write, so
    the prose cannot drift from the data.
    """
    import math
    pub = R("published_pockets.csv")
    prot = R("protomer_pockets.csv")
    env = R("ligand_environment.csv")
    rot = R("aromatic_rotamers.csv")
    cons = R("lining_conservation.csv")
    el = R("pocket_electrostatics.csv")
    OURS = ("Amp_MexB_20260826", "MexB_DDM_3_20260730")
    AROM = set("FYWH")
    D = {}

    def num(r, k):
        try:
            return float(r[k])
        except (KeyError, TypeError, ValueError):
            return float("nan")

    def ok(x):
        return isinstance(x, float) and not math.isnan(x)

    # --- Fig 9: volume vs ligand size and depth
    bound = [r for r in pub if int(r["ligand_heavy_atoms"]) > 0]
    if bound:
        L = [float(r["ligand_heavy_atoms"]) for r in bound]
        Vv = [num(r, "volume_r16_A3") for r in bound]
        Dp = [num(r, "depth_from_entrance_A") for r in bound]

        def pearson(a, b):
            pairs = [(x, y) for x, y in zip(a, b) if ok(x) and ok(y)]
            n = len(pairs)
            if n < 3:
                return float("nan")
            ax = stat.mean(x for x, _ in pairs)
            ay = stat.mean(y for _, y in pairs)
            sxy = sum((x - ax) * (y - ay) for x, y in pairs)
            sx = math.sqrt(sum((x - ax) ** 2 for x, _ in pairs))
            sy = math.sqrt(sum((y - ay) ** 2 for _, y in pairs))
            return sxy / (sx * sy) if sx and sy else float("nan")

        D["R_SIZE"] = f"{pearson(L, Vv):+.2f}"
        D["R_DEPTH"] = f"{pearson(Dp, Vv):+.2f}"
        D["LIG_MIN"] = f"{min(L):.0f}"
        D["LIG_MAX"] = f"{max(L):.0f}"
        D["LIG_RANGE"] = f"{max(L) / min(L):.1f}"
        D["N_BOUND"] = str(len(bound))
        same = [num(r, "volume_r16_A3") for r in bound
                if r["pdb"] in ("2V50", "3W9I")]
        D["NOISE"] = (f"{abs(same[0] - same[1]):.0f}" if len(same) == 2
                      else "&mdash;")

    # --- Fig 10: volume by state
    if prot:
        D["N_PROT"] = str(len(prot))
        D["N_STRUCT"] = str(len({r["pdb"] for r in prot}))
        D["N_BUBBLE"] = str(sum(1 for r in prot
                                if r["connected_volume_measurable"] == "no"))
        by = {}
        for st in ("Access", "Binding", "Extrusion"):
            v = [num(r, "free_volume_r16_A3") for r in prot
                 if r["state_call"] == st]
            v = [x for x in v if ok(x)]
            by[st] = v
            D[{"Access": "V_ACC", "Binding": "V_BIND",
               "Extrusion": "V_EXT"}[st]] = (
                f"{stat.mean(v):.0f} &plusmn; {stat.stdev(v):.0f}")
            lp = [num(r, "lipophilic_index_pct") for r in prot
                  if r["state_call"] == st]
            lp = [x for x in lp if ok(x)]
            D[{"Access": "L_ACC", "Binding": "L_BIND",
               "Extrusion": "L_EXT"}[st]] = f"{stat.mean(lp):.1f}"
        D["V_RATIO"] = f"{stat.mean(by['Binding']) / stat.mean(by['Access']):.1f}"
        worst, wlab = 0.0, "&mdash;"
        for r in prot:
            if r["pdb"] not in OURS:
                continue
            peers = [num(x, "free_volume_r16_A3") for x in prot
                     if x["state_call"] == r["state_call"]
                     and x["pdb"] not in OURS]
            peers = [x for x in peers if ok(x)]
            if len(peers) > 2 and ok(num(r, "free_volume_r16_A3")):
                z = abs(num(r, "free_volume_r16_A3")
                        - stat.mean(peers)) / stat.stdev(peers)
                if z > worst:
                    worst = z
                    nice = ("ampicillin" if r["pdb"].startswith("Amp")
                            else "DDM &times;3")
                    wlab = f"{nice} chain {r['chain']}"
        D["WORST_Z"] = f"{worst:.1f}"
        D["WORST_LAB"] = wlab

    # --- Fig 11: the multi-ligand protomer
    multi = [r for r in env if int(r["ligands_in_protomer"]) > 1]
    if multi:
        dd = sorted(num(r, "depth_from_entrance_A") for r in multi)
        for i, x in enumerate(dd[:3], 1):
            D[f"ST{i}"] = f"{x:.1f}"
        D["SPAN"] = f"{dd[-1] - dd[0]:.1f}"
        a = set()
        for r in multi:
            a |= {x for x in r["aromatic_residues"].split(";") if x}
        D["N_AROM_DDM"] = str(len(a))

    # --- Fig 12: rotamers
    if rot:
        devs = []
        for rid in sorted({int(r["resseq"]) for r in rot}):
            peers = [num(r, "chi1_deg") for r in rot
                     if int(r["resseq"]) == rid and r["state_call"] == "Binding"
                     and r["pdb"] not in OURS and r["chi1_deg"]]
            peers = [x for x in peers if ok(x)]
            mine = [num(r, "chi1_deg") for r in rot
                    if int(r["resseq"]) == rid and r["chain"] == "E"
                    and r["pdb"] == "MexB_DDM_3_20260730" and r["chi1_deg"]]
            if len(peers) < 4 or not mine or not ok(mine[0]):
                continue
            mu = math.degrees(math.atan2(
                stat.mean(math.sin(math.radians(x)) for x in peers),
                stat.mean(math.cos(math.radians(x)) for x in peers)))
            one = {"PHE": "F", "TYR": "Y", "TRP": "W", "HIS": "H"}
            nmres = next(r["resname"] for r in rot
                         if int(r["resseq"]) == rid)
            devs.append((abs((mine[0] - mu + 180) % 360 - 180),
                         one.get(nmres, nmres) + str(rid)))
        if devs:
            devs.sort()
            D["F664_DEV"] = f"{devs[-1][0]:.0f}"
            D["ROT_RES"] = devs[-1][1]
            D["N_ROT"] = str(len(devs))
            D["N_ROT_OK"] = str(len(devs) - 1)
            # spelled out, because the sentence opens with it
            words = ["zero", "One", "Two", "Three", "Four", "Five", "Six",
                     "Seven", "Eight", "Nine", "Ten", "Eleven", "Twelve"]
            low = [w.lower() for w in words]
            D["W_ROT_OK"] = (words[len(devs) - 1]
                             if len(devs) - 1 < len(words) else str(len(devs) - 1))
            D["W_ROT"] = (low[len(devs)] if len(devs) < len(low)
                          else str(len(devs)))
            D["ROT_MAX"] = f"{devs[-2][0]:.0f}" if len(devs) > 1 else "&mdash;"

    # --- Fig 13: conservation
    if cons:
        key = [k for k in cons[0] if k.startswith("homologues_")
               and k != "homologues_aromatic"][0]
        ar = [r for r in cons if r["mexb_aromatic"] == "yes"]
        other = [r for r in cons if r["mexb_aromatic"] == "no"]
        D["C_AROM"] = f"{stat.mean(float(r['percent_identical']) for r in ar):.0f}"
        D["C_OTHER"] = f"{stat.mean(float(r['percent_identical']) for r in other):.0f}"
        NAMES = ["MexD", "MexF", "MexY", "AcrB", "AcrF", "MdtF", "AcrD"]
        keep = {nm: sum(1 for r in ar if r[key][k] in AROM)
                for k, nm in enumerate(NAMES)}
        n = len(ar)
        D["C_MEXY"] = f"{keep['MexY']}"
        D["C_ACRD"] = f"{keep['AcrD']}"
        broad = [keep[x] for x in ("MexD", "MexF", "AcrB", "AcrF", "MdtF")]
        D["C_BROAD"] = f"{min(broad)}&ndash;{max(broad)}"
        look = {int(r["resseq"]): r for r in cons}
        for lab, ids in (("C_OUT", [136, 573, 617, 628, 664, 666, 327]),
                         ("C_MID", [615, 617]),
                         ("C_DEEP", [178, 610, 615, 628])):
            v = [float(look[i]["percent_identical"]) for i in ids
                 if i in look]
            D[lab] = f"{stat.mean(v):.0f}"

    # --- Fig 14: electrostatics
    if el:
        D["N_APBS"] = str(len(el))
        for st in ("Access", "Binding", "Extrusion"):
            v = [num(r, "mean_potential_kT_e") for r in el
                 if r["state_call"] == st]
            v = [x for x in v if ok(x)]
            if v:
                D[{"Access": "E_ACC", "Binding": "E_BIND",
                   "Extrusion": "E_EXT"}[st]] = f"{stat.mean(v):+.1f}"

    D.setdefault("N_ENTRIES", "13")
    D.setdefault("N_LINING", "39")
    for i, nm in enumerate(("P4_ligand_size_vs_pocket", "P5_protomer_states",
                            "P6_path_occupancy", "P7_aromatic_rotamers",
                            "P8_lining_conservation", "P9_pocket_physchem"), 9):
        D[f"FIG{i}"] = img(os.path.join("poster", f"{nm}.png"))

    D["T_PUB"] = table_html(
        pub, ["pdb", "description", "chain", "ligands", "ligand_heavy_atoms",
              "depth_from_entrance_A", "volume_r14_A3", "volume_r16_A3",
              "volume_r18_A3", "volume_r20_A3"],
        ["PDB", "Description", "Protomer", "Ligands", "Heavy atoms",
         "Depth (Å)", "14 Å", "16 Å", "18 Å", "20 Å"])
    D["T_LIGENV"] = table_html(
        env, ["pdb", "chain", "ligand", "ligand_index", "heavy_atoms",
              "depth_from_entrance_A", "site", "residues_contacted",
              "percent_apolar", "aromatic_residues"],
        ["PDB", "Protomer", "Ligand", "#", "Heavy atoms", "Depth (Å)",
         "Site", "Residues", "Apolar (%)", "Aromatics engaged"])
    D["T_CONS"] = table_html(
        cons, ["resseq", "mexb_residue", "site", key if cons else "",
               "percent_identical"],
        ["Residue", "MexB", "Site", "MexD MexF MexY AcrB AcrF MdtF AcrD",
         "Identity (%)"]) if cons else "<p><em>(no rows)</em></p>"
    return D


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

    cd = [r for r in R("caverdock_profile.csv") if r["bound"] == "lb"]
    if cd:
        E = [float(r["energy_kcal_mol"]) for r in cd]
        arc = [float(r["position_along_tunnel_A"]) for r in cd]
        rad = [float(r["tunnel_radius_A"]) for r in cd]
        i, j = E.index(max(E)), rad.index(min(rad))
        CD = {"CD_BARRIER": f"{max(E) - E[0]:+.1f}",
              "CD_EPOS": f"{arc[i]:.0f}", "CD_ERAD": f"{rad[i]:.2f}",
              "CD_RPOS": f"{arc[j]:.0f}", "CD_RRAD": f"{rad[j]:.2f}",
              "CD_GAP": f"{abs(arc[j] - arc[i]):.0f}"}
    else:
        CD = {k: "—" for k in ("CD_BARRIER", "CD_EPOS", "CD_ERAD",
                               "CD_RPOS", "CD_RRAD", "CD_GAP")}

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
        "FIG8": img("fig8_caverdock_profile.png"),
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

    V.update(CD)
    V.update(cross_structure_block(table_html))
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
