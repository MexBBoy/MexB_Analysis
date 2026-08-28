#!/usr/bin/env python3
"""Build the data payload for the interactive 3D viewer.

Emits work/viewer_data.js: for each structure a trimmed PDB (backbone for
cartoon, full side chains for the pocket and switch-loop residues, plus the
ligands) and the tunnel traces as [x, y, z, radius] arrays.
"""
from __future__ import annotations

import csv
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import CXDIR, DBP, PBP, REPO, SWITCH_LOOP, TABLES, \
    WORK_DIR, load_structures

# 3Dmol builds its cartoon from CA (path) and O (ribbon orientation);
# N and C are dead weight here and roughly double the parse cost.
BACKBONE = {"CA", "O"}
KEEP_SIDECHAIN = set(DBP) | set(PBP) | set(SWITCH_LOOP)
SHORT = {"Amp_MexB_20260826": "amp", "MexB_DDM_3_20260730": "ddm"}


def trim_pdb(s):
    # HELIX/SHEET records let 3Dmol draw a proper cartoon from CA+O alone;
    # without them it can only draw a featureless trace.
    out = []
    for line in open(s.path):
        if line.startswith(("HELIX", "SHEET")):
            out.append(line.rstrip("\n"))
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[76:78].strip().upper() == "H":
            continue
        name = line[12:16].strip()
        if line.startswith("ATOM"):
            r = int(line[22:26])
            if name not in BACKBONE and r not in KEEP_SIDECHAIN:
                continue
        out.append(line.rstrip("\n"))
    return "\n".join(out) + "\nEND"


def traces(struct_name, mode):
    got = {}
    for f in sorted(glob.glob(os.path.join(
            CXDIR, f"{struct_name}_{mode}_*_tunnel.pdb"))):
        base = os.path.basename(f)
        m = re.match(rf"{re.escape(struct_name)}_{re.escape(mode)}_"
                     r"(?P<chain>[^_]+)_(?P<seed>.+)_t(?P<rank>\d+)"
                     r"_tunnel\.pdb$", base)
        if not m:
            continue
        chain, seed, rank = m["chain"], m["seed"], m["rank"]
        pts = []
        for line in open(f):
            if line.startswith("HETATM"):
                pts.append([round(float(line[30:38]), 2),
                            round(float(line[38:46]), 2),
                            round(float(line[46:54]), 2),
                            round(float(line[60:66]), 2)])
        got.setdefault(f"{chain}", []).append(
            {"chain": chain, "seed": seed, "rank": rank, "pts": pts})
    return got


def rows(name):
    p = os.path.join(TABLES, name)
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh))


def main():
    states = rows("states.csv")
    tun = rows("tunnels.csv")
    vol = rows("pocket_volumes.csv")
    ligs = rows("ligand_summary.csv")
    data = {}
    for s in load_structures():
        key = SHORT[s.name]
        st = {r["chain"]: r for r in states if r["structure"] == s.name}
        vv = {r["chain"]: r for r in vol if r["structure"] == s.name}
        data[key] = {
            "name": s.name,
            "pdb": trim_pdb(s),
            "ligands": sorted({n for (_, _, n, _) in s.ligands()}),
            "ligandResi": [[c, r, n] for (c, r, n, _) in s.ligands()],
            "tunnels": {"protein": traces(s.name, "protein"),
                        "withlig": traces(s.name, "withlig")},
            "states": {c: {"pn": r["PN1_PN2_sep"], "pc": r["PC1_PC2_sep"],
                           "call": r["state_call"],
                           "agree": r["diagnostics_agree"]}
                       for c, r in st.items()},
            "volumes": {c: {"free": r["volume_ligands_stripped_A3"],
                            "occ": r["volume_with_ligands_A3"],
                            "pct": r["occluded_pct"]}
                        for c, r in vv.items()},
            "tunnelRows": [r for r in tun if r["structure"] == s.name
                           and not r["note"]],
            "ligandRows": [r for r in ligs if r["structure"] == s.name],
        }
    out = os.path.join(WORK_DIR, "viewer_data.js")
    with open(out, "w") as fh:
        fh.write("const DATA = " + json.dumps(data, separators=(",", ":"))
                 + ";\n")
    mb = os.path.getsize(out) / 1e6
    print(f"wrote {out} ({mb:.2f} MB)")

    # inline the payload into the viewer page
    tpl = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "viewer_template.html")
    html = open(tpl).read()
    payload = open(out).read()
    assert "/*__DATA__*/" in html, "template placeholder missing"
    html = html.replace("/*__DATA__*/", payload)
    page = os.path.join(REPO, "results", "viewer", "mexb_tunnels.html")
    os.makedirs(os.path.dirname(page), exist_ok=True)
    with open(page, "w") as fh:
        fh.write(html)
    print(f"wrote {page} ({os.path.getsize(page)/1e6:.2f} MB)")
    for k, v in data.items():
        nt = sum(len(x) for x in v["tunnels"]["protein"].values())
        print(f"  {k}: {len(v['pdb'].splitlines())} atom lines, "
              f"{nt} protein-mode tunnels, ligands {v['ligands']}")


if __name__ == "__main__":
    main()
