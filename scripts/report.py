#!/usr/bin/env python3
"""Stage 8 - assemble results/REPORT.md from the generated tables."""
from __future__ import annotations

import csv
import datetime
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mexb_common import REPO, TABLES, WORK_DIR

PREV = os.path.join(WORK_DIR, "prev_tables")
OUT = os.path.join(REPO, "results", "REPORT.md")


def rows(name):
    p = os.path.join(TABLES, name)
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh))


def table(name, cols=None, where=None, limit=None):
    rs = rows(name)
    if where:
        rs = [r for r in rs if where(r)]
    if not rs:
        return "_(no rows)_\n"
    cols = cols or list(rs[0].keys())
    out = ["| " + " | ".join(cols) + " |",
           "|" + "|".join("---" for _ in cols) + "|"]
    for r in (rs[:limit] if limit else rs):
        out.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
    if limit and len(rs) > limit:
        out.append(f"\n_{len(rs) - limit} further rows in "
                   f"`results/tables/{name}`._")
    return "\n".join(out) + "\n"


def changed_since_previous():
    """Which tables changed since the last run."""
    if not os.path.isdir(PREV):
        return None
    diffs = []
    for f in sorted(os.listdir(TABLES)):
        if not f.endswith(".csv"):
            continue
        a, b = os.path.join(TABLES, f), os.path.join(PREV, f)
        if not os.path.exists(b):
            diffs.append((f, "new"))
            continue
        if open(a).read() != open(b).read():
            ra, rb = rows(f), []
            with open(b) as fh:
                rb = list(csv.DictReader(fh))
            n = sum(1 for x, y in zip(ra, rb) if x != y)
            extra = abs(len(ra) - len(rb))
            diffs.append((f, f"{n} row(s) changed"
                             + (f", {extra} row(s) added/removed"
                                if extra else "")))
    for f in sorted(os.listdir(PREV)):
        if f.endswith(".csv") and not os.path.exists(
                os.path.join(TABLES, f)):
            diffs.append((f, "removed"))
    return diffs


def flags():
    out = []
    for f in sorted(os.listdir(TABLES)):
        if f.startswith("flags_") and f.endswith(".txt"):
            for line in open(os.path.join(TABLES, f)):
                line = line.strip()
                if line:
                    out.append(line)
    seen, uniq = set(), []
    for f in out:
        if f not in seen:
            seen.add(f)
            uniq.append(f)
    return uniq


def tool_status():
    fp = os.path.join(WORK_DIR, "fpocket", "bin", "fpocket")
    have_fp = os.path.exists(fp)
    import importlib.util
    have_kv = importlib.util.find_spec("pyKVFinder") is not None
    caver = any(os.path.exists(os.path.join(p, "caver"))
                for p in os.environ.get("PATH", "").split(":") if p)
    cx = any(os.path.exists(os.path.join(p, b))
             for p in os.environ.get("PATH", "").split(":") if p
             for b in ("chimerax", "ChimeraX"))
    return have_kv, have_fp, caver, cx


def main():
    have_kv, have_fp, caver, cx = tool_status()
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    try:
        rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=REPO, capture_output=True,
                             text=True).stdout.strip() or "uncommitted"
    except Exception:
        rev = "unknown"

    val = rows("validation.csv")
    npass = sum(1 for r in val if r["status"] == "PASS")
    fails = [r for r in val if r["status"] == "FAIL"]
    fl = flags()
    diffs = changed_since_previous()

    L = []
    A = L.append
    A(f"# MexB substrate-bound structures - tunnel and pocket analysis\n")
    A(f"Generated {now} from commit `{rev}` by `bash run.sh`. "
      f"Every table under `results/` is reproducible from a clean `work/`; "
      f"nothing here is hand-edited.\n")

    A("## Tools\n")
    A("| tool | status | used for |")
    A("|---|---|---|")
    A(f"| pyKVFinder | {'available' if have_kv else 'NOT AVAILABLE'} | "
      f"stage 6 cavity detection (guided and unguided) |")
    A(f"| fpocket | {'available (built from source)' if have_fp else 'NOT AVAILABLE'} | "
      f"stage 6 pocket volume and druggability |")
    A(f"| CAVER 3.0 | NOT AVAILABLE | the academic build is behind a "
      f"registration wall and could not be fetched in this environment; "
      f"tunnels come from `scripts/tunnels.py` only |")
    A(f"| ChimeraX | {'available' if cx else 'NOT AVAILABLE'} | "
      f"`.defattr` and tunnel-trace files are still written for viewing "
      f"locally |")
    A("")
    A("Because CAVER could not be run, **every tunnel number in this report "
      "comes from a single implementation and has not been cross-checked "
      "against the tool reviewers expect.** That is the largest single "
      "caveat here.\n")

    A("## Validation against PROTOCOL section 6\n")
    A(f"**{npass}/{len(val)} checks pass.**\n")
    if fails:
        A("Failing checks:\n")
        A("| check | got | expected |")
        A("|---|---|---|")
        for r in fails:
            A(f"| {r['check']} | {r['value']} | {r['expected']} |")
        A("")
    A(table("validation.csv", ["check", "value", "expected", "status"]))

    A("## Flags raised this run\n")
    if fl:
        for f in fl:
            A(f"- {f}")
    else:
        A("_None._")
    A("")

    A("## Stage 1 - ingest and validation\n")
    A(table("inventory.csv"))
    A("Ligand inventory:\n")
    A(table("ligand_inventory.csv"))
    A("Ligand B-factors are group-refined - a single value per molecule - "
      "in both models, so they carry no per-atom information and must be "
      "described as such wherever quoted (known issue 5).\n")

    A("## Stage 2 - state assignment\n")
    A(table("states.csv"))
    A(table("proton_relay.csv", limit=20))

    A("## Stage 3 - inter-protomer comparison\n")
    A(table("rmsd_protomer_pairs.csv"))
    A("Key mechanical quantities (TM-frame fit):\n")
    A(table("rmsd_by_region.csv",
            where=lambda r: (r["fit_frame"] == "TM_trimmed"
                             and r["region"] in ("TM2", "Ialpha", "TM7-12"))))
    A("Rigid-body domain motions:\n")
    A(table("rigid_body.csv",
            where=lambda r: r["fit_frame"] == "TM_trimmed"))

    A("## Stage 4 - cross-structure comparison\n")
    A(table("cross_structure.csv"))

    A("## Stage 5 - ligand environment\n")
    A(table("ligand_summary.csv"))
    A("Cross-ligand contact matrix (residues contacted by more than one "
      "ligand):\n")
    A(table("contact_matrix.csv",
            where=lambda r: int(r["n_ligands"]) > 1))

    A("## Stage 6 - pockets and cavities\n")
    A(table("cavities.csv"))
    A("Grid-based substrate-site volumes. **These are internally comparable "
      "across these structures only and are not drop-in replacements for "
      "fpocket or CASTp volumes.** The 'stripped' column removes the "
      "structure's own ligands, including the three DDM molecules, so the "
      "sites can be compared like for like; the difference is occlusion.\n")
    A(table("pocket_volumes.csv"))
    if os.path.exists(os.path.join(TABLES, "fpocket.csv")):
        A("fpocket, on the ligand-stripped trimer:\n")
        A(table("fpocket.csv", limit=18))
    A("Pocket composition (known issue 4):\n")
    A(table("pocket_composition.csv", limit=12))

    A("## Stage 7 - tunnels\n")
    A("Run on the trimer in every case. Seeds are the bound ligand centroid "
      "where the chain has a ligand, otherwise the DBP/PBP midpoint of that "
      "chain. `protein` mode strips all ligands; `withlig` mode keeps them "
      "as obstructions, so the difference measures occlusion of the exit "
      "route.\n")
    A(table("tunnels.csv",
            ["structure", "mode", "chain", "seed", "tunnel_rank",
             "bottleneck_radius_A", "geodesic_path_length_A",
             "constriction_lining_clearance_A", "channel_call",
             "assignment_confidence"]))

    A("## Changes since the previous run\n")
    if diffs is None:
        A("_No previous run to compare against._")
    elif not diffs:
        A("No table changed.")
    else:
        for f, what in diffs:
            A(f"- `{f}`: {what}")
    A("")
    open(OUT, "w").write("\n".join(L))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
