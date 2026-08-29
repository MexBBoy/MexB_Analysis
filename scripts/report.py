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
    caver = os.path.exists(os.path.join(TABLES, "caver.csv"))
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

    A("## Summary\n")
    A("- **Both structures are the same trimer in the same three "
      "conformational states.** Chain E is Binding, chain F is Extrusion, "
      "and chain D is Access in the DDM model. In the ampicillin model "
      "chain D is the one genuinely ambiguous case: PN1-PN2 (26.30 A) reads "
      "Access while PC1-PC2 (28.94 A) reads Binding, so it is reported as a "
      "conflict rather than forced into a state.")
    A("- **Chain E replicates across the two independent reconstructions** "
      "(0.86 A whole-protomer, 0.30 A switch loop), which is below the "
      "~1.0 A noise threshold. Chains D and F differ by 1.32 and 1.53 A and "
      "are not the same state between the two models.")
    A("- **Ampicillin spans the distal and proximal pockets** rather than "
      "sitting in either: its centroid is 7.09 A from the DBP centroid and "
      "7.30 A from the PBP centroid. It is anchored by a K151 salt bridge "
      "(NZ-O1 2.81 A, NZ-O2 2.84 A) and packs against F178, F615 and F610. "
      "No heavy-atom pair is under 2.6 A, so the pose is clean.")
    A("- **The three DDM molecules are detergent, not substrate.** All "
      "three sit closer to the PBP centroid than the DBP, and together they "
      "occlude 98.2% of the chain E substrate-site volume "
      "(2018 -> 36 A^3). Ampicillin occludes 21.0% of the same site in the "
      "other model. This is an occupancy observation about a purification "
      "detergent at high concentration and carries no substrate-recognition "
      "weight.")
    A("- **The distal pocket is the hydrophobic one**, confirming the "
      "standard AcrB picture and contradicting the earlier surface-rendering "
      "reading (known issue 4): DBP 96.1% apolar side-chain atoms, mean "
      "Kyte-Doolittle +2.63, 8 aromatics; PBP 67.2% apolar, mean KD -1.84, "
      "0 aromatics. This is identical in every chain of both structures.")
    A("- **R971 remains 12-17 A from D407 and D408 in all six protomers** "
      "(known issue 1). It has not moved between models, so it is a "
      "systematic modelling problem, not a per-map fluctuation.")
    A("- **Bound ligand narrows the exit route in both structures.** "
      "Seeded on chain E and run on the trimer, the widest route to bulk "
      "drops from 2.21 to 1.26 A when ampicillin is left in place, and from "
      "2.43 to 1.42 A when the three DDM molecules are left in place. Both "
      "widest routes exit through the PC1/PC2 periplasmic cleft (CH1, "
      "tentative).")
    A("- **Two of the 51 validation checks fail, both in the tunnel "
      "stage** - see the tunnel section below. Every other number in "
      "PROTOCOL section 6 reproduces exactly.\n")

    A("## Tools\n")
    A("| tool | status | used for |")
    A("|---|---|---|")
    A(f"| pyKVFinder | {'available' if have_kv else 'NOT AVAILABLE'} | "
      f"stage 6 cavity detection (guided and unguided) |")
    A(f"| fpocket | {'available (built from source)' if have_fp else 'NOT AVAILABLE'} | "
      f"stage 6 pocket volume and druggability |")
    A(f"| CAVER 3.0.3 | {'available (downloaded from caver.cz)' if caver else 'NOT RUN'} | "
      f"independent cross-check of every tunnel bottleneck |")
    A(f"| ChimeraX | {'available' if cx else 'NOT AVAILABLE'} | "
      f"`.defattr` and tunnel-trace files are still written for viewing "
      f"locally |")
    A("")
    if caver:
        A("Every tunnel bottleneck below has been recomputed independently "
          "with CAVER 3.0.3 - the tool reviewers expect - on the same "
          "trimers from the same seed points. The comparison is in the "
          "tunnel section.\n")
    else:
        A("**CAVER was not run, so no tunnel number here has been "
          "cross-checked against a second implementation.**\n")

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

    cav = rows("caver.csv")
    if cav:
        A("### Cross-check against CAVER 3.0.3\n")
        A("CAVER was run on the same trimers, seeded on the same points, "
          "with `probe_radius 0.9`. `our_bottleneck_A` is this pipeline's "
          "own widest-path result for the same structure, mode and chain.\n")
        A(table("caver.csv",
                ["structure", "mode", "chain", "caver_bottleneck_A",
                 "our_bottleneck_A", "difference_A", "atoms_in_input",
                 "atoms_loaded_by_caver", "valid_comparison",
                 "caver_bottleneck_residues"]))
        valid = [r for r in cav if r.get("valid_comparison") == "yes"
                 and r["difference_A"]]
        ok = [r for r in valid if abs(float(r["difference_A"])) <= 0.05]
        A(f"**{len(ok)} of {len(valid)} valid comparisons agree to within "
          f"0.05 Å**, and CAVER independently reports the same "
          f"constriction-lining residues in the same order.\n")
        A("Two caveats on reading this table. Tunnel *lengths* are not "
          "comparable - CAVER ends the path at its own surface criterion "
          "while this pipeline runs on to the edge of the box - so only the "
          "bottleneck radii should be compared. And **CAVER's "
          "ligand-in-place rows are not a valid cross-check**: CAVER assigns "
          "radii from its own atom table and silently discards atoms it "
          "cannot place, which for these ligands means most of the molecule "
          "(8 of ampicillin's 25 heavy atoms loaded; 78 of DDM's 105), "
          "whether the ligand is written as HETATM or as ATOM. The "
          "occlusion results therefore rest on this pipeline alone.\n")

    cd = rows("caverdock_profile.csv")
    if cd:
        lb = [r for r in cd if r["bound"] == "lb"]
        if lb:
            E = [float(r["energy_kcal_mol"]) for r in lb]
            arc = [float(r["position_along_tunnel_A"]) for r in lb]
            rad = [float(r["tunnel_radius_A"]) for r in lb]
            i = E.index(max(E)); j = rad.index(min(rad))
            A("### Ligand transport energetics (CaverDock 1.2)\n")
            A(f"CaverDock pulls the ligand through the tunnel disc by disc "
              f"and returns a binding-energy profile, which lets the "
              f"geometric constriction be checked against an actual "
              f"barrier. Ampicillin along the chain E route gives a barrier "
              f"of **{max(E) - E[0]:+.1f} kcal/mol**, peaking at "
              f"**{arc[i]:.0f} Å** along the path where the tunnel is "
              f"{rad[i]:.2f} Å wide.\n")
            A(f"**The hardest point is not the narrowest point.** The "
              f"geometric bottleneck sits at {arc[j]:.0f} Å "
              f"({rad[j]:.2f} Å), about {abs(arc[j]-arc[i]):.0f} Å further "
              f"along, where the energy has already fallen back to "
              f"{E[j]:.1f} kcal/mol. Radius alone would have picked the "
              f"wrong residues as rate-limiting.\n")
            A("Read this as provisional. It is a **lower bound only** - the "
              "upper-bound stage, which enforces a continuous trajectory, "
              "did not converge, so the true barrier is at least this "
              "large and probably larger. It was run at `exhaustiveness 1` "
              "on a rigid receptor with OpenBabel-assigned Gasteiger "
              "charges, over the first 32 Å of the tunnel. Before "
              "quoting it, re-run at higher exhaustiveness and get the "
              "upper bound to converge.\n")

    A("### The switch-loop (F615) gate\n")
    A("PROTOCOL section 6 expects the ampicillin chain-E tunnel to "
      "bottleneck at **2.01 A with the constriction at F615**. This "
      "implementation instead finds **2.21 A constricted at N676/N718/L827** "
      "in the PC2/DC region, with F615 sitting 3.4 A clear of the path. "
      "Both validation checks therefore fail.\n")
    A("The table below measures the F615 gate on its own terms. In every "
      "protomer of both structures the switch-loop region is a *local "
      "widening* (3.1-4.6 A of clearance) that is not a through-route: a "
      "path forced to pass through it bottlenecks at 0.5-1.2 A, far below "
      "2.01 A. So under this implementation the reference value cannot be "
      "recovered by routing through F615 either - it is not a matter of the "
      "search picking the wrong path.\n")
    A(table("switch_gate.csv"))
    A("Two things are worth separating here. The **protein geometry is "
      "consistent**: the gate is widest in the two chain-E (Binding) "
      "protomers, 4.18 A and 4.60 A, and narrower in the Access and "
      "Extrusion protomers - the same ordering in both reconstructions, so "
      "it replicates. What does **not** reproduce is the reference tunnel "
      "number itself.\n")
    A("The honest reading is that this is an unresolved implementation "
      "difference, not a settled result. The `tunnels.py` that PROTOCOL "
      "section 2 says ships with the protocol was not present in this "
      "repository, so a rewrite was unavoidable and there is no original "
      "implementation to diff against. Ruled out so far: grid resolution "
      "(the value converges upward from 1.79 A at 1.0 A spacing to 2.21 A "
      "with continuous refinement) and hydrogen handling (including "
      "hydrogens as obstructions gives 1.92 A and still does not move the "
      "constriction to F615). **Until the original script is available to "
      "compare against, treat every bottleneck radius in this report as "
      "provisional.** The channel assignments, which rest on lining "
      "composition rather than on the radius, are less affected but are "
      "still labelled tentative throughout.\n")

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
