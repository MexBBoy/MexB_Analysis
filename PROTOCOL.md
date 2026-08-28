# MexB structural analysis pipeline — protocol

Hand this file to Claude Code at the root of the analysis repo. It describes
what to build, the constants to build it against, and the checks that must
pass before any number is believed.

Project context: cryo-EM structures of the *Pseudomonas aeruginosa* MexB
efflux transporter in substrate-bound states. Several structures exist or are
in progress; the pipeline must run identically over all of them so that
adding structure *n+1* is one command, not an afternoon.

---

## 0. Ground truth — do not re-derive these

**Sequence.** MexB is UniProt **P52002** (gene *mexB*, PA0426), 1046 residues.
Model residue numbering matches UniProt exactly with **zero offset**. This has
been verified against two structures with no sequence mismatches. Any new
model that does not satisfy this must be flagged, not silently renumbered.

**Domain ranges** (MexB numbering):

| domain | residues |
| --- | --- |
| porter | 40-130, 136-180, 278-281, 284-328, 571-666, 679-717, 814-818, 821-858 |
| docking / funnel | 181-277, 718-813 |
| TM (trimmed) | 10-35, 337-495, 516-565, 876-1030 |
| TM (untrimmed) | 10-35, 337-359, 361-565, 862-1030 |

The trimmed TM set drops the 496-515 cytoplasmic loop and the 859-875
junction. Both are highly mobile and inflate every RMSD they appear in. Use
the trimmed set for superposition; report the untrimmed one only if asked.

**Porter subdomains** (calibrated — reproduces previously published centroid
separations to within 0.4 A):

| subdomain | residues |
| --- | --- |
| PN1 | 40-130 |
| PN2 | 136-180, 278-281, 284-328 |
| PC1 | 571-666 |
| PC2 | 679-717, 814-818, 821-858 |

**Binding pockets:**

- Distal (DBP, hydrophobic trap): 136, 139, 178, 277, 279, 327, 573, 610, 612, 615, 617, 626, 628, 630
- Proximal (PBP): 79, 128, 151, 152, 176, 180, 273, 274, 276, 668, 672, 674, 676, 717, 819, 825, 828
- Switch loop: 613-622

The DBP set is the standard AcrB hydrophobic trap mapped onto MexB and is
reliable. **The PBP set is less certain** — it was assembled from AcrB
literature rather than taken from a MexB-specific source. Flag any conclusion
that depends sensitively on it.

**Proton relay:** D407, D408, K939, R971, T976. Measure minimum distance
between side-chain functional atoms: ASP OD1/OD2, LYS NZ, ARG NE/NH1/NH2,
THR OG1.

**Reference state diagnostics** (published MexB values, A):

| state | PN1-PN2 | PC1-PC2 |
| --- | --- | --- |
| Access | 26.5 | 26.2 |
| Binding | 28.3 | 29.1 |
| Extrusion | 29.8 | 24.5 |

---

## 1. Repository layout

```
structures/          input models, one PDB per structure (trimer, with ligands)
  Amp_MexB_<date>.pdb
  MexB_DDM_<n>_<date>.pdb
  ...
maps/                optional: sharpened maps and local-resolution maps, same stem
scripts/             analysis code (see section 2)
work/                intermediate files — safe to delete, never committed
results/
  tables/            CSV, one per analysis, all structures stacked
  figures/           PNG/SVG
  chimerax/          .defattr and pseudo-atom PDBs for viewing
  REPORT.md          generated summary, regenerated on every run
run.sh               regenerates everything from scratch
```

Rule: every file under `results/` must be reproducible by running `run.sh`
from a clean `work/`. Nothing hand-edited goes in there.

---

## 2. Scripts provided

Two working scripts ship with this protocol and should be moved into
`scripts/`:

- **`mexb_analysis.py`** — sliding-window per-residue Ca RMSD (with ChimeraX
  `.defattr` output), proton relay distance matrix, PC1-PC2 pseudo-contact
  counts. Subcommands `rmsd`, `relay`, `cleft`.
- **`tunnels.py`** — tunnel bottleneck radius by threshold connectivity, with
  geodesic path tracing and constriction residue identification. Writes a
  pseudo-atom tunnel trace with local radius in the B-factor column.

Extend these rather than rewriting them. If a rewrite is genuinely needed,
reproduce the validation numbers in section 6 first.

---

## 3. Environment

```bash
pip install numpy scipy pyKVFinder --break-system-packages
```

Optional, install if available and prefer it where it overlaps:

- **CAVER 3.0** (caver.cz, free academic licence) — the tunnel tool reviewers
  expect. If installed, run it alongside `tunnels.py` and report both.
- **fpocket** (github.com/Discngine/fpocket) — pocket volume and druggability.
  Build from source with `make`.
- **ChimeraX** — for headless figure generation via `chimerax --nogui --script`.

If a tool is unavailable, say so in the report rather than silently
substituting.

---

## 4. Stages

Run every stage over every structure and every chain. Results stack into one
table per analysis with `structure` and `chain` columns, so that adding a
structure adds rows rather than requiring new code.

### Stage 1 — Ingest and validate

For each input PDB:

- Parse chains, residue ranges, gaps, HETATM records, alternate locations.
- Translate the sequence and compare to P52002. Report mismatches and the
  numbering offset. **Halt on any mismatch or non-zero offset.**
- List ligands with residue name, chain, residue number, heavy-atom count.
- Report B-factors per ligand: min, median, max. Note whether they are
  per-atom or a single group value.

Output: `results/tables/inventory.csv`, plus a per-structure block in the
report.

### Stage 2 — State assignment

Per chain, compute:

- PN1-PN2 and PC1-PC2 centroid separations; nearest reference state and the
  deviation from it.
- PN1-PN2 and PC1-PC2 pseudo-contact counts (Ca-Ca within 10 A). Report raw
  count and count normalised by subdomain size.
- Proton relay distance matrix, all pairs.

Assign each chain to Access / Binding / Extrusion, and state which diagnostics
agree and which disagree. Do not force an assignment when diagnostics
conflict — report the conflict.

Output: `results/tables/states.csv`, `results/tables/proton_relay.csv`.

### Stage 3 — Inter-protomer conformational analysis

For each structure, all three protomer pairs, under two superposition frames
(fit on trimmed TM, and fit on porter):

- Whole-protomer Ca RMSD, and per-domain RMSD for porter, docking, TM.
- Sliding-window per-residue deviation, written as CSV and ChimeraX
  `.defattr`.
- Per-region summary (RMS and max) over: TM1, PN1, PN2, DN, TM2, Ialpha,
  TM3-6, cytoplasmic loop 496-515, TM6b, PC1, PC2, DC, junction 859-875,
  TM7-12.
- Rigid-body transform (rotation angle, axis, shift along axis) per domain.

The mechanistically important quantities are **TM2 and Ialpha displacement
between access and binding** and **the R2 / TM7-12 swing angle**. Report these
explicitly rather than leaving them buried in a table.

Output: `results/tables/rmsd_by_region.csv`, `results/chimerax/*.defattr`.

### Stage 4 — Cross-structure comparison

Same protomer, different structures. For each chain present in more than one
structure: whole-protomer RMSD, porter deviation after TM fit, TM deviation
after porter fit, switch-loop RMSD after porter fit.

The purpose is to distinguish real conformational differences from
reconstruction-to-reconstruction noise. Anything under ~1.0 A should be
treated as the same state.

Output: `results/tables/cross_structure.csv`.

### Stage 5 — Ligand environment

For every ligand in every structure:

- Contact residues at 4.5 A: minimum distance, contact count, count of
  contacts to aromatic ring atoms.
- Candidate hydrogen bonds: N/O to N/O within 3.6 A, listed atom-to-atom.
- Steric clashes: any heavy-atom pair under 2.6 A. **Any clash is a hard
  flag** — the pose or a side-chain rotamer needs revisiting before the
  structure is used for figures.
- Pocket assignment: distance from ligand centroid to DBP and PBP centroids,
  and a per-ligand-atom nearest-pocket call, so that a molecule spanning both
  is described as spanning rather than forced into one.
- Contacts to DBP side-chain atoms, as a count.

Then build the **cross-ligand contact matrix**: rows are residues, columns are
ligands across all structures, cells are minimum distance. Include a column
counting in how many ligands each residue appears. This is a key paper figure
and must regenerate automatically whenever a pose changes.

Output: `results/tables/contacts_<structure>_<ligand>.csv`,
`results/tables/contact_matrix.csv`.

### Stage 6 — Pocket and cavity analysis

Using pyKVFinder (and fpocket if available):

- Ligand-guided detection for each bound ligand: cavity volume, area, max and
  average depth, average hydropathy, lining residues.
- Unguided detection per protomer; identify the cavity containing the
  substrate site by DBP residue overlap; report the same metrics.
- Report whether PBP and DBP fall in one continuous cavity or two.
- Compute hydropathy separately for the PBP-lining and DBP-lining subsets.

Also compute a grid-based volume of the substrate site with each structure's
own ligands stripped, so that sites are comparable across structures, and
again with ligands present, to quantify occlusion.

Output: `results/tables/cavities.csv`, `results/tables/pocket_volumes.csv`.

### Stage 7 — Tunnels

Per protomer, seeded on the bound ligand where present and on the transferred
site position where absent (superpose on the pocket-lining residues to
transfer):

- Bottleneck radius from the site to bulk solvent.
- Geodesic path length and the coordinates of the narrowest point.
- Residues lining the constriction, with clearance for each.
- A pseudo-atom trace of the tunnel with local radius in the B-factor column.

Run twice for structures with bound detergent: once protein-only, once with
detergent as an obstruction. The difference quantifies occlusion of the exit
route.

**Run on the trimer, never on an isolated protomer.** An isolated protomer has
an artificially open protomer-protomer interface and the path finder escapes
through it, producing meaningless bottlenecks.

Where possible, assign each path to CH1, CH2 or CH3 by its exit location and
lining residues, and say so explicitly when the assignment is uncertain.

Output: `results/tables/tunnels.csv`, `results/chimerax/*_tunnel.pdb`.

### Stage 8 — Report

Generate `results/REPORT.md` containing every table, the flags raised, and a
short prose summary per stage. The report must state which tools ran, which
were unavailable, and which numbers changed since the previous run.

---

## 5. Known issues to check on every run

These are live problems, not settled facts. Re-check each run and report
status.

1. **R971 is 12-17 A from both D407 and D408 in every protomer of both
   structures analysed so far.** R971 is essential by mutagenesis (Guan &
   Nakae) and acts as the electrostatic switch in AcrB, so this is almost
   certainly a side-chain rotamer problem rather than biology. Report the
   distance every run and flag if it stays this large.

2. **Chain F proton relay is inconsistent between structures.** D407-K939 was
   4.24 A in one model and 2.82 A in another; D408-K939 was 9.26 versus 5.05.
   Treat chain F relay geometry as unresolved until checked against density.

3. **Chain F sits further from both other protomers than they do from each
   other.** State assignment and model quality for F both need verification.

4. **Hydrophobicity of the two pockets.** Atom-composition scoring gives the
   proximal pocket 67% apolar (mean Kyte-Doolittle -1.93, no aromatics) and
   the distal pocket 99% apolar (+2.63, eight aromatics) — the standard AcrB
   picture. An earlier reading from surface renderings suggested the reverse.
   Report the computed values and note the discrepancy until resolved.

5. **Ligand B-factors are group-refined** (a single value per molecule in the
   models seen so far). Say so whenever B-factors are reported.

6. **Poses change between model versions.** An earlier ampicillin pose had a
   1.74 A interpenetration with F610 and a different contact set. Never carry
   a contact list forward from a previous model version — always regenerate.

---

## 6. Validation — the pipeline must reproduce these

Before trusting any new or rewritten code, check it against these known
values. They come from the ampicillin-bound and DDM-bound structures.

Numbering and sequence: 0 mismatches against P52002, offset 0. Chain D
1-1032 complete; chain E 1-1030 with gaps at 229 and 360; chain F 1-1033 with
a gap at 663.

State diagnostics, ampicillin structure: PN1-PN2 D 26.30 / E 27.82 / F 29.96;
PC1-PC2 D 28.94 / E 28.38 / F 24.22.

Pseudo-contacts, ampicillin structure: PC1-PC2 D 5 / E 5 / F 16; PN1-PN2
D 21 / E 1 / F 4.

Proton relay, ampicillin structure, chain D: D407-K939 2.71, D408-K939 2.70.
Chain F: K939-T976 4.81.

Ampicillin contacts: K151 2.81 (NZ-O1) and 2.84 (NZ-O2); F178 3.25 with 37
contacts; F615 3.35; F610 3.76; no pair under 2.6 A.

Tunnel bottleneck, ampicillin trimer, chain E seeded on the ligand: 2.01 A,
constriction at F615.

Cross-structure chain E, ampicillin versus DDM: 0.86 A whole-protomer,
0.30 A switch loop.

If a rewritten implementation disagrees with any of these by more than
rounding, the implementation is wrong until proven otherwise.

---

## 7. Standards for reporting

- Distinguish what the data show from what they suggest. Conformational
  claims rest on protein geometry; ligand-based claims rest on a detergent or
  a substrate and are weaker.
- Detergent is not substrate. DDM is the purification detergent at high
  concentration and will occupy any groove that fits it. Keep protein-geometry
  arguments separate from detergent-occupancy arguments.
- Grid-based volumes computed here are internally comparable but are not
  drop-in replacements for fpocket or CASTp. Label them as such.
- A result that reproduces across two independent reconstructions is much
  stronger than one from a single map. Say explicitly when something
  replicates.
- Never fill a missing measurement with a plausible value. Report it as
  missing.

---

## 8. First run

```bash
bash run.sh
```

If `run.sh` does not exist yet, build it as part of the first run, then verify
against section 6 and write `results/REPORT.md`.

Priorities in order: stages 1 and 2 first, since everything downstream depends
on correct state assignment; then stage 5, since the contact matrix is the
most immediately useful paper output; then the rest.
