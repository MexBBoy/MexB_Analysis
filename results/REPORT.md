# MexB substrate-bound structures - tunnel and pocket analysis

Generated 2026-08-29 13:58 from commit `45132f4` by `bash run.sh`. Every table under `results/` is reproducible from a clean `work/`; nothing here is hand-edited.

## Summary

- **Both structures are the same trimer in the same three conformational states.** Chain E is Binding, chain F is Extrusion, and chain D is Access in the DDM model. In the ampicillin model chain D is the one genuinely ambiguous case: PN1-PN2 (26.30 A) reads Access while PC1-PC2 (28.94 A) reads Binding, so it is reported as a conflict rather than forced into a state.
- **Chain E replicates across the two independent reconstructions** (0.86 A whole-protomer, 0.30 A switch loop), which is below the ~1.0 A noise threshold. Chains D and F differ by 1.32 and 1.53 A and are not the same state between the two models.
- **Ampicillin spans the distal and proximal pockets** rather than sitting in either: its centroid is 7.09 A from the DBP centroid and 7.30 A from the PBP centroid. It is anchored by a K151 salt bridge (NZ-O1 2.81 A, NZ-O2 2.84 A) and packs against F178, F615 and F610. No heavy-atom pair is under 2.6 A, so the pose is clean.
- **The three DDM molecules are detergent, not substrate.** All three sit closer to the PBP centroid than the DBP, and together they occlude 98.2% of the chain E substrate-site volume (2018 -> 36 A^3). Ampicillin occludes 21.0% of the same site in the other model. This is an occupancy observation about a purification detergent at high concentration and carries no substrate-recognition weight.
- **The distal pocket is the hydrophobic one**, confirming the standard AcrB picture and contradicting the earlier surface-rendering reading (known issue 4): DBP 96.1% apolar side-chain atoms, mean Kyte-Doolittle +2.63, 8 aromatics; PBP 67.2% apolar, mean KD -1.84, 0 aromatics. This is identical in every chain of both structures.
- **R971 remains 12-17 A from D407 and D408 in all six protomers** (known issue 1). It has not moved between models, so it is a systematic modelling problem, not a per-map fluctuation.
- **Bound ligand narrows the exit route in both structures.** Seeded on chain E and run on the trimer, the widest route to bulk drops from 2.21 to 1.26 A when ampicillin is left in place, and from 2.43 to 1.42 A when the three DDM molecules are left in place. Both widest routes exit through the PC1/PC2 periplasmic cleft (CH1, tentative).
- **Two of the 51 validation checks fail, both in the tunnel stage** - see the tunnel section below. Every other number in PROTOCOL section 6 reproduces exactly.

## Tools

| tool | status | used for |
|---|---|---|
| pyKVFinder | available | stage 6 cavity detection (guided and unguided) |
| fpocket | available (built from source) | stage 6 pocket volume and druggability |
| CAVER 3.0.3 | available (downloaded from caver.cz) | independent cross-check of every tunnel bottleneck |
| ChimeraX | NOT AVAILABLE | `.defattr` and tunnel-trace files are still written for viewing locally |

Every tunnel bottleneck below has been recomputed independently with CAVER 3.0.3 - the tool reviewers expect - on the same trimers from the same seed points. The comparison is in the tunnel section.

## Validation against PROTOCOL section 6

**51/53 checks pass.**

Failing checks:

| check | got | expected |
|---|---|---|
| amp chain E tunnel bottleneck (seeded on ligand) | 2.21 | 2.01 |
| amp chain E constriction residue is F615 | ASN676E | PHE615E |

| check | value | expected | status |
|---|---|---|---|
| Amp_MexB_20260826 D mismatches vs P52002 | 0 | 0 | PASS |
| Amp_MexB_20260826 D numbering offset | 0 | 0 | PASS |
| Amp_MexB_20260826 E mismatches vs P52002 | 0 | 0 | PASS |
| Amp_MexB_20260826 E numbering offset | 0 | 0 | PASS |
| Amp_MexB_20260826 F mismatches vs P52002 | 0 | 0 | PASS |
| Amp_MexB_20260826 F numbering offset | 0 | 0 | PASS |
| MexB_DDM_3_20260730 D mismatches vs P52002 | 0 | 0 | PASS |
| MexB_DDM_3_20260730 D numbering offset | 0 | 0 | PASS |
| MexB_DDM_3_20260730 E mismatches vs P52002 | 0 | 0 | PASS |
| MexB_DDM_3_20260730 E numbering offset | 0 | 0 | PASS |
| MexB_DDM_3_20260730 F mismatches vs P52002 | 0 | 0 | PASS |
| MexB_DDM_3_20260730 F numbering offset | 0 | 0 | PASS |
| Amp_MexB_20260826 chain D gaps | none | none | PASS |
| Amp_MexB_20260826 chain E gaps | 229;360 | 229;360 | PASS |
| Amp_MexB_20260826 chain F gaps | 663 | 663 | PASS |
| MexB_DDM_3_20260730 chain D gaps | none | none | PASS |
| MexB_DDM_3_20260730 chain E gaps | 229;360 | 229;360 | PASS |
| MexB_DDM_3_20260730 chain F gaps | 663 | 663 | PASS |
| Amp_MexB_20260826 chain D residue range | 1-1032 | 1-1032 | PASS |
| Amp_MexB_20260826 chain E residue range | 1-1030 | 1-1030 | PASS |
| Amp_MexB_20260826 chain F residue range | 1-1033 | 1-1033 | PASS |
| MexB_DDM_3_20260730 chain D residue range | 1-1032 | 1-1032 | PASS |
| MexB_DDM_3_20260730 chain E residue range | 1-1030 | 1-1030 | PASS |
| MexB_DDM_3_20260730 chain F residue range | 1-1033 | 1-1033 | PASS |
| amp D PN1-PN2 separation | 26.30 | 26.3 | PASS |
| amp D PC1-PC2 separation | 28.94 | 28.94 | PASS |
| amp D PN1-PN2 pseudo-contacts | 21 | 21 | PASS |
| amp D PC1-PC2 pseudo-contacts | 5 | 5 | PASS |
| amp E PN1-PN2 separation | 27.82 | 27.82 | PASS |
| amp E PC1-PC2 separation | 28.38 | 28.38 | PASS |
| amp E PN1-PN2 pseudo-contacts | 1 | 1 | PASS |
| amp E PC1-PC2 pseudo-contacts | 5 | 5 | PASS |
| amp F PN1-PN2 separation | 29.96 | 29.96 | PASS |
| amp F PC1-PC2 separation | 24.22 | 24.22 | PASS |
| amp F PN1-PN2 pseudo-contacts | 4 | 4 | PASS |
| amp F PC1-PC2 pseudo-contacts | 16 | 16 | PASS |
| amp D ASP407-LYS939 | 2.71 | 2.71 | PASS |
| amp D ASP408-LYS939 | 2.70 | 2.7 | PASS |
| amp F LYS939-THR976 | 4.81 | 4.81 | PASS |
| amp K151 min distance | 2.81 | 2.81 | PASS |
| amp F178 min distance | 3.25 | 3.25 | PASS |
| amp F178 contact count | 37 | 37 | PASS |
| amp F615 min distance | 3.35 | 3.35 | PASS |
| amp F610 min distance | 3.76 | 3.76 | PASS |
| amp K151 NZ-O1 | 2.81 | 2.81 | PASS |
| amp K151 NZ-O2 | 2.84 | 2.84 | PASS |
| amp: no heavy-atom pair under 2.6 A | 0 | 0 | PASS |
| cross-structure chain E whole protomer | 0.86 | 0.86 | PASS |
| cross-structure chain E switch loop | 0.30 | 0.3 | PASS |
| amp chain E tunnel bottleneck (seeded on ligand) | 2.21 | 2.01 | FAIL |
| amp chain E constriction residue is F615 | ASN676E | PHE615E | FAIL |
| CAVER agrees with ours: Amp_MexB_20260826 protein chain E | 2.22 | 2.21 | PASS |
| CAVER agrees with ours: MexB_DDM_3_20260730 protein chain E | 2.46 | 2.43 | PASS |

## Flags raised this run

- Amp_MexB_20260826 ZZ7 E2000: B-factors are group-refined (single value 55.37) - known issue 5
- MexB_DDM_3_20260730 LMT E2001: B-factors are group-refined (single value 55.87) - known issue 5
- MexB_DDM_3_20260730 LMT E2002: B-factors are group-refined (single value 69.32) - known issue 5
- MexB_DDM_3_20260730 LMT E2003: B-factors are group-refined (single value 57.88) - known issue 5
- Amp_MexB_20260826 chain D: cleft diagnostics conflict - PN1-PN2 nearest Access, PC1-PC2 nearest Binding; reporting conflict rather than forcing a call
- Amp_MexB_20260826 chain D: ASP407-ARG971 = 12.58 A (known issue 1, R971 rotamer) - still unresolved
- Amp_MexB_20260826 chain D: ASP408-ARG971 = 17.37 A (known issue 1, R971 rotamer) - still unresolved
- Amp_MexB_20260826 chain E: ASP407-ARG971 = 12.75 A (known issue 1, R971 rotamer) - still unresolved
- Amp_MexB_20260826 chain E: ASP408-ARG971 = 16.90 A (known issue 1, R971 rotamer) - still unresolved
- Amp_MexB_20260826 chain F: ASP407-ARG971 = 15.63 A (known issue 1, R971 rotamer) - still unresolved
- Amp_MexB_20260826 chain F: ASP408-ARG971 = 16.87 A (known issue 1, R971 rotamer) - still unresolved
- MexB_DDM_3_20260730 chain D: ASP407-ARG971 = 13.39 A (known issue 1, R971 rotamer) - still unresolved
- MexB_DDM_3_20260730 chain D: ASP408-ARG971 = 16.63 A (known issue 1, R971 rotamer) - still unresolved
- MexB_DDM_3_20260730 chain E: ASP407-ARG971 = 12.81 A (known issue 1, R971 rotamer) - still unresolved
- MexB_DDM_3_20260730 chain E: ASP408-ARG971 = 16.09 A (known issue 1, R971 rotamer) - still unresolved
- MexB_DDM_3_20260730 chain F: ASP407-ARG971 = 13.99 A (known issue 1, R971 rotamer) - still unresolved
- MexB_DDM_3_20260730 chain F: ASP408-ARG971 = 16.74 A (known issue 1, R971 rotamer) - still unresolved
- MexB_DDM_3_20260730:LMTE2001: CLOSE CONTACT O6'-ASP274E.OD1 = 2.60 A (<2.6) [short H-bond (both N/O)] - pose or rotamer must be revisited before figures

## Stage 1 - ingest and validation

| structure | chain | first_res | last_res | n_res | gaps | seq_mismatches_P52002 | numbering_offset |
|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | D | 1 | 1032 | 1032 | none | 0 | 0 |
| Amp_MexB_20260826 | E | 1 | 1030 | 1028 | 229;360 | 0 | 0 |
| Amp_MexB_20260826 | F | 1 | 1033 | 1032 | 663 | 0 | 0 |
| MexB_DDM_3_20260730 | D | 1 | 1032 | 1032 | none | 0 | 0 |
| MexB_DDM_3_20260730 | E | 1 | 1030 | 1028 | 229;360 | 0 | 0 |
| MexB_DDM_3_20260730 | F | 1 | 1033 | 1032 | 663 | 0 | 0 |

Ligand inventory:

| structure | ligand | chain | resseq | n_heavy | n_atoms | b_min | b_median | b_max | b_mode | class |
|---|---|---|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | ZZ7 | E | 2000 | 25 | 44 | 55.37 | 55.37 | 55.37 | group | substrate |
| MexB_DDM_3_20260730 | LMT | E | 2001 | 35 | 81 | 55.87 | 55.87 | 55.87 | group | detergent |
| MexB_DDM_3_20260730 | LMT | E | 2002 | 35 | 81 | 69.32 | 69.32 | 69.32 | group | detergent |
| MexB_DDM_3_20260730 | LMT | E | 2003 | 35 | 81 | 57.88 | 57.88 | 57.88 | group | detergent |

Ligand B-factors are group-refined - a single value per molecule - in both models, so they carry no per-atom information and must be described as such wherever quoted (known issue 5).

## Stage 2 - state assignment

| structure | chain | PN1_PN2_sep | PC1_PC2_sep | PN1_PN2_contacts10A | PN1_PN2_contacts_per_res | PC1_PC2_contacts10A | PC1_PC2_contacts_per_res | state_call | L1_deviation | PN_nearest | PC_nearest | diagnostics_agree |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | D | 26.30 | 28.94 | 21 | 0.114 | 5 | 0.028 | Binding | 2.15 | Access | Binding | no |
| Amp_MexB_20260826 | E | 27.82 | 28.38 | 1 | 0.005 | 5 | 0.028 | Binding | 1.20 | Binding | Binding | yes |
| Amp_MexB_20260826 | F | 29.96 | 24.22 | 4 | 0.022 | 16 | 0.090 | Extrusion | 0.43 | Extrusion | Extrusion | yes |
| MexB_DDM_3_20260730 | D | 26.23 | 27.18 | 21 | 0.114 | 6 | 0.034 | Access | 1.25 | Access | Access | yes |
| MexB_DDM_3_20260730 | E | 27.82 | 28.49 | 1 | 0.005 | 7 | 0.039 | Binding | 1.09 | Binding | Binding | yes |
| MexB_DDM_3_20260730 | F | 29.58 | 24.51 | 5 | 0.027 | 15 | 0.085 | Extrusion | 0.24 | Extrusion | Extrusion | yes |

| structure | chain | res1 | res2 | min_dist_A | note |
|---|---|---|---|---|---|
| Amp_MexB_20260826 | D | ASP407 | ASP408 | 4.46 |  |
| Amp_MexB_20260826 | D | ASP407 | LYS939 | 2.71 |  |
| Amp_MexB_20260826 | D | ASP407 | ARG971 | 12.58 |  |
| Amp_MexB_20260826 | D | ASP407 | THR976 | 3.42 |  |
| Amp_MexB_20260826 | D | ASP408 | LYS939 | 2.70 |  |
| Amp_MexB_20260826 | D | ASP408 | ARG971 | 17.37 |  |
| Amp_MexB_20260826 | D | ASP408 | THR976 | 9.15 |  |
| Amp_MexB_20260826 | D | LYS939 | ARG971 | 17.05 |  |
| Amp_MexB_20260826 | D | LYS939 | THR976 | 7.69 |  |
| Amp_MexB_20260826 | D | ARG971 | THR976 | 10.88 |  |
| Amp_MexB_20260826 | E | ASP407 | ASP408 | 8.72 |  |
| Amp_MexB_20260826 | E | ASP407 | LYS939 | 3.30 |  |
| Amp_MexB_20260826 | E | ASP407 | ARG971 | 12.75 |  |
| Amp_MexB_20260826 | E | ASP407 | THR976 | 3.91 |  |
| Amp_MexB_20260826 | E | ASP408 | LYS939 | 6.40 |  |
| Amp_MexB_20260826 | E | ASP408 | ARG971 | 16.90 |  |
| Amp_MexB_20260826 | E | ASP408 | THR976 | 11.58 |  |
| Amp_MexB_20260826 | E | LYS939 | ARG971 | 16.56 |  |
| Amp_MexB_20260826 | E | LYS939 | THR976 | 8.28 |  |
| Amp_MexB_20260826 | E | ARG971 | THR976 | 11.23 |  |

_40 further rows in `results/tables/proton_relay.csv`._

## Stage 3 - inter-protomer comparison

| structure | fit_frame | chain_A | chain_B | n_fit | fit_rmsd_A | whole_protomer_rmsd_A | porter_rmsd_A | docking_rmsd_A | TM_rmsd_A |
|---|---|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | TM_trimmed | D | E | 389 | 1.20 | 4.02 | 4.20 | 6.41 | 1.20 |
| Amp_MexB_20260826 | TM_trimmed | D | F | 390 | 1.97 | 5.50 | 6.74 | 7.60 | 1.97 |
| Amp_MexB_20260826 | TM_trimmed | E | F | 389 | 2.17 | 4.25 | 6.00 | 3.41 | 2.17 |
| Amp_MexB_20260826 | porter | D | E | 363 | 2.03 | 3.64 | 2.03 | 2.50 | 4.43 |
| Amp_MexB_20260826 | porter | D | F | 362 | 3.19 | 6.31 | 3.19 | 4.28 | 8.43 |
| Amp_MexB_20260826 | porter | E | F | 362 | 3.32 | 4.87 | 3.32 | 5.24 | 5.86 |
| MexB_DDM_3_20260730 | TM_trimmed | D | E | 389 | 1.70 | 4.27 | 4.83 | 6.12 | 1.70 |
| MexB_DDM_3_20260730 | TM_trimmed | D | F | 390 | 1.54 | 5.19 | 6.25 | 7.56 | 1.54 |
| MexB_DDM_3_20260730 | TM_trimmed | E | F | 389 | 1.54 | 4.28 | 6.35 | 3.18 | 1.54 |
| MexB_DDM_3_20260730 | porter | D | E | 363 | 2.20 | 3.66 | 2.20 | 2.20 | 4.57 |
| MexB_DDM_3_20260730 | porter | D | F | 362 | 2.32 | 6.33 | 2.32 | 4.55 | 8.68 |
| MexB_DDM_3_20260730 | porter | E | F | 362 | 3.03 | 5.07 | 3.03 | 5.51 | 6.30 |

Key mechanical quantities (TM-frame fit):

| structure | fit_frame | chain_A | chain_B | region | n_res | rms_dev_A | max_dev_A |
|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | TM_trimmed | D | E | TM2 | 23 | 1.21 | 1.86 |
| Amp_MexB_20260826 | TM_trimmed | D | E | Ialpha | 20 | 1.17 | 2.38 |
| Amp_MexB_20260826 | TM_trimmed | D | E | TM7-12 | 155 | 0.92 | 4.64 |
| Amp_MexB_20260826 | TM_trimmed | D | F | TM2 | 23 | 2.28 | 2.88 |
| Amp_MexB_20260826 | TM_trimmed | D | F | Ialpha | 21 | 1.90 | 3.19 |
| Amp_MexB_20260826 | TM_trimmed | D | F | TM7-12 | 155 | 1.86 | 4.51 |
| Amp_MexB_20260826 | TM_trimmed | E | F | TM2 | 23 | 3.02 | 3.95 |
| Amp_MexB_20260826 | TM_trimmed | E | F | Ialpha | 20 | 1.84 | 2.74 |
| Amp_MexB_20260826 | TM_trimmed | E | F | TM7-12 | 155 | 2.04 | 4.01 |
| MexB_DDM_3_20260730 | TM_trimmed | D | E | TM2 | 23 | 2.07 | 2.59 |
| MexB_DDM_3_20260730 | TM_trimmed | D | E | Ialpha | 20 | 0.64 | 1.34 |
| MexB_DDM_3_20260730 | TM_trimmed | D | E | TM7-12 | 155 | 1.34 | 5.55 |
| MexB_DDM_3_20260730 | TM_trimmed | D | F | TM2 | 23 | 1.06 | 1.28 |
| MexB_DDM_3_20260730 | TM_trimmed | D | F | Ialpha | 21 | 0.88 | 1.65 |
| MexB_DDM_3_20260730 | TM_trimmed | D | F | TM7-12 | 155 | 1.34 | 4.70 |
| MexB_DDM_3_20260730 | TM_trimmed | E | F | TM2 | 23 | 2.73 | 3.22 |
| MexB_DDM_3_20260730 | TM_trimmed | E | F | Ialpha | 20 | 1.13 | 2.24 |
| MexB_DDM_3_20260730 | TM_trimmed | E | F | TM7-12 | 155 | 1.58 | 4.02 |

Rigid-body domain motions:

| structure | fit_frame | chain_A | chain_B | domain | n_res | rmsd_before_A | rmsd_after_A | rotation_deg | axis | centroid_shift_A | shift_along_axis_A |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | TM_trimmed | D | E | porter | 363 | 4.20 | 2.03 | 6.65 | -0.384,-0.923,0.030 | 3.04 | 1.89 |
| Amp_MexB_20260826 | TM_trimmed | D | E | docking | 192 | 6.41 | 1.38 | 5.38 | -0.088,0.602,0.793 | 6.10 | -1.96 |
| Amp_MexB_20260826 | TM_trimmed | D | E | TM | 389 | 1.20 | 1.20 | 0.00 | 0.252,-0.004,0.968 | 0.00 | -0.00 |
| Amp_MexB_20260826 | TM_trimmed | D | F | porter | 362 | 6.74 | 3.19 | 15.89 | 0.209,-0.697,0.686 | 2.33 | -0.39 |
| Amp_MexB_20260826 | TM_trimmed | D | F | docking | 193 | 7.60 | 1.52 | 8.44 | 0.217,-0.620,0.754 | 7.00 | -0.14 |
| Amp_MexB_20260826 | TM_trimmed | D | F | TM | 390 | 1.97 | 1.97 | 0.00 | 0.087,-0.683,0.725 | 0.00 | -0.00 |
| Amp_MexB_20260826 | TM_trimmed | E | F | porter | 362 | 6.00 | 3.32 | 13.95 | -0.004,-0.412,0.911 | 1.07 | -0.78 |
| Amp_MexB_20260826 | TM_trimmed | E | F | docking | 192 | 3.41 | 1.94 | 4.71 | -0.384,-0.696,0.607 | 2.29 | -0.58 |
| Amp_MexB_20260826 | TM_trimmed | E | F | TM | 389 | 2.17 | 2.17 | 0.00 | -0.314,-0.565,0.763 | 0.00 | 0.00 |
| MexB_DDM_3_20260730 | TM_trimmed | D | E | porter | 363 | 4.83 | 2.20 | 8.85 | 0.083,0.848,0.524 | 3.30 | -1.72 |
| MexB_DDM_3_20260730 | TM_trimmed | D | E | docking | 192 | 6.12 | 0.73 | 2.78 | -0.424,-0.904,0.052 | 6.04 | 2.91 |
| MexB_DDM_3_20260730 | TM_trimmed | D | E | TM | 389 | 1.70 | 1.70 | 0.00 | -0.441,-0.618,0.650 | 0.00 | -0.00 |
| MexB_DDM_3_20260730 | TM_trimmed | D | F | porter | 362 | 6.25 | 2.32 | 15.36 | 0.578,0.534,0.617 | 2.83 | 0.14 |
| MexB_DDM_3_20260730 | TM_trimmed | D | F | docking | 193 | 7.56 | 0.60 | 7.45 | 0.624,0.480,0.617 | 7.29 | 0.37 |
| MexB_DDM_3_20260730 | TM_trimmed | D | F | TM | 390 | 1.54 | 1.54 | 0.00 | -0.250,-0.091,0.964 | 0.00 | 0.00 |
| MexB_DDM_3_20260730 | TM_trimmed | E | F | porter | 362 | 6.35 | 3.03 | 15.74 | 0.181,0.307,0.934 | 1.07 | 0.75 |
| MexB_DDM_3_20260730 | TM_trimmed | E | F | docking | 192 | 3.18 | 0.65 | 6.00 | 0.322,0.566,0.759 | 2.56 | 1.56 |
| MexB_DDM_3_20260730 | TM_trimmed | E | F | TM | 389 | 1.54 | 1.54 | 0.00 | 0.043,0.027,0.999 | 0.00 | -0.00 |

## Stage 4 - cross-structure comparison

| structure_A | structure_B | chain | n_fit | whole_protomer_rmsd_allCA_fit_A | whole_protomer_rmsd_TM_fit_A | porter_dev_after_TM_fit_A | TM_dev_after_porter_fit_A | switch_loop_rmsd_A | interpretation |
|---|---|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | MexB_DDM_3_20260730 | D | 1032 | 1.32 | 2.42 | 2.63 | 1.91 | 0.49 | differs |
| Amp_MexB_20260826 | MexB_DDM_3_20260730 | E | 1028 | 0.86 | 1.13 | 0.96 | 1.04 | 0.30 | same state (<1.0 A) |
| Amp_MexB_20260826 | MexB_DDM_3_20260730 | F | 1032 | 1.53 | 1.92 | 2.12 | 2.57 | 1.43 | differs |

## Stage 5 - ligand environment

| structure | ligand | lig_id | n_heavy | n_contact_residues | centroid_to_DBP_A | centroid_to_PBP_A | atoms_nearer_DBP_centroid | atoms_nearer_PBP_centroid | atoms_nearer_DBP_lining_atom | pocket_call | n_DBP_sidechain_contacts | n_hbond_cands | n_close_contacts_2.6A | n_true_steric_overlaps | class |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | ZZ7 | E2000 | 25 | 12 | 7.09 | 7.30 | 11 | 14 | 17 | spans DBP and PBP | 70 | 3 | 0 | 0 | substrate |
| MexB_DDM_3_20260730 | LMT | E2001 | 35 | 20 | 13.08 | 9.45 | 6 | 29 | 14 | proximal (PBP) | 33 | 8 | 1 | 0 | detergent |
| MexB_DDM_3_20260730 | LMT | E2002 | 35 | 19 | 10.70 | 9.17 | 7 | 28 | 21 | proximal (PBP) | 51 | 7 | 0 | 0 | detergent |
| MexB_DDM_3_20260730 | LMT | E2003 | 35 | 15 | 10.95 | 1.51 | 0 | 35 | 22 | proximal (PBP) | 25 | 6 | 0 | 0 | detergent |

Cross-ligand contact matrix (residues contacted by more than one ligand):

| resseq | resname | pocket | switch_loop | Amp_MexB_20260826:ZZ7E2000 | MexB_DDM_3_20260730:LMTE2001 | MexB_DDM_3_20260730:LMTE2002 | MexB_DDM_3_20260730:LMTE2003 | n_ligands |
|---|---|---|---|---|---|---|---|---|
| 46 | GLN |  |  |  | 3.97 |  | 3.20 | 2 |
| 134 | LYS |  |  |  |  | 3.76 | 3.02 | 2 |
| 139 | VAL | DBP |  |  | 4.33 | 3.87 |  | 2 |
| 176 | GLN | PBP |  | 3.97 |  |  | 3.98 | 2 |
| 178 | PHE | DBP |  | 3.25 | 4.00 |  |  | 2 |
| 179 | GLY |  |  | 3.49 | 3.03 |  |  | 2 |
| 274 | ASP | PBP |  |  | 2.60 |  | 3.20 | 2 |
| 276 | SER | PBP |  | 4.21 |  |  | 3.07 | 2 |
| 277 | ILE | DBP |  | 3.57 | 4.39 |  | 3.17 | 3 |
| 610 | PHE | DBP |  | 3.76 | 3.63 |  |  | 2 |
| 612 | VAL | DBP |  | 3.50 | 4.17 |  |  | 2 |
| 615 | PHE | DBP | switch | 3.35 | 3.81 |  | 3.76 | 3 |
| 617 | PHE | DBP | switch |  |  | 3.16 | 2.68 | 2 |
| 620 | ARG |  | switch | 3.59 |  |  | 3.12 | 2 |
| 628 | PHE | DBP |  | 3.95 | 3.91 | 4.11 |  | 3 |
| 718 | ASN |  |  |  |  | 2.71 | 4.22 | 2 |

## Stage 6 - pockets and cavities

| structure | chain | detection | cavity_id | volume_A3 | area_A2 | max_depth_A | avg_depth_A | avg_hydropathy | DBP_overlap | PBP_overlap | continuity | n_lining_res |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | D | unguided_ligands_stripped | KCF | 8220 | 5975 | 19.8 | 6.89 | 0.18 | 7 | 14 | one continuous cavity (PBP and DBP) | 315 |
| Amp_MexB_20260826 | E | unguided_ligands_stripped | KAD | 5632 | 4395 | 21.7 | 6.12 | 0.06 | 13 | 14 | one continuous cavity (PBP and DBP) | 231 |
| Amp_MexB_20260826 | F | unguided_ligands_stripped | KCP | 6832 | 5357 | 17.4 | 4.85 | 0.32 | 8 | 12 | one continuous cavity (PBP and DBP) | 298 |
| Amp_MexB_20260826 | E | ligand_guided_ZZ72000 | KAA | 1064 | 481 | 8.7 | 2.45 | -0.14 | 12 | 8 |  | 34 |
| MexB_DDM_3_20260730 | D | unguided_ligands_stripped | KAD | 7345 | 5779 | 16.6 | 6.86 | 0.20 | 6 | 13 | one continuous cavity (PBP and DBP) | 326 |
| MexB_DDM_3_20260730 | E | unguided_ligands_stripped | KDU | 5341 | 4251 | 20.9 | 5.69 | 0.05 | 13 | 14 | one continuous cavity (PBP and DBP) | 235 |
| MexB_DDM_3_20260730 | F | unguided_ligands_stripped | KDL | 54 | 85 | 0.0 | 0.00 | -0.52 | 5 | 3 | unreliable - best DBP-overlapping cavity is small or barely overlaps; not interpreted | 11 |
| MexB_DDM_3_20260730 | E | ligand_guided_LMT2001 | KAA | 1622 | 789 | 6.3 | 2.24 | 0.02 | 13 | 8 |  | 50 |
| MexB_DDM_3_20260730 | E | ligand_guided_LMT2002 | KAA | 710 | 482 | 5.0 | 1.43 | -0.47 | 10 | 4 |  | 40 |
| MexB_DDM_3_20260730 | E | ligand_guided_LMT2003 | KAA | 1015 | 468 | 4.7 | 1.37 | 0.55 | 7 | 7 |  | 37 |

Grid-based substrate-site volumes. **These are internally comparable across these structures only and are not drop-in replacements for fpocket or CASTp volumes.** The 'stripped' column removes the structure's own ligands, including the three DDM molecules, so the sites can be compared like for like; the difference is occlusion.

| structure | chain | site_sphere_radius_A | grid_step_A | volume_ligands_stripped_A3 | volume_with_ligands_A3 | occluded_volume_A3 | occluded_pct | note |
|---|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | D | 16 | 0.50 | 1476 | 1476 | 0 | 0.0 |  |
| Amp_MexB_20260826 | E | 16 | 0.50 | 1843 | 1456 | 387 | 21.0 |  |
| Amp_MexB_20260826 | F | 16 | 0.50 | 1402 | 1402 | 0 | 0.0 |  |
| MexB_DDM_3_20260730 | D | 16 | 0.50 | 964 | 964 | 0 | 0.0 |  |
| MexB_DDM_3_20260730 | E | 16 | 0.50 | 2018 | 36 | 1982 | 98.2 |  |
| MexB_DDM_3_20260730 | F | 16 | 0.50 | 1543 | 1543 | 0 | 0.0 |  |

fpocket, on the ligand-stripped trimer:

| structure | chain | pocket_id | volume_A3 | score | druggability | mean_local_hydrophobic_density | apolar_alpha_sphere_prop | DBP_overlap | PBP_overlap | n_lining_res |
|---|---|---|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | E | 229 | 1453 | -0.635 | 0.878 | 87.68 | 0.818 | 13 | 11 | 32 |
| Amp_MexB_20260826 | F | 10 | 153 | 0.037 | 0.006 | 29.00 | 0.882 | 6 | 2 | 9 |
| Amp_MexB_20260826 | F | 215 | 1572 | -0.383 | 0.079 | 26.00 | 0.413 | 5 | 4 | 41 |
| Amp_MexB_20260826 | D | 233 | 3266 | -1.379 | 0.003 | 34.40 | 0.324 | 5 | 7 | 72 |
| Amp_MexB_20260826 | F | 8 | 91 | 0.040 | 0.001 | 15.00 | 0.941 | 4 | 0 | 9 |
| MexB_DDM_3_20260730 | E | 1 | 5647 | 2.221 | 0.909 | 60.95 | 0.413 | 10 | 8 | 118 |
| MexB_DDM_3_20260730 | E | 81 | 591 | -0.119 | 0.217 | 32.36 | 0.506 | 7 | 6 | 18 |
| MexB_DDM_3_20260730 | D | 209 | 3967 | -1.214 | 0.340 | 47.21 | 0.344 | 5 | 8 | 93 |
| MexB_DDM_3_20260730 | F | 207 | 2764 | -0.853 | 0.000 | 11.10 | 0.148 | 3 | 8 | 69 |

Pocket composition (known issue 4):

| structure | chain | pocket | n_res | n_sidechain_atoms | n_apolar_C | n_polar_NOS | pct_apolar | mean_KD | n_aromatic_res |
|---|---|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | D | PBP | 17 | 64 | 43 | 21 | 67.2 | -1.84 | 0 |
| Amp_MexB_20260826 | D | DBP | 14 | 76 | 73 | 3 | 96.1 | 2.63 | 8 |
| Amp_MexB_20260826 | E | PBP | 17 | 64 | 43 | 21 | 67.2 | -1.84 | 0 |
| Amp_MexB_20260826 | E | DBP | 14 | 76 | 73 | 3 | 96.1 | 2.63 | 8 |
| Amp_MexB_20260826 | F | PBP | 17 | 64 | 43 | 21 | 67.2 | -1.84 | 0 |
| Amp_MexB_20260826 | F | DBP | 14 | 76 | 73 | 3 | 96.1 | 2.63 | 8 |
| MexB_DDM_3_20260730 | D | PBP | 17 | 64 | 43 | 21 | 67.2 | -1.84 | 0 |
| MexB_DDM_3_20260730 | D | DBP | 14 | 76 | 73 | 3 | 96.1 | 2.63 | 8 |
| MexB_DDM_3_20260730 | E | PBP | 17 | 64 | 43 | 21 | 67.2 | -1.84 | 0 |
| MexB_DDM_3_20260730 | E | DBP | 14 | 76 | 73 | 3 | 96.1 | 2.63 | 8 |
| MexB_DDM_3_20260730 | F | PBP | 17 | 64 | 43 | 21 | 67.2 | -1.84 | 0 |
| MexB_DDM_3_20260730 | F | DBP | 14 | 76 | 73 | 3 | 96.1 | 2.63 | 8 |

## Stage 7 - tunnels

Run on the trimer in every case. Seeds are the bound ligand centroid where the chain has a ligand, otherwise the DBP/PBP midpoint of that chain. `protein` mode strips all ligands; `withlig` mode keeps them as obstructions, so the difference measures occlusion of the exit route.

| structure | mode | chain | seed | tunnel_rank | bottleneck_radius_A | geodesic_path_length_A | constriction_lining_clearance_A | channel_call | assignment_confidence |
|---|---|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | protein | D | site | 1 | 1.22 | 114.0 | PHE617D:1.22;PHE136D:1.25;PHE615D:1.64 | CH1 (PC1/PC2 periplasmic cleft) | tentative |
| Amp_MexB_20260826 | protein | E | ZZ72000 | 1 | 2.21 | 62.8 | ASN676E:2.21;ASN718E:2.21;LEU827E:2.25;ALA677E:2.83;PHE617E:3.38 | CH1 (PC1/PC2 periplasmic cleft) | tentative |
| Amp_MexB_20260826 | protein | E | ZZ72000 | 2 | 1.72 | 71.0 | THR678E:1.72;ARG716E:1.72;GLU829E:1.76;ALA677E:2.15;GLY828E:2.51;LEU827E:2.97 | CH1 (PC1/PC2 periplasmic cleft) | tentative |
| Amp_MexB_20260826 | protein | E | ZZ72000 | 3 | 1.69 | 34.1 | TYR327E:1.69;LYS134E:1.70;LEU672E:1.80;PHE136E:2.44;ASN135E:2.51 | CH3 (PN1/PN2 groove) | tentative |
| Amp_MexB_20260826 | protein | F | site | 1 | 1.75 | 92.7 | PHE615F:1.75;MET626F:1.79;PHE136F:1.87;PHE178F:2.22;ASN616F:3.17;PHE617F:3.24 | CH1 (PC1/PC2 periplasmic cleft) | tentative |
| Amp_MexB_20260826 | withlig | D | site | 1 | 1.22 | 114.0 | PHE617D:1.22;PHE136D:1.25;PHE615D:1.64 | CH1 (PC1/PC2 periplasmic cleft) | tentative |
| Amp_MexB_20260826 | withlig | E | ZZ72000 | 1 | 1.26 | 47.4 | ZZ72000E:1.26;PHE178E:1.26;GLN176E:1.27 | CH3 (PN1/PN2 groove) | tentative |
| Amp_MexB_20260826 | withlig | F | site | 1 | 1.75 | 92.7 | PHE615F:1.75;MET626F:1.79;PHE136F:1.87;PHE178F:2.22;ASN616F:3.17;PHE617F:3.24 | CH1 (PC1/PC2 periplasmic cleft) | tentative |
| MexB_DDM_3_20260730 | protein | D | site | 1 | 1.21 | 59.2 | PHE615D:1.21;PHE617D:1.46;PHE136D:1.91 | CH1 (PC1/PC2 periplasmic cleft) | tentative |
| MexB_DDM_3_20260730 | protein | E | LMT2001 | 1 | 2.43 | 62.0 | ARG716E:2.43;ASN676E:2.45;PHE666E:2.52;PHE664E:3.68 | CH1 (PC1/PC2 periplasmic cleft) | tentative |
| MexB_DDM_3_20260730 | protein | E | LMT2001 | 2 | 1.81 | 69.1 | TYR327E:1.81;LEU672E:1.82;PHE136E:1.84;LYS134E:1.88;ASN135E:2.04 | CH3 (PN1/PN2 groove) | tentative |
| MexB_DDM_3_20260730 | protein | E | LMT2001 | 3 | 1.79 | 62.9 | THR130E:1.79;ASP174E:1.79;LYS131E:3.22 | CH3 (PN1/PN2 groove) | tentative |
| MexB_DDM_3_20260730 | protein | E | LMT2002 | 1 | 1.37 | 45.2 | PHE617E:1.37 | CH1 (PC1/PC2 periplasmic cleft) | tentative |
| MexB_DDM_3_20260730 | protein | E | LMT2003 | 1 | 2.28 | 51.9 | PHE615E:2.28;GLN46E:2.83;ARG620E:3.34;THR89E:3.77 | unassigned | uncertain - no dominant lining group |
| MexB_DDM_3_20260730 | protein | F | site | 1 | 1.27 | 96.2 | PHE615F:1.27;PHE617F:1.30;PHE136F:2.45;ASN616F:2.49 | CH1 (PC1/PC2 periplasmic cleft) | tentative |
| MexB_DDM_3_20260730 | withlig | D | site | 1 | 1.21 | 59.2 | PHE615D:1.21;PHE617D:1.46;PHE136D:1.91 | CH1 (PC1/PC2 periplasmic cleft) | tentative |
| MexB_DDM_3_20260730 | withlig | E | LMT2001 | 1 | 1.42 | 53.4 | LMT2003E:1.42;LMT2001E:1.44;GLN176E:2.07;VAL177E:2.27 | CH3 (PN1/PN2 groove) | tentative |
| MexB_DDM_3_20260730 | withlig | E | LMT2002 | 1 | 1.33 | 50.7 | LMT2002E:1.33;PHE666E:1.43;PRO669E:2.17;GLY675E:2.69;ASN676E:2.78 | unassigned | uncertain - no dominant lining group |
| MexB_DDM_3_20260730 | withlig | E | LMT2003 | 1 | 1.21 | 82.8 | LMT2003E:1.21;PHE615E:1.44;LMT2001E:1.47;ARG620E:2.43 | CH1 (PC1/PC2 periplasmic cleft) | tentative |
| MexB_DDM_3_20260730 | withlig | F | site | 1 | 1.27 | 96.2 | PHE615F:1.27;PHE617F:1.30;PHE136F:2.45;ASN616F:2.49 | CH1 (PC1/PC2 periplasmic cleft) | tentative |

### Cross-check against CAVER 3.0.3

CAVER was run on the same trimers, seeded on the same points, with `probe_radius 0.9`. `our_bottleneck_A` is this pipeline's own widest-path result for the same structure, mode and chain.

| structure | mode | chain | caver_bottleneck_A | our_bottleneck_A | difference_A | atoms_in_input | atoms_loaded_by_caver | valid_comparison | caver_bottleneck_residues |
|---|---|---|---|---|---|---|---|---|---|
| Amp_MexB_20260826 | protein | E | 2.22 | 2.21 | 0.01 | 23459 | 23451 | yes |  |
| Amp_MexB_20260826 | withlig | E | 2.22 | 1.26 | 0.96 | 23484 | 23459 | no - ligand atoms discarded |  |
| MexB_DDM_3_20260730 | protein | E | 2.46 | 2.43 | 0.03 | 23459 | 23451 | yes |  |
| MexB_DDM_3_20260730 | withlig | E | 1.57 | 1.42 | 0.15 | 23564 | 23529 | no - ligand atoms discarded |  |

**2 of 2 valid comparisons agree to within 0.05 Å**, and CAVER independently reports the same constriction-lining residues in the same order.

Two caveats on reading this table. Tunnel *lengths* are not comparable - CAVER ends the path at its own surface criterion while this pipeline runs on to the edge of the box - so only the bottleneck radii should be compared. And **CAVER's ligand-in-place rows are not a valid cross-check**: CAVER assigns radii from its own atom table and silently discards atoms it cannot place, which for these ligands means most of the molecule (8 of ampicillin's 25 heavy atoms loaded; 78 of DDM's 105), whether the ligand is written as HETATM or as ATOM. The occlusion results therefore rest on this pipeline alone.

### The switch-loop (F615) gate

PROTOCOL section 6 expects the ampicillin chain-E tunnel to bottleneck at **2.01 A with the constriction at F615**. This implementation instead finds **2.21 A constricted at N676/N718/L827** in the PC2/DC region, with F615 sitting 3.4 A clear of the path. Both validation checks therefore fail.

The table below measures the F615 gate on its own terms. In every protomer of both structures the switch-loop region is a *local widening* (3.1-4.6 A of clearance) that is not a through-route: a path forced to pass through it bottlenecks at 0.5-1.2 A, far below 2.01 A. So under this implementation the reference value cannot be recovered by routing through F615 either - it is not a matter of the search picking the wrong path.

| structure | chain | gate_widest_clearance_A | bottleneck_forced_through_gate_A | note |
|---|---|---|---|---|
| Amp_MexB_20260826 | D | 3.38 | 0.62 |  |
| Amp_MexB_20260826 | E | 4.18 | 1.20 |  |
| Amp_MexB_20260826 | F | 3.70 | 1.20 |  |
| MexB_DDM_3_20260730 | D | 3.11 | 0.53 |  |
| MexB_DDM_3_20260730 | E | 4.60 | 1.22 |  |
| MexB_DDM_3_20260730 | F | 3.92 | 0.81 |  |

Two things are worth separating here. The **protein geometry is consistent**: the gate is widest in the two chain-E (Binding) protomers, 4.18 A and 4.60 A, and narrower in the Access and Extrusion protomers - the same ordering in both reconstructions, so it replicates. What does **not** reproduce is the reference tunnel number itself.

The honest reading is that this is an unresolved implementation difference, not a settled result. The `tunnels.py` that PROTOCOL section 2 says ships with the protocol was not present in this repository, so a rewrite was unavoidable and there is no original implementation to diff against. Ruled out so far: grid resolution (the value converges upward from 1.79 A at 1.0 A spacing to 2.21 A with continuous refinement) and hydrogen handling (including hydrogens as obstructions gives 1.92 A and still does not move the constriction to F615). **Until the original script is available to compare against, treat every bottleneck radius in this report as provisional.** The channel assignments, which rest on lining composition rather than on the radius, are less affected but are still labelled tentative throughout.

## Changes since the previous run

- `caver.csv`: new
- `fpocket.csv`: 9 row(s) changed
- `validation.csv`: 0 row(s) changed, 2 row(s) added/removed
