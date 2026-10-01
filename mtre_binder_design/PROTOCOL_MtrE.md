# MtrE extracellular binder design — protocol

Design campaign for mini-protein and macrocycle binders against **MtrE**, the
outer-membrane channel of the *Neisseria gonorrhoeae* MtrCDE efflux pump.

**Scope constraint: extracellular only.** Nothing in this plan requires delivery
across the outer membrane. Every site defined here is reachable by a molecule
applied to intact cells.

Method template: Clement *et al.* 2026, *A complete RXFP1–relaxin interaction
model unlocks the design of potent mini-protein modulators*
(bioRxiv 2026.06.19.733483).

---

## 0. Ground truth — do not re-derive these

Everything below was measured from the coordinates by
`scripts/prepare_mtre_targets.py`, which re-checks all of it on every run and
exits rather than writing targets if any check fails.

**Structure.** PDB **4MT0**, biological assembly 1 — MtrE open state, 3.29 Å,
trimeric, 445 modelled residues per protomer.

**Numbering — this is the trap.** 4MT0 author numbering runs **21–465** and is
**precursor numbering, including the 20-residue signal peptide**.

```
auth = mature + 20
```

The PLOS paper describing this structure quotes the aspartate ring as
**D402/D405**; the deposited coordinates call the same two residues
**D422/D425**. Both are correct in their own frame. Every residue identifier in
this directory, in every config file, is **auth numbering**, matching the PDB.
A design run launched against the wrong frame targets a site 20 residues away
and will produce confident, worthless designs.

**Frame.** The trimer is centred on its centroid with the crystallographic
three-fold on *z*, oriented so the β-barrel is at **+z** (extracellular) and the
periplasmic tip at **−z**. Total span 128.5 Å.

| feature | auth residues | z (Å) |
| --- | --- | --- |
| extracellular apex | 116–118, 328–329 | +66 to +68 |
| outer-leaflet aromatic girdle | F326 | +59.5 |
| β-barrel strands | 106–112, 119–129, 308–319, 331–340 | +40 to +65 |
| periplasmic-leaflet girdle | F137, F305, F306, F344, W346 | +34 to +38 |
| equatorial domain | — | 0 to +20 |
| **aspartate ring D422/D425** | 422, 425 | **−46 to −51** |
| periplasmic tip | 416 | −61 |

**Membrane slab: z = +35 to +60.** Set by the two aromatic girdles, 23.7 Å
hydrophobic core. Anything a binder does below z = +58 is happening in lipid.

**Extracellular loops.**

| loop | auth | mature | sequence |
| --- | --- | --- | --- |
| Loop 1 | 113–118 | 93–98 | SLSGGN |
| **Loop 2** | **319–331** | **299–311** | **SVELGGLFKSGTG** |

Loop 2 is the 13-residue epitope conserved in >98% of gonococcal isolates,
recognised by bactericidal antibodies on live cells, and protective as a MAP
vaccine and as a mAb. **Surface accessibility of this loop on intact gonococci
is experimentally established** — that is the single biggest de-risking fact
available for this project.

**Pore profile (measured, not assumed).**

```
z=+66   r=9.3   extracellular mouth
z=+62   r=7.4
z=+50   r=8.9   mid-barrel
z=+34   r=9.4   deepest point reachable from outside
...
z=-46   r=4.3   D425 — THE GATE, periplasmic, unreachable
```

**The barrel lumen is an open funnel.** Pore radius never drops below 7.4 Å
anywhere between the extracellular mouth and z = +34. There is **no
extracellular constriction**: a binder can insert ~32 Å into the channel
without passing a gate. The only constriction in the whole 128 Å channel is the
aspartate ring, and it is at the far periplasmic end.

Mouth inner diameter is ~15 Å (Q111 ring, CB–CB 21.1 Å across protomers). That
admits **one α-helix with side chains, not a helical hairpin** — which
independently agrees with the colicin E1/TolC cryo-EM structure (PDB 6WXI),
where colicin plugs TolC as a single-pass folded helix.

**Inter-protomer crown groove.** Loop 2 of one protomer packs against Loop 1 of
its neighbour:

```
A/W333 CB — C/L114 CB = 6.4 Å
A/W333 CB — C/V119 CB = 9.2 Å
A/F326 CB — C/L114 CB = 11.2 Å
```

**Anchor residues by exposure** (Shrake–Rupley on the trimer, Å²):
F326 **149**, K327 126, S328 125, L322 114, N118 112, E321 106, V320 105,
V119 102, R110 100, W333 78.

---

## 1. Three campaigns

| | A — lumen plug | B — macrocycles | C — crown groove |
| --- | --- | --- | --- |
| site | barrel lumen, z +34 to +66 | Loop 2 groove **and** lumen mouth | inter-protomer Loop1/Loop2 groove |
| hotspots (auth) | A123, A312, A311, A132 | A326/A322/A327; A111/A120/A123 | A326, A333, C114, C119 |
| length | 60–110 aa | 8–16 (cyclic) | 55–100 aa |
| precedent | colicin E1/TolC 6WXI; KlebC/TolC | MtrCDE self-inhibitory peptides | anti-Loop2 mAb, MAP vaccine |
| expected function | efflux block | efflux block (mouth) / surface bind (Loop 2) | surface targeting; efflux block unlikely |
| risk | high | medium | low |

Run **C first** (accessibility proven, validates the whole assay cascade),
**B in parallel** (cheap, synthesisable), **A once the cascade works** (highest
payoff, hardest geometry).

**Mechanistic honesty.** Only campaigns that occlude the lumen — A, and the
mouth arm of B — have a route to inhibiting efflux. A Loop 2 surface binder may
bind beautifully and do nothing to efflux, because the loops are not the
conduit. Decide which outcome you are buying before you order genes. If you
want both, fuse a C-type binder to an A-type plug.

---

## 2. Generator settings

From the RXFP1 paper, on the same target with the same experimental readout:

| | BindCraft v1.5.0 | RFdiffusion + ProteinMPNN + AF2 |
| --- | --- | --- |
| designs tested | 48 | 96 |
| **experimental hit rate** | **50%** | **24%** |

So: **BindCraft is the primary generator, RFdiffusion the fold-diversity arm.**

Other settings worth copying verbatim from that paper:
- Only **2–3 hotspot residues** per campaign. Minimal hotspot specification.
- Binder lengths 65–120 aa (blockers), 40–120 aa (conformation-trappers).
- RFdiffusion: 3,000 backbones → AF2 initial-guess **interface pAE < 10** →
  ~200 → **five rounds of sequence–structure recycling** → ~495.
- **Manual curation for fold diversity** before gene synthesis. Not optional —
  it is what stops you ordering 96 variations on one helical bundle.
- Twist genes → pET28a(+)/pET29b(+), His tag placed per design to avoid the
  interface → *E. coli* C41(DE3) → 96-well batch expression, 18 h at 30 °C →
  B-PER lysis → Ni-agarose in filter plates. >80% expression success.

---

## 3. False positives — the filtering cascade

Boltz and AF3-class co-folding models are structure predictors, not binary
binding discriminators; the BindCraft authors tested AF3 as a filter and still
found a large proportion of false positives. On a membrane protein it is worse,
because the model has no bilayer and will happily dock binders onto the
lipid-facing belt.

**Never use Boltz as a generator or primary scorer.** Order of operations:

1. **Geometry** (`scripts/filter_designs.py`) — deterministic, free, run first.
   Rejects membrane-belt and periplasmic docking, and enforces real insertion
   depth for plug designs.
2. **AF2 initial-guess** — interface pAE < 10, the paper's criterion.
3. **Rosetta** — ddG, shape complementarity ≥ 0.62, unsatisfied buried polars ≤ 2.
4. **Boltz-2** — third orthogonal vote only, never alone.
5. **Decoy panel** — rescore every survivor against *E. coli* TolC (1EK9) and
   *P. aeruginosa* OprM. Discard anything scoring as well on a decoy as on
   MtrE. This also generates the selectivity data you need anyway.
6. **Scramble control** — shuffle interface residues, rescore. If scrambles
   score nearly as well, the metric is reading fold quality, not interface.

---

## 4. Experimental cascade

**Target production is the hard part.** Express MtrE with a signal peptide for
OM targeting, purify in DDM or LDAO, then **reconstitute into nanodiscs** —
this is what the colicin E1/TolC work used, and detergent micelles give
artefact-prone BLI. Biotinylate the **scaffold, not MtrE**.

1. BLI on MtrE-nanodisc, 200 s association / 500 s dissociation, 96-well.
2. **Empty-nanodisc counter-screen on every binder.** The single most important
   false-positive control in the project.
3. TolC-nanodisc counter-screen → selectivity.
4. Whole-cell flow cytometry, live FA1090 + WHO K/P/X, **with isogenic ΔmtrE**.
5. Nile red or ethidium bromide efflux/accumulation assay.
6. Antibiotic checkerboard MIC: azithromycin, ciprofloxacin, ceftriaxone,
   erythromycin — matching the MtrCDE self-inhibitory-peptide benchmark
   (2–64-fold potentiation, no human-cell toxicity) so numbers are comparable.
7. Prometheus back-reflection aggregation + CD thermal melt.
8. Cryo-EM of binder–MtrE–nanodisc for the lead.

Two orthogonal primary assays, as in the RXFP1 paper where BLI and competition
correlated at r = −0.55 to −0.72. One assay is not a screen.

*N. gonorrhoeae* is BSL-2. Frame these as **efflux-pump adjuvants**, not
standalone antimicrobials — that is what the precedents actually deliver.

---

## 5. Known caveats

- **LOS shielding.** The crown protrudes only ~7–10 Å above the outer-leaflet
  girdle. Lipooligosaccharide inner core extends ~10–15 Å above the headgroups,
  so the z = 58–62 shell is sterically contested on a live cell. The filter
  flags it rather than rejecting it, because the anti-Loop2 antibody data says
  the region is reachable in practice. Watch for designs that bind purified
  protein and fail on whole cells — that is the signature.
- **Gating state.** 4MT0 is the open state. MtrC binding stabilises it. If a
  plug only engages the open channel, the molecule is conditionally active on
  actively-effluxing cells. Acceptable for an adjuvant, but it changes assay
  design.
- **Proteolysis.** Mini-proteins at a mucosal surface are exposed. This argues
  for the macrocycle arm as the translational route.
- **No experimental MtrCDE assembly structure exists.** Stoichiometry is
  MtrD₃–MtrC₆–MtrE₃. Model by homology to AcrAB–TolC (5O66) if needed, and use
  it only to reason about gating, never as a design target.

---

## 6. Files

```
scripts/prepare_mtre_targets.py   fetch 4MT0, frame, self-check, write targets
scripts/filter_designs.py          geometry filters (run before any scoring)
targets/mtre_sites.json            residue sets, planes, pore profile
targets/target_crown.pdb           campaign C
targets/target_lumen.pdb           campaign A
targets/target_loop2.pdb           campaign B
configs/bindcraft/*.json           three target definitions + strict filters
configs/rfdiffusion/*.sh           diversity arm
configs/rfpeptides/*.sh            macrocycle arm
```

Config files carry `_comment` blocks explaining each choice. Key names are
flagged where they need verifying against your installed version — BindCraft
and RFpeptides have both moved their schemas between releases.
