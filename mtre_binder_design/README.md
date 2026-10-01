# MtrE extracellular binder design

De novo binders against **MtrE**, the outer-membrane channel of the
*N. gonorrhoeae* MtrCDE efflux pump. Extracellular sites only — nothing here
needs to cross the outer membrane.

Read **[PROTOCOL_MtrE.md](PROTOCOL_MtrE.md)** first. It holds the ground-truth
geometry, the numbering convention, and the reasoning behind every hotspot.

## Quick start

```bash
pip install numpy biopython

# Fetch 4MT0, frame it, self-check, write the three targets.
python3 scripts/prepare_mtre_targets.py --outdir targets

# Campaign C (start here) — BindCraft
bindcraft --settings configs/bindcraft/target_crown_loops.json \
          --filters  configs/bindcraft/filters_mtre_strict.json \
          --advanced <your_install>/settings_advanced/default_4stage_multimer.json

# Geometry filter before spending any scoring compute
python3 scripts/filter_designs.py --mode crown designs/mtre_crown/*.pdb
```

## The one thing that will break your run

4MT0 uses **precursor numbering** (auth = mature + 20). The aspartate ring is
**D402/D405** in the papers and **D422/D425** in the coordinates — same
residues. Everything in this directory uses auth numbering, matching the PDB.

## Two corrections to the earlier plan

Both came out of measuring the structure rather than reasoning from the
literature:

1. **Do not target the aspartate ring.** It sits at z = −46, at the periplasmic
   tip, ~110 Å below the extracellular surface. It is the channel's only
   constriction, which makes it an appealing "cork" site, and it is completely
   unreachable from outside the cell.

2. **Do not use RFdiffusion symmetric mode for a C3 binder.** That mode builds
   symmetric oligomers de novo; it does not design a C3 binder onto a C3
   target. The site you want is the *inter-protomer* groove, which is
   asymmetric in the binder's frame. For three-fold avidity, design a monomeric
   groove binder and trimerise it experimentally with a T4 foldon fusion.

## What the structure actually says

- The β-barrel lumen is an **open funnel** — pore radius 7.4–9.3 Å from the
  extracellular mouth at z = +66 down to z = +34, with **no constriction**. A
  binder can insert ~32 Å into the channel from outside.
- The mouth is ~15 Å across: **one helix fits, a hairpin does not**. This
  independently matches the colicin E1/TolC structure (6WXI), where colicin
  plugs TolC as a single-pass folded helix.
- **F326 is the best anchor on the surface** — an exposed Phe with 149 Å² SASA,
  sitting at the Loop1/Loop2 inter-protomer seam alongside W333.
- The crown protrudes only ~7–10 Å above the lipid. LOS will contest part of
  that; the anti-Loop2 antibody data says it is reachable anyway.

## Campaigns

| | site | hotspots (auth) | length | risk |
| --- | --- | --- | --- | --- |
| **A** lumen plug | barrel lumen | A123, A312, A311, A132 | 60–110 aa | high |
| **B** macrocycles | Loop 2 groove + lumen mouth | A326/A322/A327; A111/A120/A123 | 8–16 cyclic | medium |
| **C** crown groove | Loop1/Loop2 seam | A326, A333, C114, C119 | 55–100 aa | low |

Run C first, B in parallel, A once the assay cascade is validated.
