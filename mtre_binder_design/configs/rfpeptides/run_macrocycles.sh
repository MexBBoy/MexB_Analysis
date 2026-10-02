#!/usr/bin/env bash
# APPROACH B - cyclic peptides under 40 aa, via the RFpeptides protocol.
#
# CORRECTION to the first version of this file: RFpeptides is NOT a separate
# install. It is built into mainline RFdiffusion, and the flags are
# `inference.cyclic=True` plus `inference.cyc_chains`, not the
# `inference.cyclize` that earlier draft guessed at. Verified against
# RosettaCommons/RFdiffusion at 86507b6 (examples/design_macrocyclic_binder.sh
# and the "Macrocyclic peptide design with RFpeptides" section of the README).
#
# Two things differ from classic binder design and both matter:
#   1. The PEPTIDE CONTIG COMES FIRST. `inference.cyc_chains='a'` cyclises
#      output chain A, so the diffused peptide has to be the first contig
#      segment. Put the target after it. This is the reverse of the ordering in
#      run_lumen_plug.sh and run_crown_groove.sh, and getting it backwards
#      silently cyclises part of your target.
#   2. `diffuser.T=50` and `--config-name base` are both part of the published
#      protocol. Keep them.
#
# TWO SITES, run both - they have different shapes and it is not yet known
# which is engageable by a macrocycle.
#
#   site 1  LUMEN MOUTH    the extracellular opening of the barrel, pore radius
#                          7.4 A at z=+62 widening to 9.3 A at z=+66. A 10-16mer
#                          is dimensionally a good cork for a 15 A aperture, and
#                          unlike the Loop 2 site a peptide bound here plausibly
#                          blocks efflux rather than just decorating the surface.
#
#   site 2  LOOP2 GROOVE   shallow, partly polar, centred on the exposed F326.
#                          A classic macrocycle site - small hydrophobic anchor
#                          in a groove - but binding here is not expected to
#                          stop efflux on its own.
#
# Do NOT target the D422/D425 aspartate ring. It sits at z=-46, at the
# PERIPLASMIC tip, ~110 A below the extracellular surface and on the far side of
# the outer membrane. It is unreachable by anything applied to intact cells.
#
# Permeability is a non-issue for all of this: the sites are extracellular, so
# the constraint that usually kills macrocycle programmes does not apply.

set -euo pipefail
RFDIFFUSION="${RFDIFFUSION:?set RFDIFFUSION to your RFdiffusion install root}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$HERE/designs/rfpep"
mkdir -p "$OUT"

# --- site 1: lumen mouth, 10-16mers -----------------------------------------
# Hotspots line the mouth: Q111, S120, Y123 (4MT0 auth numbering).
for N in 10 12 14 16; do
  python "$RFDIFFUSION/scripts/run_inference.py" \
    --config-name base \
    inference.output_prefix="$OUT/mouth_${N}mer" \
    inference.input_pdb="$HERE/targets/target_lumen.pdb" \
    inference.num_designs=500 \
    inference.cyclic=True \
    inference.cyc_chains='a' \
    diffuser.T=50 \
    "contigmap.contigs=[${N}-${N} A105-133/0 A310-340/0 B105-133/0 B310-340/0 C105-133/0 C310-340/0]" \
    ppi.hotspot_res=[\'A111\',\'A120\',\'A123\']
done

# --- site 2: Loop 2 groove, 8-14mers ----------------------------------------
# Hotspots are the exposed anchors: F326 (SASA 149 A^2), L322, K327.
for N in 8 10 12 14; do
  python "$RFDIFFUSION/scripts/run_inference.py" \
    --config-name base \
    inference.output_prefix="$OUT/loop2_${N}mer" \
    inference.input_pdb="$HERE/targets/target_loop2.pdb" \
    inference.num_designs=500 \
    inference.cyclic=True \
    inference.cyc_chains='a' \
    diffuser.T=50 \
    "contigmap.contigs=[${N}-${N} A315-335/0 B315-335/0]" \
    ppi.hotspot_res=[\'A326\',\'A322\',\'A327\']
done

# Downstream (RFpeptides paper protocol):
#   1. ProteinMPNN with cyclic offset, ~48 sequences per backbone
#   2. RoseTTAFold2 with cyclic relative position encoding for validation
#   3. Rosetta FastRelax; score by ddG and macrocycle strain
#   4. keep designs where predicted and designed backbones agree < 1.5 A RMSD
#   5. scripts/filter_designs.py, then order 30-60 by SPPS, head-to-tail
#      cyclised, all-L in generation 1
#
# For the 20-40 aa regime, RFpeptides is the wrong tool - it is built for
# 7-20mers. Use AfCycDesign (cyclic-offset hallucination), or RFdiffusion with
# a 25-40 residue contig plus a designed disulfide, instead.
