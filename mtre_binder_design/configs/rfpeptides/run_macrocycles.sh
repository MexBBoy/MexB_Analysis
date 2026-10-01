#!/usr/bin/env bash
# APPROACH B - cyclic peptides under 40 aa.
#
# TWO SITES, run both. They have very different shapes and you do not yet know
# which is engageable by a macrocycle.
#
#   site 1  LOOP2 GROOVE   shallow, partly polar, centred on the exposed F326.
#                          Classic macrocycle site: a small hydrophobic anchor
#                          in a groove. Hotspots A326, A322, A327.
#
#   site 2  LUMEN MOUTH    the extracellular opening of the barrel, pore radius
#                          7.4 A at z=+62 widening to 9.3 A at z=+66. A 10-14mer
#                          macrocycle is dimensionally a good cork for a 15 A
#                          aperture, and unlike the Loop 2 site a bound peptide
#                          here plausibly blocks efflux rather than just sticking
#                          to the surface. Hotspots A111, A120, A123.
#
# IMPORTANT CORRECTION, carried over from the first draft of this plan: do NOT
# target the D422/D425 aspartate ring. It sits at z=-46, at the PERIPLASMIC tip,
# ~110 A below the extracellular surface and on the wrong side of the outer
# membrane. It is unreachable by anything applied to intact cells.
#
# Permeability is a non-issue for all of this - the sites are extracellular, so
# the usual killer constraint on macrocycles does not apply here.
#
# Set RFPEPTIDES to your RFdiffusion-with-cyclic-support install root.

set -euo pipefail
RFPEPTIDES="${RFPEPTIDES:?set RFPEPTIDES to your RFpeptides install root}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
mkdir -p "$HERE/designs/rfpep"

# --- site 1: Loop 2 groove, 8-14mers -----------------------------------------
for N in 8 10 12 14; do
  python "$RFPEPTIDES/scripts/run_inference.py" \
    inference.output_prefix="$HERE/designs/rfpep/loop2_${N}mer" \
    inference.input_pdb="$HERE/targets/target_loop2.pdb" \
    inference.num_designs=500 \
    inference.cyclic=True \
    inference.cyclize=True \
    "contigmap.contigs=[A315-335/0 B315-335/0 ${N}-${N}]" \
    'ppi.hotspot_res=[A326,A322,A327]' \
    denoiser.noise_scale_ca=0 \
    denoiser.noise_scale_frame=0
done

# --- site 2: lumen mouth, 10-16mers ------------------------------------------
for N in 10 12 14 16; do
  python "$RFPEPTIDES/scripts/run_inference.py" \
    inference.output_prefix="$HERE/designs/rfpep/mouth_${N}mer" \
    inference.input_pdb="$HERE/targets/target_lumen.pdb" \
    inference.num_designs=500 \
    inference.cyclic=True \
    inference.cyclize=True \
    "contigmap.contigs=[A105-133/0 A310-340/0 B105-133/0 B310-340/0 C105-133/0 C310-340/0 ${N}-${N}]" \
    'ppi.hotspot_res=[A111,A120,A123]' \
    denoiser.noise_scale_ca=0 \
    denoiser.noise_scale_frame=0
done

# Downstream (RFpeptides paper protocol):
#   1. ProteinMPNN with cyclic offset, ~48 sequences per backbone
#   2. RoseTTAFold2 with cyclic relative position encoding for validation
#   3. Rosetta FastRelax, score by ddG and macrocycle strain
#   4. keep designs where predicted and designed backbones agree < 1.5 A RMSD
#   5. order 30-60 by SPPS, head-to-tail cyclised, all-L in generation 1
#
# FLAG: the inference.cyclic / inference.cyclize key names differ between the
# public RFpeptides release and the internal version used in the paper. Check
# your install's config before launching 4000 jobs.
#
# For the 20-40 aa regime the user asked about, RFpeptides is the wrong tool -
# it is built for 7-20mers. Use AfCycDesign (cyclic-offset hallucination) or
# RFdiffusion with a 25-40 residue contig plus a designed disulfide instead.
