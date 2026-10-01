#!/usr/bin/env bash
# APPROACH A - extracellular lumen plug, RFdiffusion arm.
#
# Runs alongside the BindCraft arm (configs/bindcraft/target_lumen_plug.json).
# The RXFP1 paper got 24% experimental hit rate from RFdiffusion vs 50% from
# BindCraft on the same target, so treat this as the fold-diversity arm: it
# explores topologies BindCraft will not, at a lower per-design success rate.
#
# Target geometry (measured from 4MT0, see targets/mtre_sites.json):
#   extracellular mouth        z = +66
#   deepest reachable point    z = +34
#   pore radius throughout     7.4 - 9.3 A
#   mouth inner diameter       ~15 A  -> a SINGLE helix fits, a hairpin does not
#
# Hotspots are placed at four depths so the sampler has to build something that
# actually inserts. Hotspot chain IDs refer to target_lumen.pdb (A/B/C trimer).
#
# Set RFDIFFUSION to your install root before running.

set -euo pipefail
RFDIFFUSION="${RFDIFFUSION:?set RFDIFFUSION to your RFdiffusion install root}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TARGET="$HERE/targets/target_lumen.pdb"
OUT="$HERE/designs/rfd_lumen_plug/plug"
mkdir -p "$(dirname "$OUT")"

# ---------------------------------------------------------------------------
# Run 1: unconstrained binder design, 60-100 residues.
# 3000 backbones, matching the scale the RXFP1 paper used per campaign.
# ---------------------------------------------------------------------------
python "$RFDIFFUSION/scripts/run_inference.py" \
  inference.output_prefix="$OUT" \
  inference.input_pdb="$TARGET" \
  inference.num_designs=3000 \
  'contigmap.contigs=[A105-133/0 A310-340/0 B105-133/0 B310-340/0 C105-133/0 C310-340/0 60-100]' \
  'ppi.hotspot_res=[A123,A312,A311,A132]' \
  denoiser.noise_scale_ca=0 \
  denoiser.noise_scale_frame=0

# ---------------------------------------------------------------------------
# Run 2 (recommended): fold-conditioned, forcing a long single helix.
# The colicin E1 / TolC structure (PDB 6WXI) shows colicin plugging TolC as a
# single-pass folded helix, and the 15 A mouth here says the same thing. Giving
# RFdiffusion that prior up front beats filtering for it afterwards.
#
# Build the secondary-structure / block-adjacency inputs first:
#   python $RFDIFFUSION/helper_scripts/make_secstruc_adj.py \
#       --input_pdb scaffolds/single_helix_40.pdb --out_dir scaffolds/plug_ss
# where single_helix_40.pdb is an ideal 40-residue poly-ALA helix (~60 A long,
# enough to span mouth to floor with the cap sitting above the crown).
# ---------------------------------------------------------------------------
# python "$RFDIFFUSION/scripts/run_inference.py" \
#   --config-name symmetry \
#   inference.output_prefix="${OUT}_ss" \
#   inference.input_pdb="$TARGET" \
#   inference.num_designs=1000 \
#   scaffoldguided.scaffoldguided=True \
#   scaffoldguided.target_pdb=True \
#   scaffoldguided.target_path="$TARGET" \
#   scaffoldguided.scaffold_dir="$HERE/scaffolds/plug_ss" \
#   'ppi.hotspot_res=[A123,A312,A311,A132]' \
#   scaffoldguided.target_ss="$HERE/scaffolds/plug_ss/target_ss.pt" \
#   scaffoldguided.target_adj="$HERE/scaffolds/plug_ss/target_adj.pt" \
#   denoiser.noise_scale_ca=0 \
#   denoiser.noise_scale_frame=0

# ---------------------------------------------------------------------------
# Downstream, matching the RXFP1 paper's protocol:
#   1. ProteinMPNN on each backbone (8 seqs/backbone, soluble weights)
#   2. AF2 initial-guess; keep interface pAE < 10
#   3. five rounds of sequence-structure recycling on survivors
#   4. scripts/filter_designs.py  (membrane-belt + decoy-panel consensus)
#   5. manual curation for fold diversity before ordering genes
# ---------------------------------------------------------------------------
