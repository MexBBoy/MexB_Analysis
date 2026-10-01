#!/usr/bin/env bash
# APPROACH C - extracellular crown / inter-protomer groove, RFdiffusion arm.
#
# Target: the groove where Loop 2 of one protomer meets Loop 1 of the next.
#   A/W333 CB -- C/L114 CB = 6.4 A
#   A/W333 CB -- C/V119 CB = 9.2 A
# Designing across two protomers rather than at one loop roughly doubles the
# buried surface and removes most of the single-loop flexibility problem.
#
# ON SYMMETRY: do not use RFdiffusion's symmetric mode here. That mode builds
# symmetric oligomers de novo; it does not design a C3 binder onto a C3 target,
# and the three-fold site you actually want to engage is the inter-protomer
# groove, which is itself asymmetric in the binder's frame. If you want
# three-fold avidity, design a monomeric groove binder first and then trimerise
# it experimentally by fusing to T4 foldon or a coiled-coil - that is one
# cloning step and it does not compromise the design run.

set -euo pipefail
RFDIFFUSION="${RFDIFFUSION:?set RFDIFFUSION to your RFdiffusion install root}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TARGET="$HERE/targets/target_crown.pdb"
OUT="$HERE/designs/rfd_crown/crown"
mkdir -p "$(dirname "$OUT")"

python "$RFDIFFUSION/scripts/run_inference.py" \
  inference.output_prefix="$OUT" \
  inference.input_pdb="$TARGET" \
  inference.num_designs=3000 \
  'contigmap.contigs=[A110-123/0 A317-333/0 B110-123/0 B317-333/0 C110-123/0 C317-333/0 55-95]' \
  'ppi.hotspot_res=[A326,A333,C114,C119]' \
  denoiser.noise_scale_ca=0 \
  denoiser.noise_scale_frame=0

# The crown protrudes only ~7-10 A above the outer-leaflet aromatic girdle
# (F326 at z=+59.5). Many designs will wrap around and down the barrel into
# what is really lipid and LOS. scripts/filter_designs.py rejects any design
# placing atoms below z=+58; expect it to remove a large fraction of this run.
