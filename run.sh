#!/usr/bin/env bash
# Regenerate everything under results/ from structures/ .
# Safe to run from a clean work/. Nothing under results/ is hand-edited.
set -euo pipefail
cd "$(dirname "$0")"

STEP="${TUNNEL_STEP:-0.6}"

mkdir -p work results/tables results/figures results/chimerax

# snapshot the previous tables so the report can say what moved
rm -rf work/prev_tables
if [ -d results/tables ]; then cp -r results/tables work/prev_tables; fi
rm -f results/tables/*.csv results/tables/flags_*.txt

echo "### stages 1-5, plus pocket composition"
python3 scripts/mexb_analysis.py all

echo
echo "### stage 6 - pockets and cavities"
python3 scripts/pockets.py

echo
echo "### stage 7 - tunnels (trimer, grid step ${STEP} A)"
python3 scripts/tunnels.py --step "${STEP}"

echo
echo "### switch-loop gate diagnostic"
python3 scripts/switch_gate.py

echo
echo "### section 6 validation"
python3 scripts/validate.py || echo "(validation reported failures - see results/tables/validation.csv)"

echo
echo "### stage 8 - report"
python3 scripts/report.py

echo
echo "done. see results/REPORT.md"
