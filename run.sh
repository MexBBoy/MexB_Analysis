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
rm -f results/chimerax/* results/figures/*

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
echo "### CAVER cross-check (set SKIP_CAVER=1 to skip; ~25 min)"
if [ "${SKIP_CAVER:-0}" = "1" ]; then
  echo "(skipped)"
else
  python3 scripts/run_caver.py --chains E || echo "(CAVER cross-check failed)"
fi

echo
echo "### section 6 validation"
python3 scripts/validate.py || echo "(validation reported failures - see results/tables/validation.csv)"

echo
echo "### stage 8 - report"
python3 scripts/report.py

echo
echo "### figures"
python3 scripts/figures.py

echo
echo "### poster panels"
python3 scripts/poster_figures.py

echo
echo "### viewer exports (ChimeraX sessions, painted PDBs, 3D viewer)"
python3 scripts/viewer_exports.py
python3 scripts/build_viewer.py

echo
echo "### manuscript draft"
python3 scripts/build_paper.py

echo
echo "### combined workbook"
python3 scripts/make_workbook.py

echo
echo "done. see results/REPORT.md, results/MexB_analysis_results.xlsx,"
echo "          results/figures/, results/chimerax/, results/viewer/,"
echo "          results/paper/mexb_paper.html"
