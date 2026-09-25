#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- end-to-end reproduction (assumes data_raw/ is already populated
#               by ./download_data.sh).
# =============================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

echo "== [1/9] cleaning.R  (merge NHANES files -> data/analysis_df.rds)"
Rscript cleaning.R

echo "== [2/9] analysis_reanalysis_A.R  (single-metal Cox, quartiles, sensitivity,"
echo "        subgroups, Table 1, -> data/analysis_df_plus.RDS)"
Rscript analysis_reanalysis_A.R

echo "== [3/9] analysis_reanalysis_B.R  (Se/Mn, period stability, cause-specific,"
echo "        Se:Hg ratio)"
Rscript analysis_reanalysis_B.R

echo "== [4/9] rcs3.R  (restricted cubic spline curves)"
Rscript rcs3.R

echo "== [5/9] analysis_neversmoker.R  (never-smoker sensitivity analyses)"
Rscript analysis_neversmoker.R

echo "== [6/9] lod_table.R  (limit-of-detection table per cycle)"
Rscript lod_table.R

echo "== [7/9] make_figs.py  (Figures 1-5 + FigS1; PNG 600 dpi / TIFF / EPS)"
python3 make_figs.py

echo "== [8/9] flow_figure.py  (FigS1 flow diagram)"
python3 flow_figure.py

echo "== [9/9] gen_tables.py  (markdown Tables 1-4)"
python3 gen_tables.py

echo
echo "ALL DONE -- results in results/, figures in figures/, tables in tables.md"
