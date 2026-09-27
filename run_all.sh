#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- end-to-end reproduction.
#   Assumes data_raw/ has been populated by ./download_data.sh.
#   Set NHANES_ROOT to run from another location.
# =============================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
: "${NHANES_ROOT:=$(pwd)}"; export NHANES_ROOT

echo "== [1/17] cleaning_v6.R            (merge NHANES files -> data/analysis_df.rds)"
Rscript cleaning_v6.R

echo "== [2/17] analysis_reanalysis_A.R  (single-metal Cox M1-M3, quartiles, scale metrics)"
Rscript analysis_reanalysis_A.R

echo "== [3/17] analysis_reanalysis_B.R  (Se/Mn, period stability, cause-specific, Se:Hg)"
Rscript analysis_reanalysis_B.R

echo "== [4/17] rcs_fix_contrast.R       (restricted cubic splines, CORRECT contrast variance)"
Rscript rcs_fix_contrast.R

echo "== [5/17] analysis_neversmoker.R   (never-smoker sensitivity analyses)"
Rscript analysis_neversmoker.R

echo "== [6/17] lod_table.R              (limit-of-detection / assay floor table)"
Rscript lod_table.R

echo "== [7/17] rev5_sensitivity.R       (delayed entry, assay floor, E-values, Se-Hg)"
Rscript rev5_sensitivity.R

echo "== [8/17] rev6_A.R / rev6_B.R      (Fine-Gray, stepwise, pack-years, missingness)"
Rscript rev6_A.R
Rscript rev6_B.R

echo "== [9/17] rev7_mi_ipw.R            (multiple imputation + inverse probability weighting)"
Rscript rev7_mi_ipw.R 20 10

echo "== [10/17] rev8_direction_stability.R (direction stability with unified cut-points)"
Rscript rev8_direction_stability.R

echo "== [11/17] rev9_supplementary.R   (alternative trend scores, winsorising, MI variance)"
Rscript rev9_supplementary.R

echo "== [12/17] rev10_review9.R        (Zhang-Yu E-value conversion, unweighted cause-specific Cox)"
Rscript rev10_review9.R

echo "== [13/17] rev11_addenda.R         (pairwise metal product terms; mercury-cancer landmark)"
Rscript rev11_addenda.R

echo "== [14/17] rev12_selection.R      (selection-bias comparison of participants without a usable blood specimen; Table S24)"
Rscript rev12_selection.R

echo "== [15/17] analysis_core / mixture / extra / supp / fish / se_mn.R (table inputs: single_metal_cox, qgcomp, sensitivity, Se/Mn, n-3)"
Rscript analysis_core.R
Rscript analysis_mixture.R
Rscript analysis_extra.R
Rscript analysis_supp.R
Rscript analysis_fish.R
Rscript analysis_se_mn.R

echo "== [16/17] make_figs.py            (Figures 1-5 and Figure S1; the ONLY figure entry point)"
python3 make_figs.py

echo "== [17/17] gen_tables.py / gen_supplementary.py (markdown and DOCX tables)"
python3 gen_tables.py
python3 gen_supplementary.py

echo
echo "ALL DONE -- results in results/, figures in figures/, tables in tables.md"
