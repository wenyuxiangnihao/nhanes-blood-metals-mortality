# Opposing directions in blood metal mixtures and mortality in US adults

Analysis code for a pooled, nationally representative NHANES study (10 cycles,
1999–2018) of blood **lead, cadmium and mercury** (and, for 2011–2018,
**selenium and manganese**) in relation to **all-cause and cardiovascular
mortality** (NCHS 2019 Linked Mortality File). We use survey-weighted Cox
models, quantile g-computation with a directional decomposition, restricted
cubic splines and a large set of sensitivity / subgroup / cause-specific
analyses.

Key point: because cadmium and lead act in the opposite direction to mercury,
the single net mixture effect conceals both — reporting the positive- and
negative-direction partial effects is more informative than ψ alone.

> **This repository contains code and documentation only. No NHANES data are
> redistributed here** (all source data are public and re-downloaded by the
> scripts below).

---

## 1. Study at a glance

| | |
|---|---|
| Design | Repeated cross-sectional NHANES cycles with mortality linkage |
| Cycles | 1999–2000, 2001–2002, **2003–2004**, 2005–2006 … 2017–2018 (10 cycles) |
| Analytic sample | 33,104 adults ≥20 y (4,260 all-cause, 1,332 CVD deaths) |
| Exposures | Blood Pb (`LBXBPB`), Cd (`LBXBCD`), Hg (`LBXTHG`); Se (`LBXBSE`), Mn (`LBXBMN`) 2011–2018 |
| Outcomes | All-cause (`event_all`) and CVD (`event_cvd`) mortality |
| Main methods | Survey-weighted Cox; quantile g-computation + directional decomposition; RCS |
| Software | R 4.5.3, Python 3 (see §4) |

---

## 2. Repository contents

```
├── .gitignore                   # ignores data_raw/, data/, results/, figures/
├── .zenodo.json                 # Zenodo deposition metadata (the concept DOI follows the latest version)
├── CITATION.cff                 # machine-readable citation metadata (version + concept DOI)
├── LICENSE                      # MIT licence
├── README.md                    # this file
├── R_packages.txt               # R package versions used
├── analysis_core.R              # table inputs: single-metal Cox (M1-M3), quartiles, scale metrics
├── analysis_extra.R             # table inputs: interaction and stratified analyses
├── analysis_fish.R              # table inputs: marine n-3 (EPA+DHA) subsample analyses
├── analysis_mixture.R           # table inputs: quantile g-computation (net effect and directional decomposition)
├── analysis_neversmoker.R       # never-smoker sensitivity analyses (Table S7)
├── analysis_reanalysis_A.R      # single-metal Cox M1-M3, quartiles, scale metrics (Table 2), also writes data/analysis_df_plus.RDS
├── analysis_reanalysis_B.R      # selenium/manganese, period stability, cause-specific, Se:Hg product term
├── analysis_se_mn.R             # table inputs: selenium and manganese
├── analysis_supp.R              # table inputs: additional sensitivity analyses
├── bkmr_local.R                 # exploratory BKMR helper (NOT reported in the manuscript; see the note under 'Running the analysis')
├── bkmr_post.R                  # regenerates the exploratory BKMR summary objects without refitting (NOT reported in the manuscript)
├── bkmr_run.R                   # exploratory BKMR fit (NOT reported in the manuscript)
├── cleaning.R                   # earlier version of the cleaning script, kept for provenance
├── cleaning_v6.R                # FINAL cleaning script: merges the raw files, derives variables, builds data/analysis_df.rds and data/all_merged_df.rds
├── config.R                     # single source of truth for the four directories (reference implementation; see section 5)
├── config.py                    # Python mirror of config.R
├── data_dictionary.csv          # variable dictionary for the derived analysis dataset
├── download_core.py             # download helper used by download_data.sh
├── download_cov.sh              # covariate-file download helper
├── download_covariates.py       # covariate-file download helper (Python)
├── download_data.sh             # downloads every NHANES file + the NCHS 2019 Linked Mortality File into data_raw/
├── download_missing.sh          # re-fetches any file that failed the first download
├── explore_vars.R               # exploratory variable inspection
├── flow_figure.py               # participant flow diagram (Figure S1)
├── gen_supplementary.py         # builds submission/supplementary_tables.docx from tables.md
├── gen_tables.py                # builds tables.md from the results/ files (writes tables_generated.md for comparison)
├── lod_table.R                  # limit-of-detection table by cycle (Table S8)
├── make_figs.py                 # Figures 1-5 and Figure S1 (png/tiff/eps). THE ONLY FIGURE ENTRY POINT
├── md_to_docx.py                # renders the markdown sources to the submission .docx files
├── ph_test_10cyc.R              # proportional-hazards tests (Table S2)
├── plot_figures.py              # earlier figure script, kept for provenance
├── qgcss_splitsample.R          # sample-splitting analysis (Table S19)
├── rcs3.R                       # restricted cubic splines, earlier implementation
├── rcs_fix_contrast.R           # restricted cubic splines with the correct contrast variance (Figures 4-5)
├── requirements.txt             # Python package versions used
├── rev10_review9.R              # Zhang-Yu E-value conversion; unweighted cause-specific column (Table S14)
├── rev11_addenda.R              # pairwise product terms (Table S22); mercury-cancer landmark analysis (Table S23)
├── rev12_selection.R            # selection-bias comparison of participants excluded for missing blood metals (Table S24)
├── rev5_sensitivity.R           # subsample weights, survey-cycle adjustment, delayed entry, Se-Hg, assay floor, E-values
├── rev6_A.R                     # additional sensitivity analyses (Table S17)
├── rev6_B.R                     # additional sensitivity analyses (Table S17)
├── rev7_mi_ipw.R                # multiple imputation and inverse probability weighting (Table S20)
├── rev7_mi_ipw_summarize.R      # summarises the imputation results without refitting
├── rev8_direction_stability.R   # stability of the direction assignment (Table S13)
├── rev9_supplementary.R         # winsorising, alternative trend tests, variance decomposition (Table S21)
├── run_all.sh                   # end-to-end pipeline (17 steps)
├── superscript_citations.py     # converts in-text citations to the superscript style used by the journal
└── test_bkmr.R                  # smoke test for the exploratory BKMR fit
```

Output directories (`data/`, `results/`, `figures/`) are created automatically
by `config.R` / `config.py` and are **not** committed.

---

## 3. Data sources and download

All inputs are **public**. The script `download_data.sh` fetches them into
`data_raw/` (override with `RAW_DIR=... ./download_data.sh`).

### 3.1 NHANES examination / questionnaire / laboratory files

URL pattern:

```
https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{startyear}/DataFiles/{FILE}.XPT
```

Per cycle (`startyear` = 1999, 2001, 2003, …, 2017; suffix `""`, `_B`, `_C`, …
`_J`), the following files are downloaded:

| Purpose | File stem (add the cycle suffix) |
|---|---|
| Blood metals (Pb, Cd, Hg, Se, Mn) | `LAB06` (1999), `L06_B` (2001), `L06BMT_C` (2003), `PBCD_D` … `PBCD_J` (2005–2017) |
| Serum cotinine | inside `LAB06`/`L06_B` (1999–2001), `L06COT_C` (2003), `COT_D` … `COT_J` (2005–2017) |
| Demographics (age, sex, race, education, PIR, weights, PSU, strata, pregnancy) | `DEMO` + suffix |
| Body measures | `BMX` + suffix |
| Blood pressure | `BPX` + suffix |
| Smoking | `SMQ` + suffix |
| Diabetes | `DIQ` + suffix |
| Hypertension | `BPQ` + suffix |
| Self-reported CVD/cancer | `MCQ` + suffix |
| Alcohol | `ALQ` + suffix |
| Physical activity | `PAQ` + suffix |
| Serum creatinine (eGFR) | `LAB18` (1999), `L40_B` (2001), `L40_C` (2003), `BIOPRO_D` … `BIOPRO_J` |
| Total cholesterol + HDL | `LAB13` (1999), `L13_B` (2001), `L13_C` (2003), `TCHOL_*` + `HDL_*` (2005–2017) |
| *(optional)* alternate heavy-metal files | `L06HM_B` (2001), `L06HM_C` (2003) |
| *(optional)* marine n-3 fatty acids | `FAS_G` (2011–2012), `FAS_H` (2013–2014) |

### 3.2 NCHS 2019 Linked Mortality File (public-use, fixed-width)

```
https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/
```

Downloaded files (one per cycle):

```
NHANES_1999_2000_MORT_2019_PUBLIC.dat ... NHANES_2017_2018_MORT_2019_PUBLIC.dat
```

Only the public-use LMF is used, so the linked mortality folder above is
sufficient. Documentation, the variable layout and the official SAS/R read-in
program are at
<https://www.cdc.gov/nchs/data-linkage/mortality-public.htm>. The columns read
by `cleaning.R` are `SEQN` (1–6), `eligstat` (15), `mortstat` (16),
`ucod_leading` (17–19), `permth_int` (43–45) and `permth_exm` (46–48).

> Note: the linked-mortality files can also be obtained through the NCHS
> Research Data Center / the "linked mortality public-use files" landing page;
> the 2019 linkage covers deaths through 31 Dec 2019 and is what this code uses.

---

## 4. Software requirements

**R 4.5.3** with:

* `survey` — survey-weighted Cox models (`svycoxph`) and design objects
* `survival` — `Surv`
* `splines` — `ns()` for restricted cubic splines
* `foreign` — `read.xport()` for NHANES `*.XPT`
* `haven` — alternative XPT reader (optional)
* `bkmr` — BKMR (**not used in the manuscript**; needed only to re-run the exploratory script, installed separately)

```r
install.packages(c("survey","survival","splines","foreign","haven"))
# `bkmr` is NOT needed to reproduce the manuscript; install it only if you want to
# re-run the exploratory (unreported) fit:  install.packages("bkmr")
```

**Python 3** with:

```bash
pip install matplotlib numpy pandas Pillow
```

(`Pillow` is only needed to write the LZW-compressed TIFFs.)

---

## 5. Configuration

Every script resolves its own project root, so the repository runs from anywhere:

* the environment variable `NHANES_ROOT` is used if it is set;
* otherwise the directory that contains the script is used (for R scripts this is
  read from `commandArgs()`; for Python scripts from `__file__`);
* R scripts additionally prepend `../Rlibs` to `.libPaths()` when that directory
  exists, which is convenient for a self-contained library but is optional.

`config.R` and `config.py` show the same four directories and their environment
overrides (`RAW_DIR`, `DATA_DIR`, `RESULTS_DIR`, `FIG_DIR`, all defaulting to
subdirectories of the project root). They are a reference implementation rather
than something every script sources.

To run against other directories:

```bash
NHANES_ROOT=/mnt/nhanes Rscript analysis_reanalysis_B.R
```

---
## 6. Running the analysis

```bash
# 0) get the data (network required; creates ./data_raw)
./download_data.sh

# 1) run everything in order
./run_all.sh
```

`run_all.sh` runs these 17 steps with `set -euo pipefail`:

```
 1  cleaning_v6.R             -> data/analysis_df.rds (+ data/exclusion_counts.csv)
 2  analysis_reanalysis_A.R   -> results/*.csv (+ data/analysis_df_plus.RDS)
 3  analysis_reanalysis_B.R   -> results/*.csv
 4  rcs_fix_contrast.R        -> results/rcs_curves*.csv  (correct contrast variance)
 5  analysis_neversmoker.R    -> results/neversmoker_*.csv, cadmium_clrd_by_smoking_*.csv
 6  lod_table.R               -> results/lod_by_cycle_10cyc.csv
 7  rev5_sensitivity.R        -> results/rev5_*.csv
 8  rev6_A.R, rev6_B.R        -> results/rev6_*.csv
 9  rev7_mi_ipw.R             -> results/rev7_mi_ipw_summary.csv (Table S20)
10  rev8_direction_stability.R-> results/rev8_direction_stability.csv (Table S13)
11  rev9_supplementary.R      -> results/rev9_*.csv (Table S21)
12  rev10_review9.R           -> results/rev10_*.csv (E-value conversion; Table S14 column)
13  rev11_addenda.R          -> results/rev11_*.csv (Tables S22-S23)
14  rev12_selection.R        -> results/rev12_tableS24.* (Table S24)
15  analysis_core / mixture / extra / supp / fish / se_mn.R -> results/*.csv (table inputs)
16  make_figs.py              -> figures/fig1..fig5, figS1 (png/tiff/eps).  THE ONLY FIGURE ENTRY POINT
17  gen_tables.py             -> tables.md  ;  gen_supplementary.py -> submission/supplementary_tables.docx
```

Each script also runs on its own (for example `Rscript rev10_review9.R`) as long as
the preceding datasets exist. No script writes a figure that another script also
writes: `make_figs.py` is the single figure entry point.

### BKMR (exploratory script only; not used in the manuscript)

**Status: exploratory only, not reported in the manuscript.** The `bkmr_*.R`
scripts are kept in this repository for full transparency, but the BKMR fit is
**not** part of the manuscript or of `run_all.sh`: unweighted mixture analyses
of complex survey data such as NHANES cannot incorporate the sampling weights,
strata and clusters, and the journal's guide for authors discourages
submissions that analyse chemical mixtures in complex sampling designs with
unweighted methods. In the manuscript, mixture effects are therefore estimated
with survey-weighted quantile g-computation with directional decomposition, and
the non-linear mixture analysis is deliberately omitted. Anyone wishing to
reproduce the exploratory fit needs the `bkmr` package; computational cost grows
roughly as O(N^3), so it was fitted on a random subsample of 800 participants.

---

## 7. Expected outputs

| Path | Description |
|---|---|
| `data/analysis_df.rds` | cleaned 10-cycle analysis dataset |
| `data/analysis_df_plus.RDS` | adds physical-activity / alcohol / quartile variables |
| `data/exclusion_counts.csv` | STROBE-style exclusion chain |
| `results/single_metal_cox*.csv` | per-SD HRs, models M1–M3 |
| `results/single_metal_quartiles*.csv` | quartile HRs, P-trend, P-non-linearity |
| `results/qgcomp_*` | mixture ψ + component weights (incl. directional) |
| `results/rcs_curves*.csv`, `rcs_se_mn*.csv` | spline curves |
| `results/selenium_manganese*.csv` | Se / Mn associations, nadir |
| `results/sensitivity*.csv`, `subgroup_*.csv` | sensitivity & interaction analyses |
| `results/cause_specific*.csv`, `neversmoker_*.csv` | cause-specific / never-smoker |
| `results/lod_by_cycle_10cyc.csv` | limit-of-detection table |
| `figures/fig1_forest.{png,tiff,eps}` … `figS1_flow.{png,tiff,eps}` | manuscript figures (600 dpi / vector) |
| `tables.md` | Tables 1–4 (markdown) |

---

## 8. Reproducibility notes

* **Survey design.** All models use `svydesign(id=~SDMVPSU, strata=~SDMVSTRA,
  weights=~wt, nest=TRUE)` with `wt = WTMEC2YR / 10` (10 pooled cycles) and
  `options(survey.lonely.psu="adjust")`.
* **2003–2004.** Unlike many pooled analyses, this study includes the 2003–2004
  cycle (NCHS file `L06BMT_C`), giving 10 cycles.
* **Limit of detection.** A per-cycle LOD summary is produced by `lod_table.R`;
  no values are imputed with a fixed constant in the main models.
* **BKMR** does not support survey weights; the BKMR fit is exploratory only.
* **Table generator.** `gen_tables.py` was adapted from the earlier 9-cycle
  version to the final 10-cycle results schema (NHANES metal codes, the
  `HR_overall` mixture column, the wide selenium/manganese file, and selection
  of the fully adjusted M3 row). These are label/column mappings only — every
  number in `tables.md` reproduces the manuscript tables exactly.
* **No data are shipped.** `data_raw/`, `data/`, `results/` and `figures/` are
  generated locally and should not be committed.

---

## 9. License and citation

Released under the **MIT License** — see `LICENSE` (© 2026 Yuxiang Wen).
Machine-readable citation metadata are in `CITATION.cff` and `.zenodo.json`.

If you use this code, please cite:

> Wen Y. *Opposing directions in blood metal mixtures and mortality in US
> adults.* (2026).

**Contact:** Yuxiang Wen — Department of Cardiology, The First Affiliated
Hospital of Yangtze University, Jingzhou, Hubei, China.

## Repository URL

The repository URL appears in `CITATION.cff` (`repository-code:`) and in the Code
availability statement of the manuscript. If you fork this repository, replace
`the repository URL of your fork` with your own URL; creating a
Zenodo release for a GitHub tag then provides the archival DOI.

## BKMR fit (exploratory; not part of the manuscript)

`bkmr_run.R` / `bkmr_local.R` fit the exploratory Bayesian kernel machine regression
(reproducing it requires the `bkmr` package; computational cost grows roughly as O(N^3),
so it was fitted on a random subsample of 800 participants), and `bkmr_post.R`
regenerates the summary objects without refitting. These scripts are retained here for
transparency but produce **no result reported in the manuscript** — see the note under
"Running the analysis" for the reason (survey-weighted mixture analysis is required for
NHANES-type designs).

## Additional analyses

`rev5_sensitivity.R` reproduces the round-3 supplementary analyses (subsample weights
WTSH2YR, survey-cycle adjustment, delayed entry at 48/60 months, selenium–mercury
interaction, assay-floor exclusion and E-values). `ph_test_10cyc.R` reproduces the
proportional-hazards tests in Table S2.

`qgcss_splitsample.R` implements the sample-splitting analysis reported in Table S19
(quantile g-computation with sample splitting: the direction of each metal is determined
in one random half of the data and the partial effects are estimated in the other half,
200 splits).

`rev7_mi_ipw.R` implements the missing-data analyses reported in Table S20. It imputes the
missing covariates by chained equations (M = 20 imputations, 10 iterations; linear,
proportional-odds, multinomial and logistic models as appropriate), analyses each imputed
dataset in the full eligible sample (n = 36,743) with the survey-weighted Cox models and
combines the results by Rubin's rules, and it also re-analyses the complete-case sample with
inverse probability of inclusion weights. `rev7_mi_ipw_summarize.R` regenerates
`results/rev7_mi_ipw_summary.csv` from the saved intermediate objects without refitting.

`rev8_direction_stability.R` recomputes the direction-assignment stability analysis
(Table S13) with the quartile cut-points defined within each mixture's own sample, so that
the full-sample coefficients correspond exactly to the directional sums of Table 4; it
supersedes the corresponding part of `rev6_A.R`.

`rcs_fix_contrast.R` recomputes the restricted cubic spline confidence bands from the
variance of the spline contrast, which is the correct variance for a dose-response curve
relative to its reference point (Figures 4 and 5).

`rev11_addenda.R` fits the pairwise product terms between the three core metals
(Table S22) and the mercury–cancer landmark analyses with delayed entry at 24, 48 and
60 months (Table S23).

`rev12_selection.R` characterises the participants who were excluded because no usable
blood specimen was available for lead, cadmium or mercury (n = 8,640) and compares them
with the analytic sample (n = 33,104) on the covariates used in the main models
(Table S24). It reads `data/all_merged_df.rds`, the full merged frame written by
`cleaning_v6.R` before the exclusion chain. Comparisons are unweighted, because the MEC
examination weight exists only for participants who attended the examination and the
excluded group includes non-attenders; a survey-weighted sensitivity comparison is
restricted to participants with a positive weight.

## Citation

If you use this code, please cite the archived version:

> Wen Y. *Opposing directions in blood metal mixtures and mortality in US adults.* Zenodo.
> https://doi.org/10.5281/zenodo.22975596 (concept DOI; resolves to the latest version).

Source repository: https://github.com/wenyuxiangnihao/nhanes-blood-metals-mortality
