# Opposing directions in blood metal mixtures and mortality in US adults

Analysis code for a pooled, nationally representative NHANES study (10 cycles,
1999–2018) of blood **lead, cadmium and mercury** (and, for 2011–2018,
**selenium and manganese**) in relation to **all-cause and cardiovascular
mortality** (NCHS 2019 Linked Mortality File). We use survey-weighted Cox
models, quantile g-computation with a directional decomposition, exploratory
Bayesian kernel machine regression (BKMR), restricted cubic splines and a large
set of sensitivity / subgroup / cause-specific analyses.

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
| Main methods | Survey-weighted Cox; quantile g-computation + directional decomposition; BKMR; RCS |
| Software | R 4.5.3, Python 3 (see §4) |

---

## 2. Repository contents

```
.
├── config.R                 # single source of truth for all paths (sourced by every R script)
├── config.py                # Python mirror of config.R (same env-var overrides)
├── download_data.sh         # download every NHANES file + NCHS 2019 Linked Mortality File
├── cleaning.R               # merge raw NHANES files, derive variables, build analysis datasets
├── analysis_reanalysis_A.R  # single-metal Cox (M1–M3), quartiles/P-trend, scale metrics,
│                            #   sensitivity analyses, subgroups, Table 1  -> analysis_df_plus.RDS
├── analysis_reanalysis_B.R  # Se & Mn, period stability, cause-specific, Se:Hg ratio
├── rcs3.R                   # restricted cubic spline dose–response curves
├── analysis_neversmoker.R   # never-smoker / cotinine-restricted sensitivity analyses
├── lod_table.R              # limit-of-detection summary per cycle
├── make_figs.py             # Figures 1–5 + Figure S1 (PNG 600 dpi / TIFF 600 dpi / vector EPS)
├── flow_figure.py           # Figure S1 flow diagram (PNG/TIFF/EPS)
├── cleaning_v6.R            # final cleaning script used for the published analysis
├── rev5_sensitivity.R       # round-3 sensitivity analyses (subsample weights, delayed entry, ...)
├── rev6_A.R, rev6_B.R       # direction-stability replicates, Fine–Gray, stepwise adjustment, ...
├── qgcss_splitsample.R      # quantile g-computation with sample splitting (Table S19)
├── rev7_mi_ipw.R            # multiple imputation + inverse probability weighting (Table S20)
├── rev7_mi_ipw_summarize.R  # re-summarises the MI/IPW output without refitting
├── rcs_fix_contrast.R       # restricted cubic spline contrast variance (Figures 4 and 5)
├── gen_tables.py            # markdown Tables 1–4 from results/*.csv
├── run_all.sh               # runs cleaning -> analyses -> figures -> tables (set -e)
├── data_dictionary.csv      # variable name -> meaning / unit
├── CITATION.cff, .zenodo.json, LICENSE (MIT)
└── README.md
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
* `bkmr` — BKMR (optional; only for the exploratory fit)

```r
install.packages(c("survey","survival","splines","foreign","haven","bkmr"))
```

**Python 3** with:

```bash
pip install matplotlib numpy pandas Pillow
```

(`Pillow` is only needed to write the LZW-compressed TIFFs.)

---

## 5. Configuration

`config.R` (and `config.py`) define exactly four directories, each overridable
by an environment variable so the code can run anywhere:

| Variable | Default | Contents |
|---|---|---|
| `RAW_DIR` | `./data_raw` | downloaded NHANES `.XPT` + `.dat` files |
| `DATA_DIR` | `./data` | cleaned `analysis_df.rds`, `analysis_df_plus.RDS` |
| `RESULTS_DIR` | `./results` | all result tables (`*.csv`) |
| `FIG_DIR` | `./figures` | all figures (`*.png`, `*.tiff`, `*.eps`) |

Every script loads these automatically (R scripts `source(config.R)`, Python
scripts `from config import ...`). To run on a cluster, e.g.:

```bash
RAW_DIR=/mnt/nhanes RESULTS_DIR=/fast/res Rscript analysis_reanalysis_B.R
```

`HM_CONFIG` can point to a `config.R` located elsewhere.

---

## 6. Running the analysis

```bash
# 0) get the data (network required; creates ./data_raw)
./download_data.sh

# 1) run everything in order
./run_all.sh
```

`run_all.sh` executes, with `set -e`:

```
download_data.sh (once)  →
cleaning.R               → data/analysis_df.rds        (+ data/exclusion_counts.csv)
analysis_reanalysis_A.R  → results/*.csv               (+ data/analysis_df_plus.RDS)
analysis_reanalysis_B.R  → results/*.csv
rcs3.R                   → results/rcs_curves*.csv
analysis_neversmoker.R   → results/neversmoker_*.csv, cadmium_clrd_by_smoking_*.csv
lod_table.R              → results/lod_by_cycle_10cyc.csv
make_figs.py             → figures/fig1..fig5, figS1 (png/tiff/eps)
flow_figure.py           → figures/figS1_flow (png/tiff/eps)
gen_tables.py            → tables.md
```

Each script is also runnable on its own (e.g. `Rscript rcs3.R`), as long as the
preceding datasets exist. All scripts assume the working directory is the
repository root, but paths are resolved from `config.R` so they work from
anywhere.

### BKMR (optional)

The exploratory BKMR fit requires the `bkmr` package and is computationally
heavy (a random subsample of 800 participants is used in the paper). The
summary figures (`bkmr_overall`, `bkmr_single`) can be redrawn from the
`results/bkmr_*_summaries.csv` tables once the fit has been run; this step is
not part of `run_all.sh`.

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

Released under the **MIT License** — see `LICENSE` (© 2025 Yuxiang Wen).
Machine-readable citation metadata are in `CITATION.cff` and `.zenodo.json`.

If you use this code, please cite:

> Wen Y. *Opposing directions in blood metal mixtures and mortality in US
> adults.* (2025).

**Contact:** Yuxiang Wen — Department of Cardiology, The First Affiliated
Hospital of Yangtze University, Jingzhou, Hubei, China.

## Repository URL

The repository URL appears in `CITATION.cff` (`repository-code:`) and in the Code
availability statement of the manuscript. If you fork this repository, replace
`https://github.com/your-org/nhanes-blood-metals-mortality` with your own URL; creating a
Zenodo release for a GitHub tag then provides the archival DOI.

## BKMR fit

The exploratory Bayesian kernel machine regression (BKMR) fit is provided in
`bkmr_run.R` / `bkmr_local.R` (run in R with the `bkmr` package; computational cost grows
roughly as O(N^3), so it was fitted on a random subsample of 800 participants). The
manuscript figures for BKMR (Figures S2 and S3) are regenerated from
`results/bkmr_*_summaries.csv` by `bkmr_post.R`, which does not require the `bkmr` package.

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

`rcs_fix_contrast.R` recomputes the restricted cubic spline confidence bands from the
variance of the spline contrast, which is the correct variance for a dose-response curve
relative to its reference point (Figures 4 and 5).
