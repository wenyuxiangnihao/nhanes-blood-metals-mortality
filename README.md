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

## 2. Repository contents

```
├── CITATION.cff                  
├── LICENSE                       
├── README.md                     
├── R_packages.txt                
├── analysis_core.R               
├── analysis_extra.R              
├── analysis_fish.R               
├── analysis_mixture.R            
├── analysis_neversmoker.R        
├── analysis_reanalysis_A.R       
├── analysis_reanalysis_B.R       
├── analysis_se_mn.R              
├── analysis_supp.R               
├── bkmr_local.R                  
├── bkmr_post.R                   
├── bkmr_run.R                    
├── build_review_package.py       
├── cleaning.R                    # earlier version of the cleaning script, kept for provenance
├── cleaning_v6.R                 # FINAL cleaning script: merges the raw files, derives variables, builds data/analysis_df.rds
├── config.R                      # single source of truth for the four directories (reference implementation; see section 5)
├── config.py                     # Python mirror of config.R
├── data_dictionary.csv           
├── download_core.py              
├── download_cov.sh               
├── download_covariates.py        
├── download_data.sh              # downloads every NHANES file + the NCHS 2019 Linked Mortality File into data_raw/
├── download_missing.sh           
├── explore_vars.R                
├── fill_author.py                
├── fill_author2.py               
├── fix_manuscript.py             
├── gen_supplementary.py          
├── gen_tables.py                 
├── gen_urls.py                   
├── lod_table.R                   
├── make_figs.py                  
├── md_to_docx.py                 
├── patch_units_verify.py         
├── ph_test_10cyc.R               
├── plot_figures.py               
├── publish_to_github.ps1         
├── pubmed_key.py                 
├── pubmed_verify.py              
├── qgcss_splitsample.R           
├── rcs_fix_contrast.R            
├── repair_files.py               
├── requirements.txt              
├── rev5_manuscript_edits.py      
├── rev5_sensitivity.R            
├── rev6_A.R                      
├── rev6_B.R                      
├── rev7_edits_1.py               
├── rev7_edits_2.py               
├── rev7_edits_3.py               
├── rev7_mi_ipw.R                 
├── rev7_mi_ipw_summarize.R       
├── rev8_direction_stability.R    
├── rev9_supplementary.R          
├── run_all.sh                    
├── test_bkmr.R                   
├── trim_abstract.py              
├── trim_abstract2.py             
├── verify_sizes.py               
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

`run_all.sh` runs these 14 steps with `set -euo pipefail`:

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
13  make_figs.py              -> figures/fig1..fig5, figS1 (png/tiff/eps).  THE ONLY FIGURE ENTRY POINT
14  gen_tables.py             -> tables.md  ;  gen_supplementary.py -> submission/supplementary_tables.docx
```

Each script also runs on its own (for example `Rscript rev10_review9.R`) as long as
the preceding datasets exist. No script writes a figure that another script also
writes: `make_figs.py` is the single figure entry point.

### BKMR (optional)
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
> adults.* (2026).

**Contact:** Yuxiang Wen — Department of Cardiology, The First Affiliated
Hospital of Yangtze University, Jingzhou, Hubei, China.

## Repository URL

The repository URL appears in `CITATION.cff` (`repository-code:`) and in the Code
availability statement of the manuscript. If you fork this repository, replace
`the repository URL of your fork` with your own URL; creating a
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

`rev8_direction_stability.R` recomputes the direction-assignment stability analysis
(Table S13) with the quartile cut-points defined within each mixture's own sample, so that
the full-sample coefficients correspond exactly to the directional sums of Table 4; it
supersedes the corresponding part of `rev6_A.R`.

`rcs_fix_contrast.R` recomputes the restricted cubic spline confidence bands from the
variance of the spline contrast, which is the correct variance for a dose-response curve
relative to its reference point (Figures 4 and 5).

## Citation

If you use this code, please cite the archived version:/n/n> Wen Y. *Opposing directions in blood metal mixtures and mortality in US adults.* Zenodo. https://doi.org/10.5281/zenodo.22975596 (concept DOI; resolves to the latest version).

Source repository: https://github.com/wenyuxiangnihao/nhanes-blood-metals-mortality
