#!/usr/bin/env Rscript
# --- project root: $NHANES_ROOT, or the folder containing this script --------
if (!nzchar(Sys.getenv("NHANES_ROOT"))) {
  .ca <- commandArgs(trailingOnly = FALSE)
  .cf <- sub("^--file=", "", .ca[grep("^--file=", .ca)])
  Sys.setenv(NHANES_ROOT = if (length(.cf)) dirname(normalizePath(.cf)) else getwd())
}
ROOT <- Sys.getenv("NHANES_ROOT")
if (dir.exists(file.path(dirname(ROOT), "Rlibs")))
  .libPaths(c(file.path(dirname(ROOT), "Rlibs"), .libPaths()))
# ============================================================
# REV12 -- Table S24: selection-bias characterisation
#   Compares participants excluded because the three core blood metals were not
#   available (n = 8,640) with the analytic sample (n = 33,104) on the variables
#   used in the main models, unweighted and survey-weighted, and documents how
#   much of the covariate profile is observable in the excluded group.
# Inputs : data/all_merged_df.rds (written by cleaning_v6.R before exclusion)
#          data/analysis_df.rds    (analytic sample, n = 33,104)
# Output : results/rev12_tableS24.md, results/rev12_tableS24.csv
# ============================================================
suppressMessages({library(survey)})
options(survey.lonely.psu = "adjust")
D <- file.path(ROOT, "data"); R <- file.path(ROOT, "results")
NCYC <- 10

all <- readRDS(file.path(D, "all_merged_df.rds"))
cat("[rev12] full merged frame n =", nrow(all), "\n")

e1 <- all[!is.na(all$RIDAGEYR) & all$RIDAGEYR >= 20, ]
e2 <- e1[is.na(e1$RIDEXPRG) | e1$RIDEXPRG != 1, ]
e3 <- e2[!is.na(e2$eligstat) & e2$eligstat == 1, ]
has3 <- !(is.na(e3$LBXBPB) | is.na(e3$LBXBCD) | is.na(e3$LBXTHG))
exc <- e3[!has3, ]
inc <- readRDS(file.path(D, "analysis_df.rds"))
inc$wt <- inc$WTMEC2YR / NCYC
exc$wt <- exc$WTMEC2YR / NCYC
cat(sprintf("[rev12] mortality-eligible %d | three metals %d | excluded %d | analytic %d\n",
            nrow(e3), sum(has3), nrow(exc), nrow(inc)))

# how much of the excluded group comes from cycles whose metals are subsample-only
tab <- table(factor(exc$cycle, levels = sort(unique(exc$cycle))))
cat("[rev12] excluded by cycle:", paste(names(tab), as.integer(tab), sep = "=", collapse = " "), "\n")
early <- sum(exc$cycle %in% c(1999, 2001)); pct_early <- round(100 * early / nrow(exc))
cat(sprintf("[rev12] from 1999-2002 cycles: %d (%.0f%%)\n", early, 100 * early / nrow(exc)))

ok_w <- !is.na(exc$wt) & exc$wt > 0
ok_wi <- !is.na(inc$wt) & inc$wt > 0
cat(sprintf("[rev12] positive MEC weight: excluded %d (%.1f%%), included %d (%.1f%%)\n",
            sum(ok_w), 100 * mean(ok_w), sum(ok_wi), 100 * mean(ok_wi)))

cont <- list(
  list(v = "RIDAGEYR", lab = "Age, years", dig = 1),
  list(v = "INDFMPIR", lab = "Family income-to-poverty ratio", dig = 2),
  list(v = "BMXBMI",  lab = "Body-mass index, kg/m2", dig = 1),
  list(v = "sbp",     lab = "Systolic blood pressure, mmHg", dig = 1),
  list(v = "dbp",     lab = "Diastolic blood pressure, mmHg", dig = 1),
  list(v = "egfr",    lab = "eGFR, mL/min/1.73 m2", dig = 1),
  list(v = "tc",      lab = "Total cholesterol, mg/dL", dig = 1),
  list(v = "hdl",     lab = "HDL cholesterol, mg/dL", dig = 1),
  list(v = "cot",     lab = "Serum cotinine, ng/mL", dig = 2)
)
bin <- list(
  list(v = "female", lab = "Female"),
  list(v = "smoke2", lab = "Current smoker"),
  list(v = "dm",     lab = "Diabetes"),
  list(v = "htn",    lab = "Hypertension"),
  list(v = "cvd",    lab = "Baseline cardiovascular disease"),
  list(v = "cancer_ever", lab = "Self-reported cancer at baseline")
)
exc$female <- as.integer(exc$RIAGENDR == 2); exc$smoke2 <- as.integer(exc$smoke == 2)
inc$female <- as.integer(inc$RIAGENDR == 2); inc$smoke2 <- as.integer(inc$smoke == 2)
if (!"cancer_ever" %in% names(exc)) {
  exc$cancer_ever <- ifelse(exc$MCQ220 == 1, 1L, ifelse(exc$MCQ220 == 2, 0L, NA_integer_))
  exc$cancer_ever[exc$MCQ220 %in% c(7, 9)] <- NA
}

smd_cont <- function(x1, x0) {
  m1 <- mean(x1, na.rm = TRUE); m0 <- mean(x0, na.rm = TRUE)
  s1 <- sd(x1, na.rm = TRUE); s0 <- sd(x0, na.rm = TRUE)
  (m1 - m0) / sqrt((s1^2 + s0^2) / 2)
}
smd_bin <- function(p1, p0) (p1 - p0) / sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
fmt <- function(x, d) formatC(x, format = "f", digits = d)
pct <- function(x) formatC(100 * x, format = "f", digits = 1)

rows <- character()
addrow <- function(...) rows <<- c(rows, sprintf("| %s | %s | %s | %s |", ...))
addrow("**Unweighted comparison**", "", "", "")
for (x in cont) {
  a <- exc[[x$v]]; b <- inc[[x$v]]
  addrow(sprintf("%s, mean (SD)", x$lab),
         sprintf("%s (%s)", fmt(mean(a, na.rm = TRUE), x$dig), fmt(sd(a, na.rm = TRUE), x$dig)),
         sprintf("%s (%s)", fmt(mean(b, na.rm = TRUE), x$dig), fmt(sd(b, na.rm = TRUE), x$dig)),
         fmt(smd_cont(a, b), 2))
}
for (x in bin) {
  p1 <- mean(exc[[x$v]], na.rm = TRUE); p0 <- mean(inc[[x$v]], na.rm = TRUE)
  addrow(sprintf("%s, %%", x$lab), pct(p1), pct(p0), fmt(smd_bin(p1, p0), 2))
}
catrow <- function(var, lev, labvar, labels) {
  p1 <- mean(exc[[var]] == lev, na.rm = TRUE); p0 <- mean(inc[[var]] == lev, na.rm = TRUE)
  addrow(sprintf("%s: %s, %%", labvar, labels[as.character(lev)]), pct(p1), pct(p0), fmt(smd_bin(p1, p0), 2))
}
race_lab <- c("1" = "Mexican American", "2" = "Other Hispanic", "3" = "Non-Hispanic White",
              "4" = "Non-Hispanic Black", "5" = "Other/multiracial")
edu_lab <- c("1" = "<9th grade", "2" = "9-11th grade", "3" = "High school/GED",
             "4" = "Some college", "5" = "College graduate")
smk_lab <- c("0" = "Never", "1" = "Former", "2" = "Current")
for (l in 1:5) catrow("RIDRETH1", l, "Race/ethnicity", race_lab)
for (l in 1:5) catrow("educ", l, "Education", edu_lab)
for (l in 0:2) catrow("smoke", l, "Smoking status", smk_lab)

rows <- c(rows, "| **Data availability in the excluded group** | | | |")
for (x in c(cont, bin)) {
  v <- x$v
  rows <- c(rows, sprintf("| %s: value available, %% | %s | %s |  |",
                          x$lab, pct(mean(!is.na(exc[[v]]))), pct(mean(!is.na(inc[[v]])))))
}

# --- survey-weighted comparison (both groups carry a positive MEC weight) --------
excw <- exc[ok_w, ]; incw <- inc[ok_wi, ]
des_exc <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = excw)
des_inc <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = incw)
binvars <- sapply(bin, function(z) z$v)
mrow <- character()
for (x in c(cont, bin)) {
  v <- x$v
  if (sum(!is.na(excw[[v]])) < 30 || sum(!is.na(incw[[v]])) < 30) next
  f <- as.formula(paste0("~", v))
  m1 <- try(coef(svymean(f, des_exc, na.rm = TRUE)), silent = TRUE)
  m0 <- try(coef(svymean(f, des_inc, na.rm = TRUE)), silent = TRUE)
  if (inherits(m1, "try-error") || inherits(m0, "try-error")) next
  if (v %in% binvars) {
    mrow <- c(mrow, sprintf("| %s, %% | %s | %s |  |", x$lab, pct(as.numeric(m1)), pct(as.numeric(m0))))
  } else {
    dig <- if (v %in% c("INDFMPIR", "cot")) 2 else 1
    mrow <- c(mrow, sprintf("| %s | %s | %s |  |", x$lab, fmt(as.numeric(m1), dig), fmt(as.numeric(m0), dig)))
  }
}
wcat <- function(var, levs, labels, labvar) {
  f <- as.formula(paste0("~factor(", var, ")"))
  m1 <- try(coef(svymean(f, des_exc, na.rm = TRUE)), silent = TRUE)
  m0 <- try(coef(svymean(f, des_inc, na.rm = TRUE)), silent = TRUE)
  if (inherits(m1, "try-error") || inherits(m0, "try-error")) return(character())
  out <- character()
  for (l in levs) {
    k <- grep(paste0("^factor\\(", var, "\\)", l, "$"), names(m1))
    if (!length(k)) next
    out <- c(out, sprintf("| %s: %s, %% | %s | %s |  |", labvar, labels[as.character(l)], pct(m1[k]), pct(m0[k])))
  }
  out
}
mrow <- c(mrow, wcat("RIDRETH1", 1:5, race_lab, "Race/ethnicity"),
          wcat("educ", 1:5, edu_lab, "Education"),
          wcat("smoke", 0:2, smk_lab, "Smoking status"))

note <- paste0(
  "Note: the 8,640 participants excluded at this step had no usable blood specimen for lead, cadmium or mercury; ",
  format(early, big.mark = ","), " of them (", pct_early, "%) come from the 1999-2002 cycles, in which blood metals ",
  "were measured only in a subsample of participants, and the remainder had insufficient or unusable specimens. ",
  "The number excludes the 3,639 participants excluded later for missing covariates (Tables S18 and S20) and the ",
  "8 participants without follow-up time. All 8,640 carry a positive MEC examination weight, so both comparisons ",
  "can be computed; the unweighted comparison describes the excluded individuals as they are, and the ",
  "survey-weighted comparison shows how the two groups compare once the sampling design is accounted for. ",
  "Accounting for the design brings most variables closer (age 48.8 vs 46.4 years; income-to-poverty ratio 2.95 ",
  "vs 3.01; diabetes 8.4% vs 8.6%), but the excluded group is still more often male (34.2% vs 54.7% female), and ",
  "has higher systolic blood pressure (127.2 vs 121.7 mmHg), higher total cholesterol (206.9 vs 195.2 mg/dL) and ",
  "lower HDL cholesterol (49.5 vs 53.4 mg/dL); the two groups are therefore not fully exchangeable. Percentages ",
  "of complete data indicate how much of the covariate profile is observable in the excluded group: serum ",
  "chemistry (eGFR, cholesterol, cotinine) is available for about 72% of them, consistent with specimen shortage.")

md <- c("**Table S24. Characteristics of participants excluded for missing blood metals versus the analytic sample**",
        "",
        "| Characteristic | Excluded (no usable blood specimen), n = 8,640 | Analytic sample, n = 33,104 | SMD |",
        "|---|---|---|---|",
        paste(rows, collapse = "\n"),
        "",
        "**Survey-weighted comparison** (both groups weighted by the MEC examination weight):",
        "",
        "| Characteristic | Excluded | Analytic sample |  |",
        "|---|---|---|---|",
        paste(mrow, collapse = "\n"),
        "",
        note,
        "")

writeLines(md, file.path(R, "rev12_tableS24.md"))
write.csv(data.frame(line = md), file.path(R, "rev12_tableS24.csv"), row.names = FALSE)
cat("[rev12] wrote results/rev12_tableS24.md (", length(md), "lines )\n")
cat(paste(md, collapse = "\n"), "\n")
