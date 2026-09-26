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
# REV9 -- 第八轮审稿意见要求的补充分析
#   (1) E-value -> superseded, see rev10_review9.R (kept here only as a pointer)
#   (2) 铅趋势检验的替代打分 (ln 中位数 / 四分位秩次)
#   (3) 极端值缩尾敏感性 (99.5th percentile truncation)
#   (4) MI 的 between/within 方差分解 (解释 FMI 为何很小)
# 输出: results/rev9_*.csv  (附加到 Table S12 / S21)
# ============================================================
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu = "adjust")
D <- file.path(ROOT, "data"); R <- file.path(ROOT, "results")
d <- readRDS(file.path(D, "analysis_df_v6.rds"))
MET3 <- c("LBXBPB", "LBXBCD", "LBXTHG")
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
CONST <- sapply(MET3, function(m) { v <- log(d[[m]]); c(mean = mean(v), sd = sd(v)) })
for (m in MET3) d[[paste0("z_", m)]] <- (log(d[[m]]) - CONST["mean", m]) / CONST["sd", m]
mkdes <- function(dd) svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = dd)

cat("=== REV9 START ===\n")

# ---------- (1) E-value ----------
# REMOVED in round 11: this block contained a 0.5-anchored approximation
# (1 - 0.5^sqrt(hr)) / (1 - 0.5^sqrt(1/hr)) that was mislabelled as the Zhang-Yu
# conversion and produced an E-value (1.46) inconsistent with the published
# Table S12b (1.57).  The E-value calculation, with the correct Zhang-Yu equation
# and the cohort's own baseline risk, now lives solely in rev10_review9.R, whose
# output results/rev10_evalue.csv is the one quoted in the manuscript.
cat("\n[1] E-value: superseded, see rev10_review9.R and results/rev10_evalue.csv\n")

# ---------- (2) Trend tests with alternative scores ----------
cat("\n[2] Trend tests for lead (all-cause), alternative scoring\n")
des <- mkdes(d)
dq <- d
for (m in MET3) {
  br <- unique(quantile(d[[m]], probs = seq(0, 1, .25), na.rm = TRUE)); br[1] <- -Inf; br[length(br)] <- Inf
  dq[[paste0("qc_", m)]] <- as.integer(cut(d[[m]], breaks = br, labels = FALSE, include.lowest = TRUE)) - 1
  med <- tapply(d[[m]], dq[[paste0("qc_", m)]], median, na.rm = TRUE)
  dq[[paste0("med_", m)]] <- med[as.character(dq[[paste0("qc_", m)]])]      # raw-scale median
  dq[[paste0("lnmed_", m)]] <- log(dq[[paste0("med_", m)]])                  # log median
  dq[[paste0("rank_", m)]] <- dq[[paste0("qc_", m)]]                         # 0-3 rank
}
cat("  lead quartile medians (ug/dL):", paste(round(tapply(d$LBXBPB, dq$qc_LBXBPB, median, na.rm = TRUE), 2), collapse = ", "), "\n")
tr <- data.frame()
for (sv in c("med", "lnmed", "rank")) {
  for (m in MET3) {
    desx <- mkdes(dq)
    f <- as.formula(paste0("Surv(time,event_all) ~ ", sv, "_", m, " + ", COVS))
    fit <- svycoxph(f, design = desx)
    co <- summary(fit)$coefficients
    V <- vcov(fit); term <- paste0(sv, "_", m)
    # 分数变量自身的 HR（每 1 单位）与 P
    tr <- rbind(tr, data.frame(metal = m, score = sv, HR_per_unit = round(exp(coef(fit)[term]), 3),
                               lo = round(exp(coef(fit)[term] - 1.96 * sqrt(V[term, term])), 3),
                               hi = round(exp(coef(fit)[term] + 1.96 * sqrt(V[term, term])), 3),
                               P = signif(co[term, ncol(co)], 3)))
  }
}
print(tr[tr$metal == "LBXBPB", ], row.names = FALSE)
write.csv(tr, file.path(R, "rev9_trend_alternative_scores.csv"), row.names = FALSE)

# ---------- (3) Winsorised sensitivity (99.5th percentile) ----------
cat("\n[3] Winsorising at the 99.5th percentile\n")
dw <- d
for (m in MET3) {
  cap <- quantile(d[[m]], 0.995, na.rm = TRUE)
  n_cap <- sum(d[[m]] > cap, na.rm = TRUE)
  dw[[m]] <- pmin(d[[m]], cap)
  dw[[paste0("z_", m)]] <- (log(dw[[m]]) - CONST["mean", m]) / CONST["sd", m]
  cat(sprintf("  %s: cap at %.2f, %d observations capped\n", m, cap, n_cap))
}
win <- data.frame()
for (evt in c("event_all", "event_cvd")) {
  desw <- mkdes(dw)
  for (m in MET3) {
    f <- as.formula(paste0("Surv(time,", evt, ") ~ z_", m, " + ", COVS))
    fit <- svycoxph(f, design = desw); V <- vcov(fit); term <- paste0("z_", m)
    win <- rbind(win, data.frame(outcome = evt, metal = m,
      HR = round(exp(coef(fit)[term]), 3),
      lo = round(exp(coef(fit)[term] - 1.96 * sqrt(V[term, term])), 3),
      hi = round(exp(coef(fit)[term] + 1.96 * sqrt(V[term, term])), 3)))
  }
}
print(win, row.names = FALSE)
write.csv(win, file.path(R, "rev9_winsorised.csv"), row.names = FALSE)

# ---------- (4) MI between/within variance ----------
cat("\n[4] MI variance decomposition (why FMI is small)\n")
mi <- readRDS(file.path(R, "rev7_mi_raw.rds"))
dec <- data.frame()
for (k in grep("^single\\|", names(mi$b), value = TRUE)) {
  bs <- mi$b[[k]]; vs <- mi$v[[k]]; M <- length(bs)
  W <- mean(vs); B <- var(bs); T <- W + (1 + 1 / M) * B
  dec <- rbind(dec, data.frame(key = k, M = M, within_W = signif(W, 4), between_B = signif(B, 6),
                               total_T = signif(T, 4), FMI = round((1 + 1 / M) * B / T, 4),
                               pct_SE_vs_CC = NA))
}
print(dec, row.names = FALSE)
write.csv(dec, file.path(R, "rev9_mi_variance_decomposition.csv"), row.names = FALSE)
cat("=== REV9 DONE ===\n")
