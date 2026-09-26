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
# REV8 -- 重算方向稳定性，统一四分位切点口径
#
# 背景: Table 4（主分析）的 5-metal 行用 13,460 子样本内的四分位切点;
#       而 REV6-A 的 Table S13 用全部 33,104 内的切点 -> 两表五金属系数不一致。
# 本脚本统一口径: 每种混合物的切点都在"该混合物全部金属非缺失"的样本内定义
#   - 3-metal: n = 33,104
#   - 5-metal: n = 13,460
# 这样 Table S13 的 b_full 与 Table 4 的 Σβ 完全对得上。
#
# 输出: results/rev8_direction_stability.csv
# ============================================================
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu = "adjust")
D <- file.path(ROOT, "data"); R <- file.path(ROOT, "results")

d <- readRDS(file.path(D, "analysis_df_v6.rds"))
MET3 <- c("LBXBPB", "LBXBCD", "LBXTHG")
MET5 <- c(MET3, "LBXBSE", "LBXBMN")
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"

score_q <- function(dat, exps) {
  dq <- dat
  for (m in exps) {
    br <- unique(quantile(dat[[m]], probs = seq(0, 1, .25), na.rm = TRUE))
    br[1] <- -Inf; br[length(br)] <- Inf
    dq[[paste0("qc_", m)]] <- as.integer(cut(dat[[m]], breaks = br, labels = FALSE, include.lowest = TRUE)) - 1
  }
  dq
}

run_one <- function(exps, event) {
  keep <- complete.cases(d[, exps, drop = FALSE])          # 子样本内定义切点
  dd <- d[keep, , drop = FALSE]
  dq <- score_q(dd, exps)
  qc <- paste0("qc_", exps)
  f <- as.formula(paste0("Surv(time,", event, ") ~ ", paste(qc, collapse = " + "), " + ", COVS))
  des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = dq)
  full <- svycoxph(f, design = des)
  bfull <- coef(full)[qc]
  repdes <- as.svrepdesign(des, type = "JKn")
  RW <- weights(repdes, type = "analysis")
  NREP <- ncol(RW)
  out <- matrix(NA_real_, NREP, length(qc), dimnames = list(NULL, qc))
  for (r in seq_len(NREP)) {
    wr <- RW[, r]; k <- wr > 0
    fit <- try(coxph(f, data = dq[k, , drop = FALSE], weights = wr[k]), silent = TRUE)
    if (!inherits(fit, "try-error")) out[r, ] <- coef(fit)[qc]
  }
  data.frame(mixture = if (length(exps) == 3) "3-metal" else "5-metal",
             metal = sub("qc_", "", qc), outcome = event, n_sub = nrow(dq), n_reps = NREP,
             b_full = as.numeric(bfull),
             prop_b_pos = sapply(seq_along(qc), function(j) mean(out[, j] > 0, na.rm = TRUE)),
             prop_sign_consistent = sapply(seq_along(qc), function(j) {
               mean(sign(out[, j]) == sign(bfull[j]), na.rm = TRUE) }),
             sd_b = apply(out, 2, sd, na.rm = TRUE),
             q025 = apply(out, 2, quantile, .025, na.rm = TRUE),
             q975 = apply(out, 2, quantile, .975, na.rm = TRUE), row.names = NULL)
}

cat("=== REV8 direction stability (unified quartile cut-points) ===\n")
res <- data.frame()
for (exps in list(MET3, MET5)) {
  for (ev in c("event_all", "event_cvd")) {
    cat("fitting", length(exps), "-metal", ev, "...\n"); flush.console()
    res <- rbind(res, run_one(exps, ev))
  }
}
write.csv(res, file.path(R, "rev8_direction_stability.csv"), row.names = FALSE)
print(res, row.names = FALSE)

# 与主分析一致性核对
cat("\n=== 与 Table 4 / qgcomp_directional 核对 (Σβ+ / Σβ- / psi) ===\n")
qgd <- read.csv(file.path(R, "qgcomp_directional_10cyc.csv"), stringsAsFactors = FALSE)
for (mix in c("3-metal", "5-metal")) for (ev in c("event_all", "event_cvd")) {
  b <- res$b_full[res$mixture == mix & res$outcome == ev]
  pos <- sum(b[b > 0]); neg <- sum(b[b < 0])
  ref <- qgd[qgd$set == mix & qgd$outcome == ev, ]
  cat(sprintf("%-8s %-10s  S13  Σb+ = %.4f (vs %.4f)   Σb- = %.4f (vs %.4f)   psi = %.4f (vs %.4f)\n",
              mix, ev, pos, ref$sum_beta_pos, neg, ref$sum_beta_neg, pos + neg,
              ref$sum_beta_pos + ref$sum_beta_neg))
}
cat("=== DONE ===\n")
