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
# REV10 -- 第九轮审稿意见要求的计算
#   (1) E-value: 两种 HR->RR 转换 (VanderWeele 近似 / Zhang-Yu with cohort p0)
#   (2) Table S14 补一列: 未加权 cause-specific Cox (分离"加权"与"尺度")
# 输出: results/rev10_evalue.csv, results/rev10_causespecific_unweighted.csv
# ============================================================
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu = "adjust")
D <- file.path(ROOT, "data"); R <- file.path(ROOT, "results")
d <- readRDS(file.path(D, "analysis_df_v6.rds"))
MET3 <- c("LBXBPB", "LBXBCD", "LBXTHG")
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
CONST <- sapply(MET3, function(m) { v <- log(d[[m]]); c(mean = mean(v), sd = sd(v)) })
for (m in MET3) d[[paste0("z_", m)]] <- (log(d[[m]]) - CONST["mean", m]) / CONST["sd", m]

cat("=== REV10 ===\n")
# ---- (1) E-value, two conversions ----
rr_vanderweele <- function(hr) (1 - 0.5^sqrt(hr)) / (1 - 0.5^sqrt(1 / hr))   # VanderWeele 2017 approximation
rr_zhangyu     <- function(hr, p0) (1 - (1 - p0)^hr) / p0                     # Zhang & Yu 1998
evalue <- function(rr) { if (rr < 1) rr <- 1 / rr; rr + sqrt(rr * (rr - 1)) }
n_all <- nrow(d); ev_all <- sum(d$event_all); ev_cvd <- sum(d$event_cvd)
p0_all <- ev_all / n_all; p0_cvd <- ev_cvd / n_all
cat(sprintf("cumulative incidence: all-cause %.4f, cardiovascular %.4f\n", p0_all, p0_cvd))
out <- data.frame()
for (s in list(c("All-cause", 0.86, 0.82, 0.89, p0_all), c("Cardiovascular", 0.90, 0.83, 0.98, p0_cvd))) {
  lab <- s[1]; hr <- as.numeric(s[2]); lo <- as.numeric(s[3]); hi <- as.numeric(s[4]); p0 <- as.numeric(s[5])
  lim <- if (hr < 1) hi else lo
  out <- rbind(out, data.frame(
    outcome = lab, HR = hr, HR_CI = paste0(lo, "\u2013", hi), p0 = round(p0, 4),
    RR_VanderWeele = round(rr_vanderweele(hr), 4), E_VanderWeele = round(evalue(rr_vanderweele(hr)), 3),
    E_VanderWeele_CI = round(evalue(rr_vanderweele(lim)), 3),
    RR_ZhangYu = round(rr_zhangyu(hr, p0), 4), E_ZhangYu = round(evalue(rr_zhangyu(hr, p0)), 3),
    E_ZhangYu_CI = round(evalue(rr_zhangyu(lim, p0)), 3),
    E_no_conversion = round(evalue(hr), 3), E_no_conversion_CI = round(evalue(lim), 3)))
}
print(out, row.names = FALSE)
write.csv(out, file.path(R, "rev10_evalue.csv"), row.names = FALSE)

# ---- (2) unweighted cause-specific Cox for the purposes of Table S14 ----
cat("\n[2] unweighted cause-specific Cox (Model 3 covariates)\n")
d$event_cancer   <- as.integer(d$mortstat == 1 & !is.na(d$ucod_leading) & d$ucod_leading == 2)
d$event_cerebro  <- as.integer(d$mortstat == 1 & !is.na(d$ucod_leading) & d$ucod_leading == 5)
d$event_clrd     <- as.integer(d$mortstat == 1 & !is.na(d$ucod_leading) & d$ucod_leading == 3)
d$event_heart    <- as.integer(d$mortstat == 1 & !is.na(d$ucod_leading) & d$ucod_leading == 1)
res <- data.frame()
for (ev in c("event_cvd", "event_heart", "event_cancer", "event_cerebro", "event_clrd")) {
  for (m in MET3) {
    f <- as.formula(paste0("Surv(time,", ev, ") ~ z_", m, " + ", COVS))
    fit <- coxph(f, data = d)                       # unweighted
    s <- summary(fit); ci <- confint(fit); term <- paste0("z_", m)
    res <- rbind(res, data.frame(outcome = ev, metal = m, n = fit$n, deaths = fit$nevent,
      HR_unweighted = round(exp(coef(fit)[term]), 3),
      lo = round(exp(ci[term, 1]), 3), hi = round(exp(ci[term, 2]), 3),
      P = signif(s$coefficients[term, ncol(s$coefficients)], 3)))
  }
}
print(res, row.names = FALSE)
write.csv(res, file.path(R, "rev10_causespecific_unweighted.csv"), row.names = FALSE)
cat("=== DONE ===\n")
