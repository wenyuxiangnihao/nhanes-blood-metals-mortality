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
# REV11 -- 豆包第三轮 N1/N3 补充分析
#   (1) 三核心金属两两乘积项（加权 Cox M3，全因+CVD），验证方向分解的可加性假设
#   (2) 汞–癌症死亡的迟入组（0/24/48/60 月），闭合癌症终点 reverse causation 质疑
# 输出: results/rev11_pairwise_interaction.csv, results/rev11_cancer_landmark.csv
# ============================================================
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu = "adjust")
D <- file.path(ROOT, "data"); R <- file.path(ROOT, "results")
d <- readRDS(file.path(D, "analysis_df_plus.RDS"))
d$event_cancer <- as.integer(d$event_all == 1 & !is.na(d$ucod_leading) & d$ucod_leading == 2)
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
mkdes <- function(dd) svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = dd)
getrow <- function(fit, term) {
  s <- summary(fit)$coefficients
  i <- grep(term, rownames(s))[1]
  est <- s[i, 1]; se <- s[i, grep("robust", colnames(s))[1]]; if (is.na(se)) se <- s[i, "se(coef)"]
  HR <- exp(est); lo <- exp(est - 1.96 * se); hi <- exp(est + 1.96 * se)
  p <- 2 * pnorm(-abs(est / se))
  c(HR = HR, lo = lo, hi = hi, p = p)
}
safe <- function(expr, tag) { r <- try(force(expr), silent = TRUE); if (inherits(r, "try-error")) { cat("  [warn]", tag, ":", as.character(attr(r, "condition")), "\n"); NULL } else r }

# ---------------------------------------------------------------- (1) pairwise
cat("==== (1) pairwise metal product terms (Model 3) ====\n")
pairs <- list(c("z_LBXBCD", "z_LBXBPB"),
              c("z_LBXBCD", "z_LBXTHG"),
              c("z_LBXBPB", "z_LBXTHG"))
out <- data.frame()
for (ev in c("event_all", "event_cvd")) {
  for (pr in pairs) {
    d$.ix <- d[[pr[1]]] * d[[pr[2]]]
    fit <- safe(svycoxph(as.formula(paste0("Surv(time,", ev, ") ~ ", pr[1], " + ", pr[2], " + .ix + ", COVS)),
                         design = mkdes(d)), paste(ev, pr[1], pr[2]))
    if (is.null(fit)) next
    r <- getrow(fit, ".ix")
    out <- rbind(out, data.frame(outcome = ev, product = paste(pr[1], "x", pr[2]),
                                 n = nrow(d), n_deaths = sum(d[[ev]]),
                                 HR = r["HR"], lo = r["lo"], hi = r["hi"], p = r["p"]))
  }
}
# BH adjust within outcome
out$p_bh <- NA
for (ev in unique(out$outcome)) {
  idx <- out$outcome == ev
  out$p_bh[idx] <- p.adjust(out$p[idx], method = "BH")
}
out <- out[order(out$outcome, out$p), ]
write.csv(out, file.path(R, "rev11_pairwise_interaction.csv"), row.names = FALSE)
print(out, row.names = FALSE)

# ------------------------------------------------- (2) mercury-cancer landmark
cat("\n==== (2) mercury-cancer delayed entry 0/24/48/60 months ====\n")
s2 <- data.frame()
for (st in c(0, 24, 48, 60)) {
  dd <- d[d$time > st, ]; dd$start <- st
  fit <- safe(svycoxph(as.formula(paste0("Surv(start,time,event_cancer) ~ z_LBXTHG + ", COVS)),
                       design = mkdes(dd)), paste("cancer landmark", st))
  if (is.null(fit)) next
  r <- getrow(fit, "z_LBXTHG")
  s2 <- rbind(s2, data.frame(entry_month = st, n = nrow(dd), cancer_deaths = sum(dd$event_cancer),
                             HR = r["HR"], lo = r["lo"], hi = r["hi"], p = r["p"]))
}
write.csv(s2, file.path(R, "rev11_cancer_landmark.csv"), row.names = FALSE)
print(s2, row.names = FALSE)
cat("\nDONE rev11_addenda\n")
