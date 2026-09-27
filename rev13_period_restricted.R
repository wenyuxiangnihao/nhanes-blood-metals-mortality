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
# 周期限定敏感性分析（2003-2018 全样本周期）
#
# 动机：1999-2002 两周期血金属仅在子样本中测定（缺失率约 72%），
#       8,640 名因无可用血样被排除者中 79% 来自这两周期（Table S24）。
#       本分析把队列限定在 2003-2018（各周期缺失率 9-18%，即全样本测定），
#       用以检验主结论是否由 1999-2002 的子样本设计驱动。
#
# 口径与主分析完全一致：同一模型 COVS、同一权重、同一事件定义；
# z 分数在限定样本内重新标准化。
# ============================================================
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu = "adjust")
D <- file.path(ROOT, "data"); R <- file.path(ROOT, "results")
dir.create(R, showWarnings = FALSE)

d0 <- readRDS(file.path(D, "analysis_df.rds"))
d <- d0[d0$cycle >= 2003, ]
d$event_cancer <- as.integer(d$mortstat == 1 & !is.na(d$ucod_leading) & d$ucod_leading == 2)

MET3 <- c("LBXBPB", "LBXBCD", "LBXTHG")
LBL  <- c(LBXBPB = "Lead", LBXBCD = "Cadmium", LBXTHG = "Mercury")
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"

for (m in MET3) {
  v <- log(d[[m]])
  d[[paste0("z_", m)]] <- (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
}

cat(sprintf("[rev13] full analysis n = %d ; cycles >= 2003 n = %d\n", nrow(d0), nrow(d)))
cat(sprintf("[rev13] deaths in the restricted set: all-cause %d, cardiovascular %d, cancer %d\n",
            sum(d$event_all), sum(d$event_cvd), sum(d$event_cancer)))

des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = d)
res <- data.frame()

# ---------- 1. single-metal survey-weighted Cox (Model 3) ----------
for (m in MET3) {
  for (ev in c("event_all", "event_cvd", "event_cancer")) {
    f <- as.formula(paste0("Surv(time,", ev, ") ~ z_", m, " + ", COVS))
    fit <- svycoxph(f, design = des)
    ci <- summary(fit)$conf.int; co <- summary(fit)$coefficients
    vn <- paste0("z_", m)
    res <- rbind(res, data.frame(
      block = "single_metal", metal = LBL[[m]], outcome = ev,
      n = nrow(d), deaths = sum(d[[ev]]),
      HR = ci[vn, "exp(coef)"], lo = ci[vn, "lower .95"], hi = ci[vn, "upper .95"],
      p = co[vn, "Pr(>|z|)"]))
  }
}

# ---------- 2. quantile g-computation, overall + directional ----------
qg <- function(dat, exps, event) {
  dq <- dat
  for (m in exps) {
    br <- unique(quantile(dat[[m]], probs = seq(0, 1, 0.25), na.rm = TRUE))
    br[1] <- -Inf; br[length(br)] <- Inf
    dq[[paste0("qc_", m)]] <- as.integer(cut(dat[[m]], breaks = br, labels = FALSE, include.lowest = TRUE)) - 1
  }
  qc <- paste0("qc_", exps)
  dq$S <- rowSums(dq[, qc])
  desq <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = dq)
  fs <- svycoxph(as.formula(paste0("Surv(time,", event, ") ~ S + ", COVS)), design = desq)
  psiS <- coef(fs)["S"]; seS <- sqrt(vcov(fs)["S", "S"])
  ff <- svycoxph(as.formula(paste0("Surv(time,", event, ") ~ ", paste(qc, collapse = "+"), " + ", COVS)), design = desq)
  b <- coef(ff)[qc]; V <- vcov(ff)[qc, qc, drop = FALSE]
  # psi and its SE are taken from the individual-component model (psi = sum of the
  # component coefficients, SE from the corresponding quadratic form), exactly as in
  # analysis_reanalysis_A.R, so that psi = sum(beta+) + sum(beta-) holds identically here
  # as well; sigma_S from the single-index model is retained for reference only.
  psi <- sum(b); se <- sqrt(sum(V))
  bp <- b[b > 0]; bn <- b[b < 0]
  sump <- if (length(bp)) sum(bp) else 0
  sumn <- if (length(bn)) sum(bn) else 0
  Vpp <- if (length(bp)) sum(V[names(bp), names(bp), drop = FALSE]) else 0
  Vnn <- if (length(bn)) sum(V[names(bn), names(bn), drop = FALSE]) else 0
  sepp <- if (Vpp > 0) sqrt(Vpp) else NA
  senn <- if (Vnn > 0) sqrt(Vnn) else NA
  list(psi = psi, se = se, psiS = psiS, seS = seS, b = b, sump = sump, sepp = sepp, sumn = sumn, senn = senn)
}

for (ev in c("event_all", "event_cvd")) {
  r <- qg(d, MET3, ev)
  res <- rbind(res, data.frame(
    block = "qgcomp_overall", metal = "3-metal mixture", outcome = ev,
    n = nrow(d), deaths = sum(d[[ev]]),
    HR = exp(r$psi), lo = exp(r$psi - 1.96 * r$se), hi = exp(r$psi + 1.96 * r$se),
    p = 2 * pnorm(-abs(r$psi / r$se))))
  res <- rbind(res,
    data.frame(block = "qgcomp_positive", metal = "3-metal mixture (positive direction)", outcome = ev,
               n = nrow(d), deaths = sum(d[[ev]]),
               HR = exp(r$sump), lo = exp(r$sump - 1.96 * r$sepp), hi = exp(r$sump + 1.96 * r$sepp), p = NA),
    data.frame(block = "qgcomp_negative", metal = "3-metal mixture (negative direction)", outcome = ev,
               n = nrow(d), deaths = sum(d[[ev]]),
               HR = exp(r$sumn), lo = exp(r$sumn - 1.96 * r$senn), hi = exp(r$sumn + 1.96 * r$senn), p = NA))
  cat(sprintf("[rev13] %s: overall %.2f (%.2f-%.2f); pos %.2f (%.2f-%.2f); neg %.2f (%.2f-%.2f)\n", ev,
              exp(r$psi), exp(r$psi - 1.96 * r$se), exp(r$psi + 1.96 * r$se),
              exp(r$sump), exp(r$sump - 1.96 * r$sepp), exp(r$sump + 1.96 * r$sepp),
              exp(r$sumn), exp(r$sumn - 1.96 * r$senn), exp(r$sumn + 1.96 * r$senn)))
  for (m in names(r$b)) cat(sprintf("      weight %-8s beta = %+.4f\n", m, r$b[[m]]))
}

write.csv(res, file.path(R, "rev13_period_restricted.csv"), row.names = FALSE)
cat("\n== period-restricted (2003-2018) sensitivity analysis ==\n"); print(res, row.names = FALSE)
