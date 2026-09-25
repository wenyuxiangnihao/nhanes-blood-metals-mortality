#!/usr/bin/env Rscript
# ============================================================
# REV7-C 汇总（可独立重跑, 不重复插补/拟合）
# 读取 results/rev7_mi_raw.rds + rev7_ipw_raw.rds, 重算 complete-case 基准,
# 重新生成 results/rev7_mi_ipw_summary.csv (修复 CC 列格式)
# ============================================================
.libPaths("/sandbox/workspace/Rlibs")
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu = "adjust")
D <- "/sandbox/workspace/heavymetal/data"; R <- "/sandbox/workspace/heavymetal/results"
MET3 <- c("LBXBPB", "LBXBCD", "LBXTHG"); OUTS <- c("event_all", "event_cvd")
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"

cc <- readRDS(file.path(D, "analysis_df_v6.rds"))
CONST <- sapply(MET3, function(m) { v <- log(cc[[m]]); c(mean = mean(v), sd = sd(v)) })
for (m in MET3) cc[[paste0("z_", m)]] <- (log(cc[[m]]) - CONST["mean", m]) / CONST["sd", m]
QCUT <- lapply(MET3, function(m) { b <- unique(quantile(cc[[m]], probs = seq(0, 1, .25), na.rm = TRUE))
  b[1] <- -Inf; b[length(b)] <- Inf; b }); names(QCUT) <- MET3

rubin <- function(bs, vs) {
  M <- length(bs); bbar <- mean(bs); W <- mean(vs); B <- var(bs)
  Tv <- W + (1 + 1 / M) * B
  df <- if (B > 0) (M - 1) * (1 + W / ((1 + 1 / M) * B))^2 else Inf
  se <- sqrt(Tv); tc <- if (is.finite(df)) qt(0.975, df) else 1.959964
  c(est = bbar, se = se, lo = bbar - tc * se, hi = bbar + tc * se, df = df, fmi = (1 + 1 / M) * B / Tv)
}
get1 <- function(fit, term) { co <- summary(fit)$coefficients; V <- vcov(fit)
  list(b = co[term, 1], v = V[term, term]) }
dirsum <- function(fit, terms) {
  co <- summary(fit)$coefficients; b <- co[terms, 1]; V <- vcov(fit)[terms, terms, drop = FALSE]
  qf <- function(S) { if (!length(S)) return(c(est = 0, v = 0))
    e <- sum(b[S]); v <- as.numeric(rep(1, length(S)) %*% V[S, S, drop = FALSE] %*% rep(1, length(S)))
    c(est = e, v = v) }
  list(pos = qf(which(b > 0)), neg = qf(which(b < 0)))
}

# --- complete-case 基准 ---
cc_b <- list(); cc_v <- list()
descc <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = cc)
for (ev in OUTS) for (m in MET3) { r <- get1(svycoxph(as.formula(paste0("Surv(time,", ev, ") ~ z_", m, " + ", COVS)), design = descc), paste0("z_", m))
  cc_b[[paste("single", ev, m, sep = "|")]] <- r$b; cc_v[[paste("single", ev, m, sep = "|")]] <- r$v }
for (ev in OUTS) {
  dd <- cc; for (m in MET3) dd[[paste0("qc_", m)]] <- as.integer(cut(dd[[m]], breaks = QCUT[[m]], labels = FALSE, include.lowest = TRUE)) - 1
  ds <- dirsum(svycoxph(as.formula(paste0("Surv(time,", ev, ") ~ qc_LBXBPB + qc_LBXBCD + qc_LBXTHG + ", COVS)),
                        design = svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = dd)),
               c("qc_LBXBPB", "qc_LBXBCD", "qc_LBXTHG"))
  cc_b[[paste("dirpos", ev, sep = "|")]] <- ds$pos["est"]; cc_v[[paste("dirpos", ev, sep = "|")]] <- ds$pos["v"]
  cc_b[[paste("dirneg", ev, sep = "|")]] <- ds$neg["est"]; cc_v[[paste("dirneg", ev, sep = "|")]] <- ds$neg["v"]
}

mi <- readRDS(file.path(R, "rev7_mi_raw.rds")); ip <- readRDS(file.path(R, "rev7_ipw_raw.rds"))
fmt <- function(b, lo, hi) sprintf("%.3f (%.3f\u2013%.3f)", b, lo, hi)

keys <- c(paste("single", rep(OUTS, each = 3), MET3, sep = "|"), paste("dirpos", OUTS, sep = "|"), paste("dirneg", OUTS, sep = "|"))
res <- data.frame()
for (k in keys) {
  kb <- sub("^single\\|", "", k)
  rmi <- rubin(mi$b[[k]], mi$v[[k]])
  rip <- rubin(ip$b[[paste("ipw", kb, sep = "|")]], ip$v[[paste("ipw", kb, sep = "|")]])
  rit <- rubin(ip$b[[paste("ipw_t", kb, sep = "|")]], ip$v[[paste("ipw_t", kb, sep = "|")]])
  b0 <- cc_b[[k]]; s0 <- sqrt(cc_v[[k]])
  # 方向分解与单金属都用 HR 尺度; 差异百分比按 log-HR 差
  dMI <- 100 * (exp(rmi["est"]) / exp(b0) - 1); dIP <- 100 * (exp(rip["est"]) / exp(b0) - 1)   # HR 尺度差异
  res <- rbind(res, data.frame(
    exposure = if (grepl("^single", k)) sub(".*\\|", "", k) else if (grepl("dirpos", k)) "Positive-direction sum" else "Negative-direction sum",
    outcome = if (grepl("event_all", k)) "All-cause" else "Cardiovascular",
    complete_case = fmt(exp(b0), exp(b0 - 1.96 * s0), exp(b0 + 1.96 * s0)),
    MI = fmt(exp(rmi["est"]), exp(rmi["lo"]), exp(rmi["hi"])), MI_fmi = round(rmi["fmi"], 3),
    IPW = fmt(exp(rip["est"]), exp(rip["lo"]), exp(rip["hi"])),
    IPW_truncated = fmt(exp(rit["est"]), exp(rit["lo"]), exp(rit["hi"])),
    pct_change_MI = round(dMI, 1), pct_change_IPW = round(dIP, 1)))
}
write.csv(res, file.path(R, "rev7_mi_ipw_summary.csv"), row.names = FALSE)
cat("=== Table S20 数据 (HR, 95% CI) ===\n"); print(res, row.names = FALSE)
cat("\nmax |%change| MI:", max(abs(res$pct_change_MI)), " IPW:", max(abs(res$pct_change_IPW)), "\n")
