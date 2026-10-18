#!/usr/bin/env Rscript
# =====================================================================
# analysis_simulation.R -- Route A addition (2026-10-19)
# Simulation: net-zero trap, direction decomposition, and survey-design
# handling in quantile g-computation for survival outcomes.
#
# DGP calibrated to the analytic sample (covariate rows resampled from
# analysis_df; log-metals = covariate-lm means + MVN residuals with the
# empirical residual covariance). Cox PH survival; baseline recalibrated
# per replicate so marginal 240-month mortality is ~13% in every scenario.
# Informative Poisson sampling (NHANES-like oversampling by age, race,
# sex), 30 strata x 2 PSUs, weights = 1/pi.
# Scenarios (per-quartile-score log-HR, read from the real 3-metal fit):
#   S1 homogeneous positive | S2 opposing (calibrated) | S3 exact net-zero
#   | S4 all-null. Estimators: unweighted / naive case-weight / full survey.
# Decomposition by estimated signs (paper) and by true signs (a-priori).
# 400 replicates per cell, seed 20261019.
# Output: results/sim_perrep.csv, results/sim_summary.csv,
#         results/sim_calibration.csv
# =====================================================================
if (!nzchar(Sys.getenv("NHANES_ROOT"))) {
  .ca <- commandArgs(trailingOnly = FALSE)
  .cf <- sub("^--file=", "", .ca[grep("^--file=", .ca)])
  Sys.setenv(NHANES_ROOT = if (length(.cf)) dirname(normalizePath(.cf)) else getwd())
}
ROOT <- Sys.getenv("NHANES_ROOT")
suppressMessages({library(survey); library(survival); library(MASS); library(parallel)})
options(survey.lonely.psu = "adjust")
RES <- file.path(ROOT, "results"); dir.create(RES, showWarnings = FALSE, recursive = TRUE)
d <- readRDS(file.path(ROOT, "data", "analysis_df.rds"))

# ---------------- calibration ----------------
d$lnPb <- log(d$LBXBPB); d$lnCd <- log(d$LBXBCD); d$lnHg <- log(d$LBXTHG)
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
lmPb <- lm(as.formula(paste("lnPb ~", COVS)), data = d)
lmCd <- lm(as.formula(paste("lnCd ~", COVS)), data = d)
lmHg <- lm(as.formula(paste("lnHg ~", COVS)), data = d)
Sres <- cov(cbind(resid(lmPb), resid(lmCd), resid(lmHg)))
Cre  <- cov2cor(Sres)

w <- read.csv(file.path(RES, "qgcomp_weights.csv"))
bw <- w[w$set == "3-metal" & w$outcome == "event_all", c("exposure","beta_unscaled")]
b_real <- setNames(bw$beta_unscaled, sub("^qc_", "", bw$exposure))   # LBXBPB LBXBCD LBXTHG

mkcols <- function(dat){
  data.frame(age10=(dat$RIDAGEYR-50)/10, male=as.integer(dat$RIAGENDR==1),
    raceMex=as.integer(dat$RIDRETH1==1), raceOH=as.integer(dat$RIDRETH1==2),
    raceNHB=as.integer(dat$RIDRETH1==4), raceOth=as.integer(dat$RIDRETH1==5),
    educ2=as.integer(dat$educ==2), educ3=as.integer(dat$educ==3),
    educ4=as.integer(dat$educ==4), educ5=as.integer(dat$educ==5),
    pir5=dat$INDFMPIR/5, bmi5=dat$BMXBMI/5,
    smokeF=as.integer(dat$smoke==1), smokeC=as.integer(dat$smoke==2),
    dm=dat$dm, htn=dat$htn, cvd=dat$cvd)
}
GAM <- c(age10=0.65, male=0.35, raceMex=-0.05, raceOH=-0.05, raceNHB=0.15, raceOth=0.00,
         educ2=0.05, educ3=-0.05, educ4=-0.15, educ5=-0.25, pir5=-0.08, bmi5=-0.02,
         smokeF=0.20, smokeC=0.60, dm=0.65, htn=0.45, cvd=0.85)
CVCOLS <- names(GAM)
BETA_M <- rbind(coef(lmPb)[-1], coef(lmCd)[-1], coef(lmHg)[-1])   # 3 x p
INT_M  <- c(coef(lmPb)[1], coef(lmCd)[1], coef(lmHg)[1])

SCEN <- list(
  S1_homog_pos = c(LBXBPB=b_real[["LBXBPB"]], LBXBCD=b_real[["LBXBCD"]], LBXTHG= abs(b_real[["LBXTHG"]])),
  S2_opposing  = c(LBXBPB=b_real[["LBXBPB"]], LBXBCD=b_real[["LBXBCD"]], LBXTHG= b_real[["LBXTHG"]]),
  S3_netzero   = c(LBXBPB=b_real[["LBXBPB"]], LBXBCD=b_real[["LBXBCD"]], LBXTHG=-(b_real[["LBXBPB"]]+b_real[["LBXBCD"]])),
  S4_null      = c(LBXBPB=0, LBXBCD=0, LBXTHG=0))
H0_TARGET <- -log(1 - 0.13) / 240
qs <- function(v){ b <- unique(quantile(v, probs=seq(0,1,length.out=5), na.rm=TRUE))
  b[1] <- -Inf; b[length(b)] <- Inf
  as.integer(cut(v, breaks=b, include.lowest=TRUE, labels=FALSE)) - 1 }

decompose <- function(b, V, plusT, minusT){
  ip <- which(b > 0); im <- which(b < 0)
  pe <- if (length(ip)) sum(b[ip]) else 0
  pe_se <- if (length(ip)) sqrt(sum(V[ip, ip])) else 0
  me <- if (length(im)) sum(b[im]) else 0
  me_se <- if (length(im)) sqrt(sum(V[im, im])) else 0
  ap <- if (any(plusT)) sum(b[plusT]) else 0
  ap_se <- if (any(plusT)) sqrt(sum(V[plusT, plusT, drop=FALSE])) else 0
  am <- if (any(minusT)) sum(b[minusT]) else 0
  am_se <- if (any(minusT)) sqrt(sum(V[minusT, minusT, drop=FALSE])) else 0
  c(pe, pe_se, me, me_se, ap, am)
}

one_rep <- function(scen_name, est, seed){
  set.seed(seed)
  bT <- SCEN[[scen_name]]
  N <- 150000
  idx <- sample.int(nrow(d), N, replace = TRUE)
  Xc <- mkcols(d[idx, ]); Xc$dm[is.na(Xc$dm)] <- 0; Xc$htn[is.na(Xc$htn)] <- 0
  eta_mu <- matrix(INT_M, nrow=N, ncol=3, byrow=TRUE) + as.matrix(Xc[, CVCOLS]) %*% t(BETA_M)
  resd <- mvrnorm(N, mu = c(0,0,0), Sigma = Sres, empirical = TRUE)
  LM <- eta_mu + resd
  qc_pop <- cbind(qs(LM[,1]), qs(LM[,2]), qs(LM[,3]))
  age <- d$RIDAGEYR[idx]
  linf <- 0.040 * (age - 50) + 0.8 * Xc$raceNHB + 0.6 * Xc$raceMex + 0.3 * Xc$male
  pincl <- exp(linf - max(linf)); pincl <- pincl / sum(pincl)
  p_i <- pmin(1, pincl * (33104 / sum(pincl)))
  samp <- runif(N) < p_i
  if (sum(samp) < 5000) return(NULL)
  sX <- Xc[samp, , drop=FALSE]; sLM <- LM[samp, , drop=FALSE]; sqc_pop <- qc_pop[samp, , drop=FALSE]
  n <- nrow(sX)
  wts <- (1 / p_i[samp]); wts <- wts * (N / sum(wts))
  strata <- cut(rank(p_i[samp], ties.method = "first"),
                breaks = seq(0, n + 1, length.out = 31), labels = FALSE)
  psu <- ave(seq_len(n), strata, FUN = function(z) sample.int(2, length(z), replace = TRUE))
  eta_s <- as.vector(sqc_pop %*% bT) + as.matrix(sX[, CVCOLS]) %*% GAM
  h0 <- H0_TARGET / mean(exp(eta_s))
  Tt <- -log(runif(n)) / (h0 * exp(eta_s))
  Tobs <- pmin(Tt, 240); ev <- as.integer(Tt <= 240)
  qc <- cbind(qs(sLM[,1]), qs(sLM[,2]), qs(sLM[,3])); colnames(qc) <- c("qPb","qCd","qHg")
  dat <- data.frame(Tobs=Tobs, ev=ev, qc, sX, w=wts, strata=factor(strata),
                    psu=factor(paste(strata, psu, sep="_")))
  F2 <- "Surv(Tobs, ev) ~ qPb + qCd + qHg + "
  FCOV <- paste(CVCOLS, collapse = " + ")
  plusT <- bT > 0; minusT <- bT < 0
  truth <- c(sum(bT), sum(bT[plusT]), sum(bT[minusT]))
  # psi is defined as Sum beta of the separate-coefficients model (paper's
  # definition); net SE from the full sub-block of the vcov matrix.
  fit2 <- switch(est,
    unwt  = coxph(as.formula(paste0(F2, FCOV)), data = dat),
    naive = coxph(as.formula(paste0(F2, FCOV)), data = dat, weights = wts),
    { des <- svydesign(id = ~psu, strata = ~strata, weights = ~w, nest = TRUE, data = dat)
      svycoxph(as.formula(paste0(F2, FCOV)), design = des) })
  V <- vcov(fit2)[c("qPb","qCd","qHg"), c("qPb","qCd","qHg"), drop = FALSE]
  b <- coef(fit2)[c("qPb","qCd","qHg")]
  o <- c(psi = sum(b), se_psi = sqrt(sum(V)), unname(b), decompose(b, V, plusT, minusT))
  out <- as.data.frame(as.list(o))
  names(out) <- c("psi","se_psi","b_Pb","b_Cd","b_Hg","plus_est","plus_est_se","minus_est","minus_est_se","plus_apri","minus_apri")
  out$scenario <- scen_name; out$estimator <- est
  out$event_rate <- mean(ev)
  out$truth_psi <- truth[1]; out$truth_plus <- truth[2]; out$truth_minus <- truth[3]
  out
}

REPS <- as.integer(Sys.getenv("SIM_REPS", "400"))
run_cell <- function(scen, est){
  res <- vector("list", REPS); errs <- character(0)
  for (r in seq_len(REPS)) {
    x <- tryCatch(one_rep(scen, est,
      seed = 20261019 + r * 97 + match(scen, names(SCEN)) * 13 + match(est, c("unwt","naive","svy")) * 7),
      error = function(e) { errs <<- c(errs, conditionMessage(e)); NULL })
    res[[r]] <- x
  }
  ans <- do.call(rbind, Filter(Negate(is.null), res))
  if (is.null(ans)) return(paste0("ALL-REP-ERROR [", scen, "/", est, "]: ",
                                  paste(unique(errs), collapse = " | ")))
  attr(ans, "errors") <- head(unique(errs), 3)
  ans
}

cells <- expand.grid(scen = names(SCEN), est = c("unwt","naive","svy"), stringsAsFactors = FALSE)
cl <- makeCluster(min(4, max(1, detectCores() - 1)))
clusterExport(cl, c("d","Sres","Cre","b_real","SCEN","GAM","CVCOLS","mkcols","qs","H0_TARGET",
                    "one_rep","decompose","run_cell","REPS","BETA_M","INT_M","COVS","RES","cells"))
clusterEvalQ(cl, {suppressMessages({library(survey); library(survival); library(MASS)})
  options(survey.lonely.psu = "adjust"); NULL})
cat("running", nrow(cells), "cells x", REPS, "reps ...\n"); t0 <- Sys.time()
allres <- parLapply(cl, seq_len(nrow(cells)), function(i) run_cell(cells$scen[i], cells$est[i]))
stopCluster(cl)
cat("cell errors:\n"); print(lapply(allres, attr, "errors"))
perrep <- do.call(rbind, Filter(Negate(is.null), allres))
perrep <- as.data.frame(perrep)
write.csv(perrep, file.path(RES, "sim_perrep.csv"), row.names = FALSE)

ci_cov <- function(est, se, truth) mean((est - 1.96*se <= truth) & (truth <= est + 1.96*se))
summ <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i){
  s <- perrep[perrep$scenario == cells$scen[i] & perrep$estimator == cells$est[i], ]
  if (!nrow(s)) return(NULL)
  data.frame(scenario=cells$scen[i], estimator=cells$est[i], nrep=nrow(s),
    ev_rate=mean(s$event_rate),
    bias_psi=mean(s$psi - s$truth_psi), rmse_psi=sqrt(mean((s$psi - s$truth_psi)^2)),
    cov_psi=ci_cov(s$psi, s$se_psi, s$truth_psi),
    bias_plus=mean(s$plus_est - s$truth_plus), cov_plus=ci_cov(s$plus_est, s$plus_est_se, s$truth_plus),
    bias_minus=mean(s$minus_est - s$truth_minus), cov_minus=ci_cov(s$minus_est, s$minus_est_se, s$truth_minus),
    bias_plus_apri=mean(s$plus_apri - s$truth_plus), cov_plus_apri=ci_cov(s$plus_apri, s$plus_est_se, s$truth_plus),
    bias_minus_apri=mean(s$minus_apri - s$truth_minus), cov_minus_apri=ci_cov(s$minus_apri, s$minus_est_se, s$truth_minus),
    falseNull=mean(abs(s$psi/s$se_psi) < 1.96 & (s$truth_plus != 0 | s$truth_minus != 0)),
    typeI_psi=mean(abs(s$psi/s$se_psi) >= 1.96),
    signrec_plus=mean(s$plus_est > 0), signrec_minus=mean(s$minus_est < 0),
    signmatch=mean(sign(s$b_Pb)==sign(SCEN[[cells$scen[i]]]["LBXBPB"]) &
                   sign(s$b_Cd)==sign(SCEN[[cells$scen[i]]]["LBXBCD"]) &
                   sign(s$b_Hg)==sign(SCEN[[cells$scen[i]]]["LBXTHG"])))
}))
write.csv(summ, file.path(RES, "sim_summary.csv"), row.names = FALSE)
cat("done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
print(summ, row.names = FALSE, digits = 3)
