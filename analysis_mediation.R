#!/usr/bin/env Rscript
# =====================================================================
# analysis_mediation.R -- Route A addition (2026-10-19)
# Survey-weighted counterfactual mediation for survival outcomes:
#   exposure: blood cadmium (log, per SD)  [primary]
#             positive-direction mixture index S+ = qc_Cd + qc_Pb [secondary]
#   mediators: eGFR (CKD-EPI 2021, all cycles)
#              log-CRP (2005-2010 CRP_D/E/F; 2015-2018 HSCRP_I/J)
#              log-uACR (2005-2018 ALB_CR_D..J) [tertiary]
#   outcomes: all-cause (primary), cardiovascular (secondary)
# Method: product method (Valeri & VanderWeele 2013) on the log-hazard
#   scale with exposure-mediator interaction; mediator model linear.
# Inference: stratified delete-1-PSU jackknife (survey.lonely.psu="adjust").
# NOTE: mediation covariate set intentionally EXCLUDES dm/htn/cvd
#   (candidate mediators themselves); see manuscript text.
# Output: results/mediation_estimates.csv, results/mediation_n.csv,
#         data/mediation_df.rds
# =====================================================================
if (!nzchar(Sys.getenv("NHANES_ROOT"))) {
  .ca <- commandArgs(trailingOnly = FALSE)
  .cf <- sub("^--file=", "", .ca[grep("^--file=", .ca)])
  Sys.setenv(NHANES_ROOT = if (length(.cf)) dirname(normalizePath(.cf)) else getwd())
}
ROOT <- Sys.getenv("NHANES_ROOT")
suppressMessages({library(survey); library(survival); library(foreign); library(parallel)})
options(survey.lonely.psu = "adjust")
RAW <- file.path(ROOT, "data_raw"); RES <- file.path(ROOT, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
d <- readRDS(file.path(ROOT, "data", "analysis_df.rds"))

rd <- function(fn){ p <- file.path(RAW, paste0(fn, ".XPT"))
  if (!file.exists(p)) return(NULL); tryCatch(read.xport(p), error = function(e) NULL) }

# ---- merge CRP / uACR ----
crp <- NULL
for (f in c("CRP_D","CRP_E","CRP_F")) { x <- rd(f)
  if (!is.null(x) && "LBXCRP" %in% names(x)) crp <- rbind(crp, data.frame(SEQN=x$SEQN, crp=x$LBXCRP)) }
for (f in c("HSCRP_I","HSCRP_J")) { x <- rd(f)
  if (!is.null(x) && "LBXHSCRP" %in% names(x)) crp <- rbind(crp, data.frame(SEQN=x$SEQN, crp=x$LBXHSCRP)) }
alb <- NULL
for (x in c("D","E","F","G","H","I","J")) { y <- rd(paste0("ALB_CR_",x))
  if (!is.null(y) && all(c("URXUMA","URXUCR") %in% names(y))) alb <- rbind(alb, data.frame(SEQN=y$SEQN, alb=y$URXUMA, ucr=y$URXUCR)) }
d <- merge(d, crp, by="SEQN", all.x=TRUE)
d <- merge(d, alb, by="SEQN", all.x=TRUE)
d$crp[d$crp <= 0] <- NA;  d$lncrp <- log(d$crp)
d$uacr <- d$alb / d$ucr * 10; d$uacr[d$uacr <= 0] <- NA; d$lnuacr <- log(d$uacr)

# ---- exposures ----
d$lncd <- log(d$LBXBCD); d$lnpb <- log(d$LBXBPB); d$lnhg <- log(d$LBXTHG)
d$cd_s <- as.numeric(scale(d$lncd))
qs <- function(v){ b <- unique(quantile(v, probs=seq(0,1,length.out=5), na.rm=TRUE))
  b[1] <- -Inf; b[length(b)] <- Inf
  as.integer(cut(v, breaks=b, include.lowest=TRUE, labels=FALSE)) - 1 }
d$qc_cd <- qs(d$LBXBCD); d$qc_pb <- qs(d$LBXBPB)
d$Splus <- d$qc_cd + d$qc_pb

COV <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke)"

# ---- point estimator for one (data, exposure, mediator, outcome) ----
est_one <- function(dat, xv, mv, ev){
  fM <- as.formula(paste0(mv, " ~ ", xv, " + ", COV))
  fY <- as.formula(paste0("Surv(time,", ev, ") ~ ", xv, " + ", mv, " + ", xv, ":", mv, " + ", COV))
  fT <- as.formula(paste0("Surv(time,", ev, ") ~ ", xv, " + ", COV))
  des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=dat)
  mmod <- svyglm(fM, design=des)
  ymod <- svycoxph(fY, design=des)
  tmod <- svycoxph(fT, design=des)
  b1m <- unname(coef(mmod)[xv])
  b0m <- mean(fitted(mmod)) - b1m * mean(dat[[xv]])   # E[M | X=0]
  th <- coef(ymod)
  th1 <- unname(th[xv]); th2 <- unname(th[mv])
  th3n <- paste0(xv, ":", mv); th3 <- if (th3n %in% names(th)) unname(th[th3n]) else 0
  TE  <- unname(coef(tmod)[xv])
  NDE <- th1 + th3 * b0m                  # contrast X: 0 -> 1
  NIE <- (th2 + th3) * (b1m)              # (theta2+theta3*a1) * delta, a1=1
  c(TE=TE, NDE=NDE, NIE=NIE, beta_M=b1m, theta1=th1, theta2=th2, theta3=th3)
}

run_combo <- function(xv, mv, ev, tag){
  need <- c("time", ev, xv, mv, "RIDAGEYR","RIAGENDR","RIDRETH1","educ","INDFMPIR","BMXBMI","smoke","SDMVPSU","SDMVSTRA","wt")
  dat <- d[complete.cases(d[, need]), ]
  n <- nrow(dat); nd <- sum(dat[[ev]]); npsu <- length(unique(interaction(dat$SDMVSTRA, dat$SDMVPSU)))
  full <- est_one(dat, xv, mv, ev)
  # jackknife: delete one PSU at a time
  dat$psu <- interaction(dat$SDMVSTRA, dat$SDMVPSU, drop=TRUE)
  dat$str <- factor(dat$SDMVSTRA)
  psus <- split(levels(dat$psu), dat$str[match(levels(dat$psu), as.character(dat$psu))])
  # safer: map each PSU to its stratum
  strmap <- tapply(as.character(dat$str), as.character(dat$psu), function(z) z[1])
  allpsu <- names(strmap)
  reps <- lapply(allpsu, function(p){
    dd <- dat[dat$psu != p, , drop=FALSE]
    out <- tryCatch(est_one(dd, xv, mv, ev), error=function(e) NULL)
    if (is.null(out)) return(NULL)
    c(psu=p, str=strmap[[p]], out)
  })
  reps <- do.call(rbind, Filter(Negate(is.null), reps))
  reps <- as.data.frame(reps, stringsAsFactors=FALSE)
  for (cc in c("TE","NDE","NIE")) reps[[cc]] <- as.numeric(reps[[cc]])
  jk_se <- function(v, s){
    v <- as.numeric(v); s <- as.character(s)
    V <- 0
    for (h in unique(s)) { vh <- v[s==h]; mh <- length(vh)
      if (mh > 1) V <- V + (mh-1)/mh * sum((vh - mean(vh))^2) }
    sqrt(V)
  }
  out <- data.frame(tag=tag, exposure=xv, mediator=mv, outcome=ev, n=n, n_deaths=nd, n_psu=npsu,
    n_jk_fail = length(allpsu) - nrow(reps),
    term=c("TE","NDE","NIE","beta_M","PM"),
    est=c(full["TE"], full["NDE"], full["NIE"], full["beta_M"], NA),
    jk_se=NA, lo=NA, hi=NA)
  pm_full <- unname(full["NIE"] / full["TE"])
  for (i in 1:4) {
    se <- jk_se(reps[[out$term[i]]], reps$str)
    out$jk_se[i] <- se; out$lo[i] <- out$est[i] - 1.96*se; out$hi[i] <- out$est[i] + 1.96*se
  }
  # PM via jackknife of NIE/TE per replicate
  pmr <- reps$NIE / reps$TE
  se_pm <- jk_se(pmr, reps$str)
  out$est[5] <- pm_full; out$jk_se[5] <- se_pm
  out$lo[5] <- pm_full - 1.96*se_pm; out$hi[5] <- pm_full + 1.96*se_pm
  out
}

combos <- list(
  c("cd_s","egfr","event_all"), c("cd_s","lncrp","event_all"),
  c("cd_s","egfr","event_cvd"), c("cd_s","lncrp","event_cvd"),
  c("Splus","egfr","event_all"), c("Splus","lncrp","event_all"),
  c("cd_s","lnuacr","event_all"))

sel <- Sys.getenv("MEDIATION_COMBOS", "")
if (nzchar(sel)) combos <- combos[as.integer(strsplit(sel, ",")[[1]])]
part_tag <- if (nzchar(sel)) paste0("_part", gsub(",", "_", sel)) else ""

res <- do.call(rbind, lapply(combos, function(cb) run_combo(cb[1], cb[2], cb[3],
  paste(cb[1], cb[2], cb[3], sep="_"))))
write.csv(res, file.path(RES, paste0("mediation_estimates", part_tag, ".csv")), row.names=FALSE)

# n table by mediator availability
ntab <- data.frame(
  egfr  = sum(complete.cases(d[, c("egfr","time","event_all","cd_s","RIDAGEYR","RIAGENDR","RIDRETH1","educ","INDFMPIR","BMXBMI","smoke")])),
  crp   = sum(complete.cases(d[, c("lncrp","time","event_all","cd_s","RIDAGEYR","RIAGENDR","RIDRETH1","educ","INDFMPIR","BMXBMI","smoke")])),
  uacr  = sum(complete.cases(d[, c("lnuacr","time","event_all","cd_s","RIDAGEYR","RIAGENDR","RIDRETH1","educ","INDFMPIR","BMXBMI","smoke")])))
write.csv(ntab, file.path(RES, "mediation_n.csv"), row.names=FALSE)
saveRDS(d, file.path(ROOT, "data", "mediation_df.rds"))
cat("== mediation done ==\n"); print(res[, c("tag","term","est","jk_se","lo","hi")], row.names=FALSE)
