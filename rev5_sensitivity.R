#!/usr/bin/env Rscript
# ============================================================
# REV5: five requested sensitivity analyses (NHANES blood metals x mortality)
#   (1) WTSH2YR subsample-weight sensitivity
#   (2) mercury delayed-entry start 48 & 60 months (vs 24)
#   (3) Hg x Se interaction (2011-2018, n=13460)
#   (4) + factor(cycle) covariate sensitivity
#   (5) E-value for mercury inverse association
#   (6) <LOD (assay floor) exclusion sensitivity
# All outputs -> results/rev5_*.csv
# ============================================================
.libPaths("/sandbox/workspace/Rlibs")
suppressMessages({library(survey); library(survival); library(foreign)})
options(survey.lonely.psu="adjust")
D <- "/sandbox/workspace/heavymetal/data"; R <- "/sandbox/workspace/heavymetal/results"
RAW <- "/sandbox/workspace/heavymetal/data_raw"
dir.create(R, showWarnings=FALSE)

d <- readRDS(file.path(D,"analysis_df.rds"))
MET3 <- c("LBXBPB","LBXBCD","LBXTHG")
for(m in MET3){ d[[paste0("ln_",m)]] <- log(d[[m]]) }
for(m in MET3){ v<-d[[paste0("ln_",m)]]; d[[paste0("z_",m)]] <- (v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE) }

COVS_M3 <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
mkdes_w <- function(dd, wv="wt") svydesign(id=~SDMVPSU, strata=~SDMVSTRA,
                                           weights=as.formula(paste0("~",wv)), nest=TRUE, data=dd)
mkdes <- function(dd) mkdes_w(dd,"wt")
getc <- function(fit,term){ capture.output(sm <- summary(fit))   # silences design auto-print
  s<-sm$conf.int; co<-sm$coefficients
  c(HR=s[term,1], lo=s[term,3], hi=s[term,4], p=co[term,ncol(co)]) }
safe <- function(e,l){ tryCatch(e, error=function(x){cat("  !!! ERROR",l,":",conditionMessage(x),"\n");NULL}) }
WR <- function(df,name) write.csv(df, file.path(R,paste0("rev5_",name,".csv")), row.names=FALSE)

cat("analytic n =", nrow(d), " all-deaths =", sum(d$event_all), " cvd =", sum(d$event_cvd), "\n")

# =====================================================================
# (1) WTSH2YR subsample-weight sensitivity
# =====================================================================
cat("\n==== (1) WTSH2YR weight sensitivity ====\n")
h <- read.xport(file.path(RAW,"PBCD_H.XPT"))[,c("SEQN","WTSH2YR")]
i <- read.xport(file.path(RAW,"PBCD_I.XPT"))[,c("SEQN","WTSH2YR")]
wtab <- rbind(h,i)
d$WTSH2YR <- wtab$WTSH2YR[match(d$SEQN, wtab$SEQN)]
cat("WTSH2YR matched for cycle 2013:", sum(!is.na(d$WTSH2YR[d$cycle==2013])),
    "/", sum(d$cycle==2013), " 2015:", sum(!is.na(d$WTSH2YR[d$cycle==2015])), "/", sum(d$cycle==2015), "\n")

# (a) full 10-cycle sample: use WTSH2YR for 2013 & 2015, WTMEC2YR otherwise; all /10
d$wt_hyb <- d$wt
sel <- d$cycle %in% c(2013,2015)
d$wt_hyb[sel] <- d$WTSH2YR[sel]/10

s1 <- data.frame()
runw <- function(dd, wv, wlabel, ev, m){
  fit <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_",m," + ",COVS_M3)),
                       design=mkdes_w(dd,wv)), paste(wlabel,ev,m))
  if(is.null(fit)) return(NULL)
  r <- getc(fit, paste0("z_",m))
  data.frame(metal=m, outcome=ev, weight=wlabel, n=nrow(dd), n_deaths=sum(dd[[ev]]),
             HR=r["HR"], lo=r["lo"], hi=r["hi"], p=r["p"])
}
for(m in MET3) for(ev in c("event_all","event_cvd")){
  r <- runw(d, "wt",    "WTMEC2YR (main)", ev, m); if(!is.null(r)) s1 <- rbind(s1,r)
  r <- runw(d, "wt_hyb","WTSH2YR(2013-16)+WTMEC2YR(other)", ev, m); if(!is.null(r)) s1 <- rbind(s1,r)
}

# (b) only 2013-2016 subsample, WTSH2YR vs WTMEC2YR
dsub <- d[d$cycle %in% c(2013,2015), ]
cat("2013-2016 subsample n =", nrow(dsub), " all-deaths =", sum(dsub$event_all), " cvd =", sum(dsub$event_cvd), "\n")
dsub$wt_sh <- dsub$WTSH2YR/10
for(m in MET3) for(ev in c("event_all","event_cvd")){
  r <- runw(dsub, "wt",    "2013-2016: WTMEC2YR", ev, m); if(!is.null(r)) s1 <- rbind(s1,r)
  r <- runw(dsub, "wt_sh", "2013-2016: WTSH2YR", ev, m); if(!is.null(r)) s1 <- rbind(s1,r)
}
WR(s1,"wtsh2yr_weight_sensitivity"); print(s1, row.names=FALSE)

# =====================================================================
# (2) mercury delayed-entry start 48 & 60 months (vs 24)
# =====================================================================
cat("\n==== (2) delayed entry 24/48/60 months ====\n")
s2 <- data.frame()
for(st in c(24,48,60)){
  dd <- d[d$time>st,]; dd$start <- st
  for(m in MET3) for(ev in c("event_all","event_cvd")){
    fit <- safe(svycoxph(as.formula(paste0("Surv(start,time,",ev,") ~ z_",m," + ",COVS_M3)),
                         design=mkdes(dd)), paste("delayed",st,ev,m))
    if(is.null(fit)) next
    r <- getc(fit, paste0("z_",m))
    s2 <- rbind(s2, data.frame(metal=m, outcome=ev, start_month=st, n=nrow(dd), n_deaths=sum(dd[[ev]]),
                               HR=r["HR"], lo=r["lo"], hi=r["hi"], p=r["p"]))
  }
}
WR(s2,"mercury_delayed_entry"); print(s2, row.names=FALSE)

# =====================================================================
# (3) Hg x Se interaction (2011-2018, n=13460)
# =====================================================================
cat("\n==== (3) Hg x Se interaction ====\n")
s <- d[!is.na(d$LBXBSE), ]
s$z_Hg <- (log(s$LBXTHG)-mean(log(s$LBXTHG),na.rm=TRUE))/sd(log(s$LBXTHG),na.rm=TRUE)
s$z_Se <- (log(s$LBXBSE)-mean(log(s$LBXBSE),na.rm=TRUE))/sd(log(s$LBXBSE),na.rm=TRUE)
br <- quantile(s$LBXBSE, probs=c(0,1/3,2/3,1), na.rm=TRUE); br[1] <- -Inf; br[4] <- Inf
s$Se_tert <- cut(s$LBXBSE, breaks=br, include.lowest=TRUE, labels=c("T1","T2","T3"))
cat("Hg-Se subset n =", nrow(s), " all-deaths =", sum(s$event_all), " cvd =", sum(s$event_cvd),
    "  tertile n:", paste(table(s$Se_tert), collapse="/"), "\n")

inter <- data.frame(); strat <- data.frame()
wtest <- function(fit, pat=":"){ b<-coef(fit); V<-vcov(fit); idx<-grep(pat, names(b), fixed=TRUE)
  if(!length(idx)) return(c(df=0,p=NA))
  stat <- as.numeric(t(b[idx]) %*% solve(V[idx,idx,drop=FALSE]) %*% b[idx])
  c(df=length(idx), p=1-pchisq(stat,length(idx))) }
for(ev in c("event_all","event_cvd")){
  # continuous product
  fc <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_Hg*z_Se + ",COVS_M3)), design=mkdes(s)), paste("HgxSe cont",ev))
  # tertile product
  ft <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_Hg*factor(Se_tert) + ",COVS_M3)), design=mkdes(s)), paste("HgxSe tert",ev))
  # overall Hg (subset, no interaction) for reference
  fo <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_Hg + ",COVS_M3)), design=mkdes(s)), paste("Hg overall",ev))
  r0 <- getc(fo,"z_Hg")
  strat <- rbind(strat, data.frame(outcome=ev, stratum="overall(2011-2018)", n=nrow(s), n_deaths=sum(s[[ev]]),
                                   HR=r0["HR"], lo=r0["lo"], hi=r0["hi"], p=r0["p"]))
  if(!is.null(fc)){ w <- wtest(fc); rc <- getc(fc,"z_Hg:z_Se")
    inter <- rbind(inter, data.frame(outcome=ev, interaction="z_Hg x z_Se (continuous)", n=nrow(s),
      n_deaths=sum(s[[ev]]), coef_HR=rc["HR"], coef_lo=rc["lo"], coef_hi=rc["hi"], coef_p=rc["p"],
      wald_df=w["df"], wald_p=w["p"])) }
  if(!is.null(ft)){ w <- wtest(ft)
    inter <- rbind(inter, data.frame(outcome=ev, interaction="z_Hg x Se_tertile (factor)", n=nrow(s),
      n_deaths=sum(s[[ev]]), coef_HR=NA, coef_lo=NA, coef_hi=NA, coef_p=NA, wald_df=w["df"], wald_p=w["p"])) }
  # stratified by Se tertile
  for(tl in c("T1","T2","T3")){
    dd <- s[!is.na(s$Se_tert) & s$Se_tert==tl, ]
    fit <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_Hg + ",COVS_M3)), design=mkdes(dd)), paste("strat",ev,tl))
    if(is.null(fit)){ strat <- rbind(strat, data.frame(outcome=ev, stratum=paste0("Se_",tl), n=nrow(dd), n_deaths=sum(dd[[ev]]),
      HR=NA,lo=NA,hi=NA,p=NA)); next }
    rr <- getc(fit,"z_Hg")
    strat <- rbind(strat, data.frame(outcome=ev, stratum=paste0("Se_",tl), n=nrow(dd), n_deaths=sum(dd[[ev]]),
                                     HR=rr["HR"], lo=rr["lo"], hi=rr["hi"], p=rr["p"]))
  }
}
WR(inter,"hg_se_interaction"); WR(strat,"hg_se_stratified")
print(inter, row.names=FALSE); print(strat, row.names=FALSE)

# =====================================================================
# (4) + factor(cycle) covariate sensitivity
# =====================================================================
cat("\n==== (4) + factor(cycle) ====\n")
COVS_M3c <- paste0(COVS_M3, " + factor(cycle)")
s4 <- data.frame()
for(m in MET3) for(ev in c("event_all","event_cvd")){
  f0 <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_",m," + ",COVS_M3)),  design=mkdes(d)), paste("base",ev,m))
  f1 <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_",m," + ",COVS_M3c)), design=mkdes(d)), paste("cyc",ev,m))
  for(f in list(list(f0,"M3 (no cycle)"), list(f1,"M3 + factor(cycle)"))){
    fit <- f[[1]]; lb <- f[[2]]; if(is.null(fit)) next
    r <- getc(fit, paste0("z_",m))
    s4 <- rbind(s4, data.frame(metal=m, outcome=ev, model=lb, n=nrow(d), n_deaths=sum(d[[ev]]),
                               HR=r["HR"], lo=r["lo"], hi=r["hi"], p=r["p"]))
  }
}
WR(s4,"cycle_covariate"); print(s4, row.names=FALSE)

# =====================================================================
# (5) E-value for mercury inverse association
# =====================================================================
cat("\n==== (5) E-value (mercury) ====\n")
evalue <- function(hr, lo, hi){
  # treat HR as RR (rare outcome). For protective (RR<1) invert first.
  hr <- as.numeric(hr); lo <- as.numeric(lo); hi <- as.numeric(hi)
  rr <- ifelse(hr<1, 1/hr, hr)
  e_pt <- rr + sqrt(rr*(rr-1))
  cl <- ifelse(hr<1, hi, lo)            # CI limit closest to the null (1)
  rr2 <- ifelse(cl<1, 1/cl, cl)
  e_ci <- if(!is.na(cl)) rr2 + sqrt(rr2*(rr2-1)) else NA
  c(rr_point=rr, E_point=e_pt, ci_limit=cl, rr_ci=ifelse(is.na(cl),NA,rr2), E_ci=e_ci)
}
s5 <- data.frame()
for(ev in c("event_all","event_cvd")){
  fit <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_LBXTHG + ",COVS_M3)), design=mkdes(d)), paste("evalue",ev))
  r <- getc(fit,"z_LBXTHG")
  e <- evalue(r["HR"], r["lo"], r["hi"])
  s5 <- rbind(s5, data.frame(outcome=ev, HR=r["HR"], lo=r["lo"], hi=r["hi"], p=r["p"],
    RR_point=e["rr_point"], E_value_point=e["E_point"],
    CI_limit_nearest_null=e["ci_limit"], RR_ci=e["rr_ci"], E_value_CI=e["E_ci"]))
}
WR(s5,"evalue_mercury"); print(s5, row.names=FALSE)

# =====================================================================
# (6) <LOD (assay floor) exclusion sensitivity
# =====================================================================
cat("\n==== (6) <LOD / assay-floor exclusion ====\n")
lod <- read.csv(file.path(R,"lod_by_cycle_10cyc.csv"), stringsAsFactors=FALSE)
cyc_num <- c("1999-2000"=1999,"2001-2002"=2001,"2003-2004"=2003,"2005-2006"=2005,"2007-2008"=2007,
             "2009-2010"=2009,"2011-2012"=2011,"2013-2014"=2013,"2015-2016"=2015,"2017-2018"=2017)
lod$cycnum <- as.numeric(cyc_num[lod$cycle])
s6 <- data.frame()
for(m in MET3){
  lk <- lod[lod$metal==m, c("cycnum","min_obs")]
  mo <- lk$min_obs[match(d$cycle, lk$cycnum)]
  keep <- !( !is.na(d[[m]]) & !is.na(mo) & abs(d[[m]]-mo) < 1e-6 )
  dd <- d[keep, ]
  cat(sprintf("  %s: excluded at floor = %d ; remaining n = %d\n", m, sum(!keep), nrow(dd)))
  for(ev in c("event_all","event_cvd")){
    f0 <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_",m," + ",COVS_M3)), design=mkdes(d)),  paste("lod base",m,ev))
    f1 <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_",m," + ",COVS_M3)), design=mkdes(dd)), paste("lod excl",m,ev))
    r0 <- getc(f0, paste0("z_",m)); r1 <- getc(f1, paste0("z_",m))
    s6 <- rbind(s6, data.frame(metal=m, outcome=ev, model="M3 (all, incl. floor)", n=nrow(d), n_deaths=sum(d[[ev]]),
      n_excluded=0, HR=r0["HR"], lo=r0["lo"], hi=r0["hi"], p=r0["p"]))
    s6 <- rbind(s6, data.frame(metal=m, outcome=ev, model="M3 (excl. assay floor)", n=nrow(dd), n_deaths=sum(dd[[ev]]),
      n_excluded=sum(!keep), HR=r1["HR"], lo=r1["lo"], hi=r1["hi"], p=r1["p"]))
  }
}
WR(s6,"lod_floor_exclusion"); print(s6, row.names=FALSE)

cat("\n=== rev5 DONE ===\n")
