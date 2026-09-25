#!/usr/bin/env Rscript
# ============================================================
# 混合暴露 (quantile g-computation, 加权Cox) + 四分位剂量反应
# ============================================================
.libPaths("/sandbox/workspace/Rlibs")
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu = "adjust")
D <- "/sandbox/workspace/heavymetal/data"; R <- "/sandbox/workspace/heavymetal/results"
d <- readRDS(file.path(D,"analysis_df.rds"))
MET3 <- c("LBXBPB","LBXBCD","LBXTHG")
MET5 <- c("LBXBPB","LBXBCD","LBXTHG","LBXBSE","LBXBMN")
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"

# ---- quantile g-computation (weighted Cox) ----
qgcomp_cox <- function(dat, expnms, event, q=4){
  dq <- dat
  for(m in expnms){
    br <- unique(quantile(dat[[m]], probs=seq(0,1,length.out=q+1), na.rm=TRUE))
    br[1] <- -Inf; br[length(br)] <- Inf
    dq[[paste0("qc_",m)]] <- as.integer(cut(dat[[m]], breaks=br, labels=FALSE, include.lowest=TRUE))-1
  }
  qcols <- paste0("qc_",expnms)
  dq$S <- rowSums(dq[,qcols])
  desq <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=dq)
  f1 <- as.formula(paste0("Surv(time,",event,") ~ S + ", COVS))
  fit1 <- svycoxph(f1, design=desq)
  psi <- coef(fit1)["S"]; se <- sqrt(vcov(fit1)["S","S"])
  f2 <- as.formula(paste0("Surv(time,",event,") ~ ", paste(qcols,collapse="+"), " + ", COVS))
  fit2 <- svycoxph(f2, design=desq); b <- coef(fit2)[qcols]
  list(psi=psi, se=se, HR=exp(psi), lo=exp(psi-1.96*se), hi=exp(psi+1.96*se), b=b)
}

out <- data.frame()
wout <- data.frame()
for(setname in c("3-metal","5-metal")){
  exps <- if(setname=="3-metal") MET3 else MET5
  dd <- if(setname=="5-metal") d[!is.na(d$LBXBSE) & !is.na(d$LBXBMN), ] else d
  for(ev in c("event_all","event_cvd")){
    r <- qgcomp_cox(dd, exps, ev)
    out <- rbind(out, data.frame(set=setname, outcome=ev, n=nrow(dd),
      HR=r$HR, lo=r$lo, hi=r$hi, p=2*pnorm(-abs(r$psi/r$se))))
    b <- r$b
    for(m in names(b)) wout <- rbind(wout, data.frame(set=setname, outcome=ev, exposure=m,
      beta_unscaled=b[[m]], weight_pos=if(b[[m]]>0) b[[m]]/sum(b[b>0]) else NA,
      weight_neg=if(b[[m]]<0) b[[m]]/sum(b[b<0]) else NA))
  }
}
write.csv(out, file.path(R,"qgcomp_mixture.csv"), row.names=FALSE)
write.csv(wout, file.path(R,"qgcomp_weights.csv"), row.names=FALSE)
cat("== quantile g-computation (per one-quartile increase in ALL metals) ==\n"); print(out, row.names=FALSE)
cat("\n== component weights ==\n"); print(wout, row.names=FALSE)

# ---- single-metal quartiles + P-trend ----
d$wt <- d$wt
res <- data.frame()
for(m in MET3){
  br <- quantile(d[[m]], probs=seq(0,1,0.25), na.rm=TRUE)
  qv <- cut(d[[m]], breaks=br, include.lowest=TRUE, labels=c("Q1","Q2","Q3","Q4"))
  d[[paste0("Q_",m)]] <- qv
  med <- tapply(d[[m]], qv, median, na.rm=TRUE)
  d[[paste0("Qmed_",m)]] <- as.numeric(med[as.character(qv)])
  for(ev in c("event_all","event_cvd")){
    desq <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=d)
    f <- as.formula(paste0("Surv(time,",ev,") ~ factor(",paste0("Q_",m),") + ", COVS))
    fit <- svycoxph(f, design=desq); s <- summary(fit)$conf.int
    for(qq in c("Q2","Q3","Q4")){
      vn <- paste0("factor(Q_",m,")",qq)
      res <- rbind(res, data.frame(metal=m, outcome=ev, level=qq, HR=s[vn,1], lo=s[vn,3], hi=s[vn,4],
        p=summary(fit)$coefficients[vn,"Pr(>|z|)"]))
    }
    ft <- as.formula(paste0("Surv(time,",ev,") ~ ",paste0("Qmed_",m)," + ", COVS))
    fitt <- svycoxph(ft, design=desq); st <- summary(fitt)$conf.int; vn2 <- paste0("Qmed_",m)
    res <- rbind(res, data.frame(metal=m, outcome=ev, level="P-trend", HR=st[vn2,1], lo=st[vn2,3], hi=st[vn2,4],
      p=summary(fitt)$coefficients[vn2,"Pr(>|z|)"]))
  }
}
write.csv(res, file.path(R,"single_metal_quartiles.csv"), row.names=FALSE)
cat("\n== quartiles ==\n"); print(res, row.names=FALSE)
