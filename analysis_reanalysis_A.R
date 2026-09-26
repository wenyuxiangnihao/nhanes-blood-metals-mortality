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
# RE-ANALYSIS (10 cycles incl. 2003-2004) -- Part A: tasks 2(a)-(g)
# ============================================================
suppressMessages({library(survey); library(survival); library(foreign); library(splines)})
options(survey.lonely.psu="adjust")
D <- file.path(ROOT, "data"); R <- file.path(ROOT, "results")
RAW <- file.path(ROOT, "data_raw")
dir.create(R, showWarnings=FALSE)
d <- readRDS(file.path(D,"analysis_df.rds"))
cat("analytic n =", nrow(d), " all-deaths =", sum(d$event_all), " cvd =", sum(d$event_cvd), "\n")

MET3 <- c("LBXBPB","LBXBCD","LBXTHG"); MET5 <- c(MET3,"LBXBSE","LBXBMN")
allMET <- MET5
for(m in allMET){ d[[paste0("ln_",m)]] <- log(d[[m]]) }
for(m in allMET){ v<-d[[paste0("ln_",m)]]; d[[paste0("z_",m)]] <- (v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE) }

COVS_M2 <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1)"
COVS_M3 <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
mkdes <- function(dd) svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=dd)
WR  <- function(df,name){ write.csv(df,file.path(R,paste0(name,"_10cyc.csv")),row.names=FALSE)
                          write.csv(df,file.path(R,paste0(name,".csv")),row.names=FALSE) }
WRo <- function(df,name) write.csv(df,file.path(R,paste0(name,"_10cyc.csv")),row.names=FALSE)
getc <- function(fit,term){ s<-summary(fit)$conf.int; co<-summary(fit)$coefficients
  c(HR=s[term,1], lo=s[term,3], hi=s[term,4], p=co[term,ncol(co)]) }
safe <- function(expr,label){ tryCatch(expr, error=function(e){ cat("  !!! ERROR in",label,":",conditionMessage(e),"\n"); NULL }) }

# ---------------- (a) single-metal Cox, 3 models ----------------
cat("\n== (a) single-metal Cox M1/M2/M3 ==\n")
sa <- data.frame()
for(m in MET3){
  for(ev in c("event_all","event_cvd")){
    vn <- paste0("z_",m)
    fits <- list(
      M1 = safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ ",vn)), design=mkdes(d)),"M1"),
      M2 = safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ ",vn," + ",COVS_M2)), design=mkdes(d)),"M2"),
      M3 = safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ ",vn," + ",COVS_M3)), design=mkdes(d)),"M3"))
    for(lb in c("M1","M2","M3")){
      fit <- fits[[lb]]; if(is.null(fit)) next
      r <- getc(fit, vn)
      sa <- rbind(sa, data.frame(metal=m, outcome=ev, model=lb, n=nrow(d),
        HR=r["HR"], lo=r["lo"], hi=r["hi"], p=r["p"]))
    }
  }
}
WR(sa,"single_metal_cox_M1M2M3"); print(sa, row.names=FALSE)

# ---------------- (c) weighted SD of ln, per-doubling, per-IQR ----------------
cat("\n== (c) weighted SD / per-doubling / per-IQR ==\n")
wvar <- function(x,w){ ok<-!is.na(x)&!is.na(w); x<-x[ok]; w<-w[ok]; mu<-sum(w*x)/sum(w); sum(w*(x-mu)^2)/(sum(w)-sum(w^2)/sum(w)) }
sc <- data.frame()
for(m in MET3){
  lnm <- d[[paste0("ln_",m)]]; w <- d$wt
  wsd <- sqrt(wvar(lnm,w)); usd <- sd(lnm,na.rm=TRUE)
  iqr <- IQR(lnm, na.rm=TRUE); wiqr <- as.numeric(diff(quantile(lnm,c(.25,.75),na.rm=TRUE)))
  for(ev in c("event_all","event_cvd")){
    fit <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ ln_",m," + ",COVS_M3)), design=mkdes(d)),"lnmodel")
    if(is.null(fit)) next
    b <- coef(fit)[paste0("ln_",m)]; V <- vcov(fit); se <- sqrt(V[paste0("ln_",m),paste0("ln_",m)])
    dbl <- exp(b*log(2)); dl<-exp((b-1.96*se)*log(2)); dh<-exp((b+1.96*se)*log(2))
    iq  <- exp(b*iqr);    il<-exp((b-1.96*se)*iqr); ih<-exp((b+1.96*se)*iqr)
    sc <- rbind(sc, data.frame(metal=m, outcome=ev, wt_sd_ln=wsd, unwt_sd_ln=usd,
      iqr_ln=iqr, HR_per_doubling=dbl, doubling_lo=dl, doubling_hi=dh,
      HR_per_IQR=iq, iqr_lo=il, iqr_hi=ih))
  }
}
WR(sc,"metal_scale"); print(sc, row.names=FALSE)

# ---------------- (b) quartiles Q2-Q4 + P-trend + P-nonlinear(quadratic) ----------
cat("\n== (b) quartiles ==\n")
sb <- data.frame()
for(m in MET3){
  br <- quantile(d[[m]], probs=seq(0,1,0.25), na.rm=TRUE)
  d[[paste0("Q_",m)]] <- cut(d[[m]], breaks=br, include.lowest=TRUE, labels=c("Q1","Q2","Q3","Q4"))
  med <- tapply(d[[m]], d[[paste0("Q_",m)]], median, na.rm=TRUE)
  d[[paste0("Qmed_",m)]] <- as.numeric(med[as.character(d[[paste0("Q_",m)]])])
  for(ev in c("event_all","event_cvd")){
    desq <- mkdes(d)
    qf <- paste0("Q_",m)
    fit <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ factor(",qf,") + ",COVS_M3)), design=desq),"Qmodel")
    if(!is.null(fit)) for(qq in c("Q2","Q3","Q4")){
      vn <- paste0("factor(",qf,")",qq); r<-getc(fit,vn)
      sb <- rbind(sb, data.frame(metal=m, outcome=ev, level=qq, HR=r["HR"], lo=r["lo"], hi=r["hi"], p=r["p"]))
    }
    ft <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ Qmed_",m," + ",COVS_M3)), design=desq),"trend")
    if(!is.null(ft)){ r<-getc(ft,paste0("Qmed_",m)); sb <- rbind(sb, data.frame(metal=m,outcome=ev,level="P-trend",HR=r["HR"],lo=r["lo"],hi=r["hi"],p=r["p"])) }
    # quadratic p-nonlinear
    dq <- d; dq$z1 <- dq[[paste0("z_",m)]]; dq$z2 <- dq$z1^2
    fq <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z1 + z2 + ",COVS_M3)), design=mkdes(dq)),"quad")
    if(!is.null(fq)){ b<-coef(fq); V<-vcov(fq); st<-unname(b["z2"]^2/V["z2","z2"]); p<-1-pchisq(st,1)
      sb <- rbind(sb, data.frame(metal=m,outcome=ev,level="P-nonlinear(quad)",HR=NA,lo=NA,hi=NA,p=p)) }
  }
}
WR(sb,"single_metal_quartiles"); print(sb, row.names=FALSE)

# ---------------- (d) weighted correlation matrix of ln(metals) ----------------
cat("\n== (d) weighted correlations ==\n")
wcor <- function(x,y,w){ ok<-!is.na(x)&!is.na(y)&!is.na(w); x<-x[ok];y<-y[ok];w<-w[ok]
  mx<-sum(w*x)/sum(w); my<-sum(w*y)/sum(w)
  cxy<-sum(w*(x-mx)*(y-my))/sum(w); vx<-sum(w*(x-mx)^2)/sum(w); vy<-sum(w*(y-my)^2)/sum(w)
  r<-cxy/sqrt(vx*vy); n<-length(x); z<-atanh(r); se<-1/sqrt(max(n-3,1)); p<-2*pnorm(-abs(z/se))
  c(r=r,p=p,n=n) }
scor <- data.frame()
for(setnm in c("3-metal","5-metal")){
  ms <- if(setnm=="3-metal") MET3 else MET5
  dd <- if(setnm=="5-metal") d[!is.na(d$LBXBSE)&!is.na(d$LBXBMN),] else d
  for(i in 1:(length(ms)-1)) for(j in (i+1):length(ms)){
    r <- wcor(dd[[paste0("ln_",ms[i])]], dd[[paste0("ln_",ms[j])]], dd$wt)
    scor <- rbind(scor, data.frame(set=setnm, v1=ms[i], v2=ms[j], r=r["r"], p=r["p"], n=r["n"]))
  }
}
WRo(scor,"metal_correlations"); print(scor, row.names=FALSE)

# ---------------- (e) quantile g-computation + directional decomposition ---------
cat("\n== (e) qgcomp (manual) + directional ==\n")
qg <- function(dat, exps, event){
  dq <- dat
  for(m in exps){ br<-unique(quantile(dat[[m]],probs=seq(0,1,0.25),na.rm=TRUE)); br[1]<- -Inf; br[length(br)]<-Inf
    dq[[paste0("qc_",m)]] <- as.integer(cut(dat[[m]],breaks=br,labels=FALSE,include.lowest=TRUE))-1 }
  qc <- paste0("qc_",exps)
  dq$S <- rowSums(dq[,qc])
  desq <- svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~wt,nest=TRUE,data=dq)
  fS <- as.formula(paste0("Surv(time,",event,") ~ S + ",COVS_M3))
  fs <- svycoxph(fS,design=desq)
  psiS <- coef(fs)["S"]; seS <- sqrt(vcov(fs)["S","S"])
  fF <- as.formula(paste0("Surv(time,",event,") ~ ",paste(qc,collapse="+")," + ",COVS_M3))
  ff <- svycoxph(fF,design=desq); b <- coef(ff)[qc]; V <- vcov(ff)[qc,qc]
  psi <- sum(b); v <- sum(V)
  bp <- b[b>0]; bn <- b[b<0]
  sump <- if(length(bp)) sum(bp) else 0; sumn <- if(length(bn)) sum(bn) else 0
  Vpp <- if(length(bp)) sum(V[names(bp),names(bp),drop=FALSE]) else 0
  Vnn <- if(length(bn)) sum(V[names(bn),names(bn),drop=FALSE]) else 0
  sepp <- if(Vpp>0) sqrt(Vpp) else NA; senn <- if(Vnn>0) sqrt(Vnn) else NA
  list(b=b, psiS=psiS, seS=seS, psi=psi, se=sqrt(v), sump=sump, sepp=sepp, sumn=sumn, senn=senn, n=nrow(dq))
}
smi <- data.frame(); swt <- data.frame(); sdir <- data.frame()
for(setnm in c("3-metal","5-metal")){
  exps <- if(setnm=="3-metal") MET3 else MET5
  dd <- if(setnm=="5-metal") d[!is.na(d$LBXBSE)&!is.na(d$LBXBMN),] else d
  for(ev in c("event_all","event_cvd")){
    r <- safe(qg(dd,exps,ev),paste0("qg ",setnm,ev)); if(is.null(r)) next
    smi <- rbind(smi, data.frame(set=setnm,outcome=ev,n=r$n,
      psi_overall=r$psi, HR_overall=exp(r$psi), lo=exp(r$psi-1.96*r$se), hi=exp(r$psi+1.96*r$se),
      p=2*pnorm(-abs(r$psi/r$se))))
    sdir <- rbind(sdir, data.frame(set=setnm,outcome=ev,
      sum_beta_pos=r$sump, HR_pos=exp(r$sump), pos_lo=exp(r$sump-1.96*r$sepp), pos_hi=exp(r$sump+1.96*r$sepp),
      sum_beta_neg=r$sumn, HR_neg=exp(r$sumn), neg_lo=exp(r$sumn-1.96*r$senn), neg_hi=exp(r$sumn+1.96*r$senn)))
    for(m in names(r$b)) swt <- rbind(swt, data.frame(set=setnm,outcome=ev,exposure=m,
      beta_unscaled=r$b[[m]],
      weight_pos=if(r$b[[m]]>0) r$b[[m]]/sum(r$b[r$b>0]) else NA,
      weight_neg=if(r$b[[m]]<0) r$b[[m]]/sum(r$b[r$b<0]) else NA))
  }
}
WRo(smi,"qgcomp_overall"); WRo(swt,"qgcomp_weights"); WRo(sdir,"qgcomp_directional")
print(smi,row.names=FALSE); print(sdir,row.names=FALSE); print(swt,row.names=FALSE)

# ---------------- (f) sensitivity analyses (all-cause) ----------------
cat("\n== (f) sensitivities ==\n")
# alcohol + physical activity covariates
alc <- d$ALQ130; alc[!is.na(d$ALQ101) & d$ALQ101==2] <- 0
d$alc <- alc
paqmap <- c("1999"="PAQ","2001"="PAQ_B","2003"="PAQ_C","2005"="PAQ_D","2007"="PAQ_E","2009"="PAQ_F",
            "2011"="PAQ_G","2013"="PAQ_H","2015"="PAQ_I","2017"="PAQ_J")
paq <- data.frame()
for(y in names(paqmap)){
  fp <- file.path(RAW,paste0(paqmap[y],".XPT")); if(!file.exists(fp)) next
  x <- tryCatch(read.xport(fp),error=function(e)NULL); if(is.null(x)) next
  if(all(c("PAQ605","PAQ620","PAQ650","PAQ665") %in% names(x))){
    mm <- x[,c("PAQ605","PAQ620","PAQ650","PAQ665")]
    act <- as.integer(rowSums(mm==1,na.rm=TRUE)>0); act[rowSums(!is.na(mm))==0] <- NA
  } else if("PAQ100" %in% names(x)) act <- ifelse(x$PAQ100==1,1L,ifelse(x$PAQ100==2,0L,NA))
  else act <- rep(NA_integer_,nrow(x))
  paq <- rbind(paq, data.frame(SEQN=x$SEQN, cycle=as.integer(y), pa_active=act))
}
d <- merge(d, paq, by=c("SEQN","cycle"), all.x=TRUE)
d$lncot <- log(pmax(d$cot,0.01))
sf <- data.frame()
COVS_nocvd <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn"
COVS_nomed <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + cvd"
runm <- function(dd, covs, ev="event_all", extra=NULL){
  out <- data.frame()
  for(m in MET3){
    f <- paste0("Surv(time,",ev,") ~ z_",m," + ",covs, if(!is.null(extra)) paste0(" + ",extra) else "")
    fit <- safe(svycoxph(as.formula(f), design=mkdes(dd)),paste0("sens ",m)); if(is.null(fit)) next
    r <- getc(fit,paste0("z_",m))
    out <- rbind(out, data.frame(metal=m, HR=r["HR"], lo=r["lo"], hi=r["hi"], p=r["p"]))
  }
  out$n <- nrow(dd); out
}
# main
t0 <- runm(d, COVS_M3); t0$model <- "main(M3)"; sf <- rbind(sf,t0)
# (i) delayed entry 24m vs old-style exclusion
dd <- d[d$time>24,]; dd$start <- 24
for(m in MET3){
  fit <- safe(svycoxph(as.formula(paste0("Surv(start,time,event_all) ~ z_",m," + ",COVS_M3)), design=mkdes(dd)),"delayed")
  if(!is.null(fit)){ r<-getc(fit,paste0("z_",m)); sf<-rbind(sf,data.frame(metal=m,HR=r["HR"],lo=r["lo"],hi=r["hi"],p=r["p"],n=nrow(dd),model="delayed_entry_24m")) }
}
t1 <- runm(d[d$time>24,], COVS_M3); t1$model <- "excl_24m(old)"; sf <- rbind(sf,t1)
# (ii) exclude baseline CVD
t2 <- runm(d[is.na(d$cvd)|d$cvd==0,], COVS_nocvd); t2$model <- "excl_baselineCVD"; sf <- rbind(sf,t2)
# (iii) omit htn & dm
t3 <- runm(d, COVS_nomed); t3$model <- "omit_htn_dm"; sf <- rbind(sf,t3)
# (iv) +cotinine
t4 <- runm(d, COVS_M3, extra="lncot"); t4$model <- "adjust_cotinine"; sf <- rbind(sf,t4)
# (v) +eGFR ; alt exclude eGFR<60
t5 <- runm(d, COVS_M3, extra="egfr"); t5$model <- "adjust_eGFR"; sf <- rbind(sf,t5)
t6 <- runm(d[is.na(d$egfr)|d$egfr>=60,], COVS_M3); t6$model <- "excl_eGFR<60"; sf <- rbind(sf,t6)
# (vi) +alcohol +PA
t7 <- runm(d, COVS_M3, extra="alc + pa_active"); t7$model <- "adjust_alc_PA"; sf <- rbind(sf,t7)
WRo(sf,"sensitivity_all"); print(sf,row.names=FALSE)

# ---------------- (g) subgroup interactions + stratified HRs ----------------
cat("\n== (g) subgroups ==\n")
d$agegrp <- ifelse(d$RIDAGEYR<60,"<60",">=60")
d$smoke2 <- ifelse(d$smoke==2,"current", ifelse(d$smoke==0,"never", ifelse(d$smoke==1,"former",NA)))
sgi <- data.frame(); sgs <- data.frame()
covs_for <- function(sv, base){
  if(sv=="RIAGENDR") return(gsub(" \\+ factor\\(RIAGENDR\\)","",base))
  if(sv=="smoke2")   return(gsub(" \\+ factor\\(smoke\\)","",base))
  base
}
for(m in MET3){
  for(sv in c("RIAGENDR","smoke2","agegrp")){
    cv <- covs_for(sv, COVS_M3)
    f <- as.formula(paste0("Surv(time,event_all) ~ z_",m,"*factor(",sv,") + ",cv))
    fit <- safe(svycoxph(f, design=mkdes(d)),paste0("int ",m,sv)); if(is.null(fit)) next
    co <- summary(fit)$coefficients; introws <- grep(paste0(":factor\\(",sv,"\\)"),rownames(co),value=TRUE)
    if(!length(introws)){ sgi <- rbind(sgi,data.frame(metal=m,strat=sv,p_interaction=NA)); next }
    b<-coef(fit); V<-vcov(fit); idx<-match(introws,names(b))
    stat <- as.numeric(t(b[idx]) %*% solve(V[idx,idx,drop=FALSE]) %*% b[idx]); p<-1-pchisq(stat,length(idx))
    sgi <- rbind(sgi, data.frame(metal=m, strat=sv, p_interaction=p))
    # stratified HRs
    for(lv in sort(unique(na.omit(d[[sv]])))){
      dd <- d[!is.na(d[[sv]]) & d[[sv]]==lv,]
      fitl <- safe(svycoxph(as.formula(paste0("Surv(time,event_all) ~ z_",m," + ",cv)), design=mkdes(dd)),"strat")
      if(!is.null(fitl)){ r<-getc(fitl,paste0("z_",m))
        sgs <- rbind(sgs, data.frame(metal=m, strat=sv, level=as.character(lv), n=nrow(dd), HR=r["HR"],lo=r["lo"],hi=r["hi"],p=r["p"])) }
    }
  }
}
sgi$p_fdr <- p.adjust(sgi$p_interaction, method="BH")
sgi <- sgi[order(sgi$p_interaction),]
WRo(sgi,"subgroup_interaction_fdr"); WRo(sgs,"subgroup_stratified")
print(sgi,row.names=FALSE); print(sgs,row.names=FALSE)

# ---------------- table1 by lead quartile (regenerate) ----------------
des <- mkdes(d); q <- quantile(d$LBXBPB, probs=c(0,.25,.5,.75,1), na.rm=TRUE)
desq <- update(des, Qpb=cut(d$LBXBPB, breaks=q, include.lowest=TRUE, labels=paste0("Q",1:4)))
t1 <- data.frame()
wmean <- function(v) tryCatch(as.numeric(svyby(as.formula(paste0("~",v)), ~Qpb, desq, svymean, na.rm=TRUE)[,2]), error=function(e) rep(NA,4))
for(v in c("RIDAGEYR","BMXBMI","INDFMPIR","sbp","dbp","egfr","tc","hdl","LBXBPB","LBXBCD","LBXTHG")) if(v %in% names(d)){
  r<-wmean(v); t1<-rbind(t1,data.frame(var=v,Q1=r[1],Q2=r[2],Q3=r[3],Q4=r[4])) }
for(vl in list(c("RIAGENDR",2),c("smoke",2),c("dm",1),c("htn",1),c("cvd",1))){
  v<-vl[1]; lv<-vl[2]
  r<-tryCatch(as.numeric(svyby(as.formula(paste0("~I(as.numeric(",v,"==",lv,"))")), ~Qpb, desq, svymean, na.rm=TRUE)[,2]), error=function(e) rep(NA,4))
  t1<-rbind(t1,data.frame(var=paste0(v,"=",lv,"(%)"),Q1=r[1],Q2=r[2],Q3=r[3],Q4=r[4])) }
for(ev in c("event_all","event_cvd")){ r<-as.numeric(svyby(as.formula(paste0("~",ev)), ~Qpb, desq, svymean, na.rm=TRUE)[,2])
  t1<-rbind(t1,data.frame(var=paste0(ev,"(%)"),Q1=100*r[1],Q2=100*r[2],Q3=100*r[3],Q4=100*r[4])) }
WR(t1,"table1_by_leadQ"); print(t1,row.names=FALSE)

saveRDS(d, file.path(D,"analysis_df_plus.RDS"))  # d with PAQ/alc/quartiles for part B
cat("\n=== Part A DONE ===\n")
