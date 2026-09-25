#!/usr/bin/env Rscript
# ============================================================
# Round-2 review, Major #4: never-smoker sensitivity analyses
#   (a) main per-SD(ln) Cox models restricted to never smokers
#   (b) cause-specific models (esp. cadmium-CLRD, mercury-cancer)
#   (c) 3-metal qgcomp + directional decomposition among never smokers
#   (d) further restriction: never smokers with serum cotinine < 10 ng/mL
#   (e) cause-specific cadmium-CLRD stratified by smoking status
# ============================================================
.libPaths("/sandbox/workspace/Rlibs")
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu="adjust")
D <- "/sandbox/workspace/heavymetal/data"; R <- "/sandbox/workspace/heavymetal/results"
d <- readRDS(file.path(D,"analysis_df.rds"))
MET3 <- c("LBXBPB","LBXBCD","LBXTHG")
for(m in MET3){ d[[paste0("ln_",m)]] <- log(d[[m]]) }
for(m in MET3){ v<-d[[paste0("ln_",m)]]; d[[paste0("z_",m)]] <- (v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE) }

COVS_NS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + dm + htn + cvd"
mkdes <- function(dd) svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=dd)
getc <- function(fit,term){ s<-summary(fit)$conf.int; co<-summary(fit)$coefficients
  c(HR=s[term,1], lo=s[term,3], hi=s[term,4], p=co[term,ncol(co)]) }
safe <- function(expr,label){ tryCatch(expr, error=function(e){ cat("  !!! ERROR in",label,":",conditionMessage(e),"\n"); NULL }) }

dns <- d[!is.na(d$smoke) & d$smoke==0, ]
cat("never-smoker analytic n =", nrow(dns), " all-deaths =", sum(dns$event_all), " cvd =", sum(dns$event_cvd), "\n")

# ---------------- (a) main models ----------------
cat("\n== (a) never-smoker main models (per SD of ln) ==\n")
res <- data.frame()
for(ev in c("event_all","event_cvd")) for(m in MET3){
  f <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_",m," + ",COVS_NS)), design=mkdes(dns)),paste(ev,m))
  if(is.null(f)) next
  cc <- getc(f,paste0("z_",m))
  res <- rbind(res, data.frame(subset="never-smoker", outcome=ev, metal=m, n=nrow(dns), n_deaths=sum(dns[[ev]]),
    HR=cc["HR"], lo=cc["lo"], hi=cc["hi"], p=cc["p"]))
}
print(res,row.names=FALSE)

# ---------------- (d) never smokers with cotinine < 10 ng/mL ----------------
cat("\n== (d) never smokers, cotinine < 10 ng/mL ==\n")
dns2 <- dns[!is.na(dns$cot) & dns$cot < 10, ]
cat("n =", nrow(dns2), " all-deaths =", sum(dns2$event_all), "\n")
res2 <- data.frame()
for(ev in c("event_all","event_cvd")) for(m in MET3){
  f <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_",m," + ",COVS_NS)), design=mkdes(dns2)),paste(ev,m))
  if(is.null(f)) next
  cc <- getc(f,paste0("z_",m))
  res2 <- rbind(res2, data.frame(subset="never-smoker+cot<10", outcome=ev, metal=m, n=nrow(dns2),
    n_deaths=sum(dns2[[ev]]), HR=cc["HR"], lo=cc["lo"], hi=cc["hi"], p=cc["p"]))
}
print(res2,row.names=FALSE)
write.csv(rbind(res,res2), file.path(R,"neversmoker_main_10cyc.csv"), row.names=FALSE)

# ---------------- (b) cause-specific among never smokers ----------------
cat("\n== (b) never-smoker cause-specific ==\n")
causes <- list("heart_disease"=1, "cancer"=2, "CLRD"=3, "cerebrovascular"=5)
for(cn in names(causes)) dns[[paste0("ev_",cn)]] <- as.integer(!is.na(dns$mortstat)&dns$mortstat==1&!is.na(dns$ucod_leading)&dns$ucod_leading==causes[[cn]])
sj <- data.frame()
for(cn in names(causes)) for(m in MET3){
  if(sum(dns[[paste0("ev_",cn)]]) < 15) { cat("  skip",cn,m,"(",sum(dns[[paste0("ev_",cn)]]),"deaths)\n"); next }
  f <- safe(svycoxph(as.formula(paste0("Surv(time,ev_",cn,") ~ z_",m," + ",COVS_NS)), design=mkdes(dns)),paste(cn,m))
  if(is.null(f)) next
  cc <- getc(f,paste0("z_",m))
  sj <- rbind(sj, data.frame(cause=cn, metal=m, n=nrow(dns), n_deaths=sum(dns[[paste0("ev_",cn)]]),
    HR=cc["HR"], lo=cc["lo"], hi=cc["hi"], p=cc["p"]))
}
print(sj,row.names=FALSE)
write.csv(sj, file.path(R,"neversmoker_cause_specific_10cyc.csv"), row.names=FALSE)

# ---------------- (e) Cd-CLRD by smoking status (all subjects) ----------------
cat("\n== (e) cadmium-CLRD stratified by smoking status ==\n")
for(cn in names(causes)) d[[paste0("ev_",cn)]] <- as.integer(!is.na(d$mortstat)&d$mortstat==1&!is.na(d$ucod_leading)&d$ucod_leading==causes[[cn]])
COVS_SM <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + dm + htn + cvd"
st <- data.frame()
for(sv in c(0,1,2)){
  ds <- d[!is.na(d$smoke)&d$smoke==sv, ]
  f <- safe(svycoxph(as.formula(paste0("Surv(time,ev_CLRD) ~ z_LBXBCD + ",COVS_SM)), design=mkdes(ds)),paste("CLRD smk",sv))
  if(is.null(f)) next
  cc <- getc(f,"z_LBXBCD")
  st <- rbind(st, data.frame(smoking=c("never","former","current")[sv+1], n=nrow(ds), n_CLRD_deaths=sum(ds$ev_CLRD),
    HR=cc["HR"], lo=cc["lo"], hi=cc["hi"], p=cc["p"]))
}
print(st,row.names=FALSE)
write.csv(st, file.path(R,"cadmium_clrd_by_smoking_10cyc.csv"), row.names=FALSE)

# ---------------- (c) qgcomp + directional, never smokers ----------------
cat("\n== (c) never-smoker qgcomp + directional (3-metal) ==\n")
qg <- function(dat, exps, event){
  dq <- dat
  for(m in exps){ br<-unique(quantile(dat[[m]],probs=seq(0,1,0.25),na.rm=TRUE)); br[1]<- -Inf; br[length(br)]<-Inf
    dq[[paste0("qc_",m)]] <- as.integer(cut(dat[[m]],breaks=br,labels=FALSE,include.lowest=TRUE))-1 }
  qc <- paste0("qc_",exps); dq$S <- rowSums(dq[,qc])
  desq <- mkdes(dq)
  fF <- as.formula(paste0("Surv(time,",event,") ~ ",paste(qc,collapse="+")," + ",COVS_NS))
  ff <- svycoxph(fF,design=desq); b <- coef(ff)[qc]; V <- vcov(ff)[qc,qc]
  psi <- sum(b); v <- sum(V)
  bp <- b[b>0]; bn <- b[b<0]
  sump <- if(length(bp)) sum(bp) else 0; sumn <- if(length(bn)) sum(bn) else 0
  Vpp <- if(length(bp)) sum(V[names(bp),names(bp),drop=FALSE]) else 0
  Vnn <- if(length(bn)) sum(V[names(bn),names(bn),drop=FALSE]) else 0
  list(psi=psi, se=sqrt(v), sump=sump, sepp=sqrt(Vpp), sumn=sumn, senn=sqrt(Vnn),
       b=b, n=nrow(dq), n_deaths=sum(dq[[event]]))
}
qgout <- data.frame()
for(ev in c("event_all","event_cvd")){
  r <- safe(qg(dns, MET3, ev), paste("qg ns",ev)); if(is.null(r)) next
  qgout <- rbind(qgout, data.frame(subset="never-smoker", outcome=ev, n=r$n, n_deaths=r$n_deaths,
    HR_overall=exp(r$psi), ov_lo=exp(r$psi-1.96*r$se), ov_hi=exp(r$psi+1.96*r$se),
    p_overall=2*pnorm(-abs(r$psi/r$se)),
    HR_pos=exp(r$sump), pos_lo=exp(r$sump-1.96*r$sepp), pos_hi=exp(r$sump+1.96*r$sepp),
    HR_neg=exp(r$sumn), neg_lo=exp(r$sumn-1.96*r$senn), neg_hi=exp(r$sumn+1.96*r$senn)))
}
print(qgout,row.names=FALSE)
write.csv(qgout, file.path(R,"neversmoker_qgcomp_10cyc.csv"), row.names=FALSE)
cat("\nDONE.\n")
