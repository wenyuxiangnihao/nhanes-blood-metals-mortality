#!/usr/bin/env Rscript
# ============================================================
# RE-ANALYSIS (10 cycles) -- Part B: tasks 2(h)-(l)
# ============================================================
.libPaths("/sandbox/workspace/Rlibs")
suppressMessages({library(survey); library(survival); library(foreign); library(splines)})
options(survey.lonely.psu="adjust")
D <- "/sandbox/workspace/heavymetal/data"; R <- "/sandbox/workspace/heavymetal/results"; RAW <- "/sandbox/workspace/heavymetal/data_raw"
d <- readRDS(file.path(D,"analysis_df_plus.RDS"))
cat("analytic n =", nrow(d), "\n")
MET3 <- c("LBXBPB","LBXBCD","LBXTHG"); MET5 <- c(MET3,"LBXBSE","LBXBMN")
for(m in MET5) d[[paste0("ln_",m)]] <- log(d[[m]])
COVS_M3 <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
mkdes <- function(dd) svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=dd)
WRo <- function(df,name) write.csv(df,file.path(R,paste0(name,"_10cyc.csv")),row.names=FALSE)
WR  <- function(df,name){ write.csv(df,file.path(R,paste0(name,"_10cyc.csv")),row.names=FALSE); write.csv(df,file.path(R,paste0(name,".csv")),row.names=FALSE) }
getc <- function(fit,term){ s<-summary(fit)$conf.int; co<-summary(fit)$coefficients; c(HR=s[term,1],lo=s[term,3],hi=s[term,4],p=co[term,ncol(co)]) }
safe <- function(e,l){ tryCatch(e, error=function(x){cat("  !!! ERROR",l,":",conditionMessage(x),"\n");NULL}) }

# rcs helper (4 knots / df=3)
rcs_fit <- function(dd, m, event, ngrid=40){
  zz <- dd[[paste0("z_",m)]]; ok <- !is.na(zz)
  bs <- ns(zz[ok], df=3); dd$s1<-NA;dd$s2<-NA;dd$s3<-NA; dd$s1[ok]<-bs[,1]; dd$s2[ok]<-bs[,2]; dd$s3[ok]<-bs[,3]
  fit <- safe(svycoxph(as.formula(paste0("Surv(time,",event,") ~ s1+s2+s3 + ",COVS_M3)), design=mkdes(dd)),"rcs")
  if(is.null(fit)) return(NULL)
  grid <- seq(quantile(zz,.01,na.rm=TRUE), quantile(zz,.99,na.rm=TRUE), length=ngrid)
  B <- predict(bs, grid)
  b <- coef(fit)[c("s1","s2","s3")]; V <- vcov(fit)[c("s1","s2","s3"),c("s1","s2","s3")]
  lp <- as.numeric(B %*% b); se <- sqrt(rowSums((B %*% V) * B)); ref <- which.min(abs(grid - median(zz,na.rm=TRUE)))
  data.frame(metal=m, z=grid, HR=exp(lp-lp[ref]), lo=exp(lp-lp[ref]-1.96*sqrt(se^2+se[ref]^2)), hi=exp(lp-lp[ref]+1.96*sqrt(se^2+se[ref]^2)))
}

# ---------------- (h) Se & Mn ----------------
cat("\n== (h) Se & Mn ==\n")
sub <- d[!is.na(d$LBXBSE) & !is.na(d$LBXBMN),]
sub$z_LBXBSE <- (sub$ln_LBXBSE-mean(sub$ln_LBXBSE,na.rm=TRUE))/sd(sub$ln_LBXBSE,na.rm=TRUE)
sub$z_LBXBMN <- (sub$ln_LBXBMN-mean(sub$ln_LBXBMN,na.rm=TRUE))/sd(sub$ln_LBXBMN,na.rm=TRUE)
cat("Se/Mn subset n =", nrow(sub), " all-deaths =", sum(sub$event_all), " cvd =", sum(sub$event_cvd), "\n")
sh <- data.frame(); rcs_all <- data.frame()
for(m in c("LBXBSE","LBXBMN")){
  mnln <- mean(sub[[paste0("ln_",m)]],na.rm=TRUE); sdln <- sd(sub[[paste0("ln_",m)]],na.rm=TRUE)
  for(ev in c("event_all","event_cvd")){
    s <- getc(safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_",m," + ",COVS_M3)), design=mkdes(sub)),"linear"),paste0("z_",m))
    sq <- sub; sq$z1<-sq[[paste0("z_",m)]]; sq$z2<-sq$z1^2
    f2 <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z1+z2 + ",COVS_M3)), design=mkdes(sq)),"quad")
    b<-coef(f2); V<-vcov(f2); p <- 1-pchisq(unname(b["z2"]^2/V["z2","z2"]),1)
    b1<-b["z1"]; b2<-b["z2"]; znad <- if(b2!=0) -b1/(2*b2) else NA
    conc_nad <- if(!is.na(znad)) exp(mnln + znad*sdln) else NA
    sh <- rbind(sh, data.frame(subset="2011-2018", metal=m, outcome=ev, n=nrow(sub),
      HR_perSD=s["HR"], lo=s["lo"], hi=s["hi"], p=s["p"], p_nonlinear=p,
      nadir_z=znad, nadir_conc=conc_nad, mean_ln=mnln, sd_ln=sdln))
  }
  rr <- rcs_fit(sub, m, "event_all"); if(!is.null(rr)) rcs_all <- rbind(rcs_all, rr)
}
WR(sh,"selenium_manganese"); WRo(rcs_all,"rcs_se_mn"); print(sh,row.names=FALSE)

# ---------------- (i) period-specific ----------------
cat("\n== (i) period-specific ==\n")
periods <- list("1999-2010"=c(1999,2001,2003,2005,2007,2009), "2011-2014"=c(2011,2013), "2015-2018"=c(2015,2017))
si <- data.frame()
for(pn in names(periods)){
  dp <- d[d$cycle %in% periods[[pn]], ]
  for(m in MET3){
    r <- getc(safe(svycoxph(as.formula(paste0("Surv(time,event_all) ~ z_",m," + ",COVS_M3)), design=mkdes(dp)),"period"),paste0("z_",m))
    si <- rbind(si, data.frame(period=pn, metal=m, n=nrow(dp), HR=r["HR"], lo=r["lo"], hi=r["hi"], p=r["p"]))
  }
}
WR(si,"period_stability"); print(si,row.names=FALSE)

# ---------------- (j) cause-specific ----------------
cat("\n== (j) cause-specific ==\n")
causes <- list("heart_disease"=1, "cancer"=2, "cerebrovascular"=5, "CLRD"=3)
for(cn in names(causes)){ d[[paste0("ev_",cn)]] <- as.integer(!is.na(d$mortstat)&d$mortstat==1&!is.na(d$ucod_leading)&d$ucod_leading==causes[[cn]]) }
sj <- data.frame()
for(cn in names(causes)){
  ev <- paste0("ev_",cn)
  cat(sprintf("  %-14s deaths = %d\n", cn, sum(d[[ev]])))
  for(m in MET3){
    r <- getc(safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_",m," + ",COVS_M3)), design=mkdes(d)),"cause"),paste0("z_",m))
    sj <- rbind(sj, data.frame(cause=cn, ucod_code=causes[[cn]], metal=m, n_deaths=sum(d[[ev]]),
      HR=r["HR"], lo=r["lo"], hi=r["hi"], p=r["p"]))
  }
}
WRo(sj,"cause_specific"); print(sj,row.names=FALSE)

# ---------------- (k) Se:Hg molar ratio ----------------
cat("\n== (k) Se:Hg molar ratio ==\n")
sk <- data.frame()
molar <- function(x) x
d$ratio_sehg <- (d$LBXBSE/78.96)/(d$LBXTHG/200.59)
d$ln_ratio <- log(d$ratio_sehg)
k <- d[!is.na(d$ln_ratio),]
k$z_ratio <- (k$ln_ratio-mean(k$ln_ratio,na.rm=TRUE))/sd(k$ln_ratio,na.rm=TRUE)
k$z_LBXTHG <- (k$ln_LBXTHG-mean(k$ln_LBXTHG,na.rm=TRUE))/sd(k$ln_LBXTHG,na.rm=TRUE)
for(ev in c("event_all","event_cvd")){
  r1 <- getc(safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_ratio + ",COVS_M3)), design=mkdes(k)),"ratio"),"z_ratio")
  sk <- rbind(sk, data.frame(subset="2011-2018", outcome=ev, model="Se:Hg_molar_ratio_perSD", n=nrow(k), HR=r1["HR"],lo=r1["lo"],hi=r1["hi"],p=r1["p"]))
  r2 <- getc(safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_LBXTHG + ",COVS_M3)), design=mkdes(k)),"Hg"),"z_LBXTHG")
  sk <- rbind(sk, data.frame(subset="2011-2018", outcome=ev, model="Hg_perSD (subset)", n=nrow(k), HR=r2["HR"],lo=r2["lo"],hi=r2["hi"],p=r2["p"]))
  r3 <- getc(safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ z_LBXTHG + z_ratio + ",COVS_M3)), design=mkdes(k)),"Hgratio"),"z_LBXTHG")
  sk <- rbind(sk, data.frame(subset="2011-2018", outcome=ev, model="Hg_perSD + adj_ratio", n=nrow(k), HR=r3["HR"],lo=r3["lo"],hi=r3["hi"],p=r3["p"]))
}
WRo(sk,"se_hg_ratio"); print(sk,row.names=FALSE)

# ---------------- (l) WTSH2YR ----------------
cat("\n== (l) WTSH2YR weight sensitivity ==\n")
sl <- data.frame()
found <- c()
for(f in c("DEMO_H","DEMO_I","DEMO_G","DEMO_J")){
  p <- file.path(RAW,paste0(f,".XPT")); if(!file.exists(p)) next
  x <- tryCatch(read.xport(p),error=function(e)NULL); if(is.null(x)) next
  if("WTSH2YR" %in% names(x)) found <- c(found,f)
}
if(!length(found)){
  sl <- data.frame(note="WTSH2YR not present in DEMO_H/DEMO_I (2013-2016): weight sensitivity not derivable")
  cat(sl$note,"\n")
} else {
  sl <- data.frame(note=paste("WTSH2YR found in:",paste(found,collapse=",")))
}
WRo(sl,"weight_sensitivity_WTSH2YR")

cat("\n=== Part B DONE ===\n")
