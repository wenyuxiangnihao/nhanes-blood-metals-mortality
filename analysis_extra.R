#!/usr/bin/env Rscript
# ============================================================
# RCS 非线性 + 亚组交互 + 敏感性分析
# ============================================================
.libPaths("/sandbox/workspace/Rlibs")
suppressMessages({library(survey); library(survival); library(splines)})
options(survey.lonely.psu="adjust")
D <- "/sandbox/workspace/heavymetal/data"; R <- "/sandbox/workspace/heavymetal/results"
d <- readRDS(file.path(D,"analysis_df.rds"))
MET <- c("LBXBPB","LBXBCD","LBXTHG")
for(m in MET){ v<-log(d[[m]]); d[[paste0("z_",m)]] <- (v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE) }
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
mkdes <- function(dd) svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=dd)

# ---------- P-nonlinear (spline df=3 vs linear) ----------
cat("== P for non-linearity (spline df=3 vs linear, Wald) ==\n")
nl <- data.frame()
for(m in MET){
  for(ev in c("event_all","event_cvd")){
    dd <- d; dd$z2 <- dd[[paste0("z_",m)]]^2
    f2 <- as.formula(paste0("Surv(time,",ev,") ~ z_",m," + z2 + ",COVS))
    fit2 <- svycoxph(f2, design=mkdes(dd))
    b <- coef(fit2); V <- vcov(fit2)
    stat <- unname(b["z2"]^2 / V["z2","z2"]); p <- 1-pchisq(stat, 1)
    nl <- rbind(nl, data.frame(metal=m, outcome=ev, p_nonlinear=p))
  }
}
print(nl, row.names=FALSE); write.csv(nl, file.path(R,"p_nonlinear.csv"), row.names=FALSE)

# ---------- RCS curve (predict vs z, ref at median) ----------
rcs <- data.frame()
for(m in MET){
  ev <- "event_all"
  f2 <- as.formula(paste0("Surv(time,",ev,") ~ ns(z_",m,", df=3) + ",COVS))
  fit2 <- svycoxph(f2, design=mkdes(d))
  zz <- d[[paste0("z_",m)]]; grid <- seq(quantile(zz,.01,na.rm=TRUE), quantile(zz,.99,na.rm=TRUE), length=40)
  nd <- d[rep(1,length(grid)), ]
  for(v in c("time","event_all","event_cvd")) nd[[v]] <- d[[v]][1]
  nd[[paste0("z_",m)]] <- grid
  pr <- tryCatch(predict(fit2, newdata=nd, type="lp", se.fit=TRUE), error=function(e) NULL)
  if(is.null(pr)){ cat("predict failed for",m,"\n"); next }
  ref <- predict(fit2, newdata=transform(nd, zz=median(zz)), type="lp")  # not used
  # reference = median of grid index
  mid <- which.min(abs(grid - median(zz)))
  lp <- pr$fit; se <- pr$se.fit
  hr <- exp(lp - lp[mid]); lo <- exp(lp - lp[mid] - 1.96*sqrt(se^2+se[mid]^2)); hi <- exp(lp - lp[mid] + 1.96*sqrt(se^2+se[mid]^2))
  rcs <- rbind(rcs, data.frame(metal=m, z=grid, HR=hr, lo=lo, hi=hi))
}
write.csv(rcs, file.path(R,"rcs_curves.csv"), row.names=FALSE)

# ---------- subgroup (interaction) ----------
d$agegrp <- ifelse(d$RIDAGEYR < 60, "<60", ">=60")
sub <- data.frame()
for(m in MET){
  for(sv in c("RIAGENDR","smoke","agegrp","dm","htn")){
    f <- as.formula(paste0("Surv(time,event_all) ~ z_",m,"*factor(",sv,") + ",COVS))
    fit <- tryCatch(svycoxph(f, design=mkdes(d)), error=function(e) NULL)
    if(is.null(fit)) next
    co <- summary(fit)$coefficients
    introws <- grep(paste0(":factor\\(",sv,"\\)"), rownames(co), value=TRUE)
    if(!length(introws)) next
    # Wald p for interaction terms
    V <- vcov(fit); b <- coef(fit)
    idx <- match(introws, names(b))
    stat <- as.numeric(t(b[idx]) %*% solve(V[idx,idx,drop=FALSE]) %*% b[idx])
    p <- 1-pchisq(stat, length(idx))
    sub <- rbind(sub, data.frame(metal=m, strat=sv, p_interaction=p))
  }
}
print(sub, row.names=FALSE); write.csv(sub, file.path(R,"subgroup_interaction.csv"), row.names=FALSE)

# ---------- sensitivity: exclude deaths within 24 months ----------
sens <- data.frame()
COVS_nocvd <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn"
for(ev in c("event_all","event_cvd")){
  for(lab in c("main","excl_24m","excl_baselineCVD")){
    dd <- d; cv <- COVS
    if(lab=="excl_24m") dd <- d[d$time > 24, ]
    if(lab=="excl_baselineCVD"){ dd <- d[is.na(d$cvd) | d$cvd==0, ]; cv <- COVS_nocvd }
    for(m in MET){
      f <- as.formula(paste0("Surv(time,",ev,") ~ z_",m," + ",cv))
      fit <- svycoxph(f, design=mkdes(dd)); s <- summary(fit)$conf.int
      vn <- paste0("z_",m)
      sens <- rbind(sens, data.frame(outcome=ev, model=lab, metal=m, n=nrow(dd),
        HR=s[vn,1], lo=s[vn,3], hi=s[vn,4]))
    }
  }
}
print(sens, row.names=FALSE); write.csv(sens, file.path(R,"sensitivity.csv"), row.names=FALSE)
cat("\nDONE\n")
