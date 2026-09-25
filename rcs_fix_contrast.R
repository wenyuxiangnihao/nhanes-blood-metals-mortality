#!/usr/bin/env Rscript
# ============================================================
# FIX (round 4, C1): restricted cubic splines with the CORRECT
# contrast variance.  Previously the CI used  sqrt(se^2 + se[ref]^2),
# which does not vanish at the reference point; the contrast
#   d(x) = B(x) - B(ref)
# must be used so that Var = d V d', which is exactly 0 at x = ref.
# ============================================================
.libPaths("/sandbox/workspace/Rlibs")
suppressMessages({library(survey); library(survival); library(splines)})
options(survey.lonely.psu="adjust")
D <- "/sandbox/workspace/heavymetal/data"; R <- "/sandbox/workspace/heavymetal/results"
d <- readRDS(file.path(D,"analysis_df.rds"))
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
mkdes <- function(dd) svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=dd)
for(m in c("LBXBPB","LBXBCD","LBXTHG","LBXBSE","LBXBMN")){ v<-log(d[[m]]); d[[paste0("z_",m)]]<-(v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE) }

rcs_fit <- function(dd, m, event, ngrid=101){
  zz <- dd[[paste0("z_",m)]]; ok <- !is.na(zz)
  bs <- ns(zz[ok], df=3)
  dd$s1<-NA; dd$s2<-NA; dd$s3<-NA
  dd$s1[ok]<-bs[,1]; dd$s2[ok]<-bs[,2]; dd$s3[ok]<-bs[,3]
  fit <- svycoxph(as.formula(paste0("Surv(time,",event,") ~ s1+s2+s3 + ",COVS)), design=mkdes(dd))
  # grid: 1st to 99th percentile of the exposure (no extrapolation)
  grid <- seq(quantile(zz,.01,na.rm=TRUE), quantile(zz,.99,na.rm=TRUE), length=ngrid)
  B <- predict(bs, grid)
  b <- coef(fit)[c("s1","s2","s3")]; V <- vcov(fit)[c("s1","s2","s3"),c("s1","s2","s3")]
  ref <- which.min(abs(grid - median(zz, na.rm=TRUE)))
  Dm  <- sweep(B, 2, B[ref,], "-")          # B(x) - B(median): contrast basis
  lp  <- as.numeric(Dm %*% b)
  se  <- sqrt(pmax(rowSums((Dm %*% V) * Dm), 0))
  data.frame(metal=m, outcome=event, z=grid, HR=exp(lp),
             lo=exp(lp-1.96*se), hi=exp(lp+1.96*se), se_logHR=se,
             ref_z=grid[ref])
}
all3 <- do.call(rbind, lapply(c("LBXBPB","LBXBCD","LBXTHG"), function(m) rcs_fit(d, m, "event_all")))
sub  <- d[!is.na(d$LBXBSE) & !is.na(d$LBXBMN), ]
se2  <- do.call(rbind, lapply(c("LBXBSE","LBXBMN"), function(m) rcs_fit(sub, m, "event_all")))
write.csv(all3, file.path(R,"rcs_curves_10cyc.csv"), row.names=FALSE)
write.csv(all3, file.path(R,"rcs_curves.csv"), row.names=FALSE)
write.csv(se2,  file.path(R,"rcs_se_mn_10cyc.csv"), row.names=FALSE)
write.csv(se2,  file.path(R,"rcs_se_mn.csv"), row.names=FALSE)

# verification: CI must be exactly 1 at the reference point
cat("\n== sanity check: HR and CI at the reference (median) ==\n")
for(df in list(all3, se2)){
  for(m in unique(df$metal)){
    r <- df[df$metal==m,]; i <- which.min(abs(r$z - r$ref_z[1]))
    cat(sprintf("  %-7s HR=%.4f  95%%CI %.4f-%.4f  SE=%.3g\n", m, r$HR[i], r$lo[i], r$hi[i], r$se_logHR[i]))
  }
}
cat("\n== widest CI (should now be modest) ==\n")
for(df in list(all3, se2)) for(m in unique(df$metal)){
  r <- df[df$metal==m,]; cat(sprintf("  %-7s max upper limit = %.2f\n", m, max(r$hi)))
}
