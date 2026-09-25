#!/usr/bin/env Rscript
.libPaths("/sandbox/workspace/Rlibs")
suppressMessages({library(survey); library(survival); library(splines)})
options(survey.lonely.psu="adjust")
D <- "/sandbox/workspace/heavymetal/data"; R <- "/sandbox/workspace/heavymetal/results"
d <- readRDS(file.path(D,"analysis_df_plus.RDS"))
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
mkdes <- function(dd) svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=dd)
rcs_fit <- function(dd,m,event,ngrid=40){
  zz<-dd[[paste0("z_",m)]]; ok<-!is.na(zz); bs<-ns(zz[ok],df=3)
  dd$s1<-NA;dd$s2<-NA;dd$s3<-NA; dd$s1[ok]<-bs[,1];dd$s2[ok]<-bs[,2];dd$s3[ok]<-bs[,3]
  fit<-svycoxph(as.formula(paste0("Surv(time,",event,") ~ s1+s2+s3 + ",COVS)),design=mkdes(dd))
  grid<-seq(quantile(zz,.01,na.rm=TRUE),quantile(zz,.99,na.rm=TRUE),length=ngrid); B<-predict(bs,grid)
  b<-coef(fit)[c("s1","s2","s3")]; V<-vcov(fit)[c("s1","s2","s3"),c("s1","s2","s3")]
  lp<-as.numeric(B%*%b); se<-sqrt(rowSums((B%*%V)*B)); ref<-which.min(abs(grid-median(zz,na.rm=TRUE)))
  data.frame(metal=m, z=grid, HR=exp(lp-lp[ref]), lo=exp(lp-lp[ref]-1.96*sqrt(se^2+se[ref]^2)), hi=exp(lp-lp[ref]+1.96*sqrt(se^2+se[ref]^2)))
}
out<-data.frame()
for(m in c("LBXBPB","LBXBCD","LBXTHG")) out<-rbind(out, rcs_fit(d,m,"event_all"))
write.csv(out, file.path(R,"rcs_curves_10cyc.csv"), row.names=FALSE)
write.csv(out, file.path(R,"rcs_curves.csv"), row.names=FALSE)
cat("rcs3 done, rows=",nrow(out),"\n")
