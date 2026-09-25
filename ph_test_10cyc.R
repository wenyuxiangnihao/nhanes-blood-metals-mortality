.libPaths("/sandbox/workspace/Rlibs")
suppressMessages({library(survey);library(survival)})
options(survey.lonely.psu="adjust")
d <- readRDS("/sandbox/workspace/heavymetal/data/analysis_df.rds")
for(m in c("LBXBPB","LBXBCD","LBXTHG")){ v<-log(d[[m]]); d[[paste0("z_",m)]]<-(v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE) }
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
mkd <- function(dd) svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~wt,nest=TRUE,data=dd)
out<-data.frame()
for(m in c("LBXBPB","LBXBCD","LBXTHG")){
  d$lt <- log(d$time+0.5); d$xt <- d[[paste0("z_",m)]]*d$lt
  f <- svycoxph(as.formula(paste0("Surv(time,event_all) ~ z_",m," + xt + ",COVS)),design=mkd(d))
  co <- summary(f)$coefficients; p_int <- co["xt",ncol(co)]
  f0 <- coxph(as.formula(paste0("Surv(time,event_all) ~ z_",m," + ",COVS)),data=d)
  z <- cox.zph(f0); pgl <- z$table["GLOBAL","p"]
  out <- rbind(out,data.frame(metal=m,p_PH_exposure=p_int,p_PH_global=pgl))
}
print(out,row.names=FALSE); write.csv(out,"/sandbox/workspace/heavymetal/results/ph_test_10cyc.csv",row.names=FALSE)
