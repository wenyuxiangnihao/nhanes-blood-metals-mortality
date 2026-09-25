#!/usr/bin/env Rscript
# ============================================================
# 汞的"鱼类/海洋n-3"混杂敏感性分析
#  2011-2014 子集：FAS(血清 EPA+DHA) 作为海产摄入代理
# ============================================================
.libPaths("/sandbox/workspace/Rlibs")
suppressMessages({library(survey); library(survival); library(foreign)})
options(survey.lonely.psu="adjust")
D <- "/sandbox/workspace/heavymetal/data"; RAW <- "/sandbox/workspace/heavymetal/data_raw"
R <- "/sandbox/workspace/heavymetal/results"
d <- readRDS(file.path(D,"analysis_df.rds"))
for(m in c("LBXBPB","LBXBCD","LBXTHG")){ v<-log(d[[m]]); d[[paste0("z_",m)]]<-(v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE) }

# merge FAS (2011-2012 = G, 2013-2014 = H)
fas <- rbind(
  transform(read.xport(file.path(RAW,"FAS_G.XPT"))[, c("SEQN","LBXEPA","LBXDHA","LBXDP6")], cycle=2011),
  transform(read.xport(file.path(RAW,"FAS_H.XPT"))[, c("SEQN","LBXEPA","LBXDHA","LBXDP6")], cycle=2013))
dd <- merge(d[d$cycle %in% c(2011,2013), ], fas, by=c("SEQN","cycle"), all.x=TRUE)
dd$n3 <- dd$LBXEPA + dd$LBXDHA
dd$log_n3 <- log(dd$n3)
cat("2011-2014 subset n =", nrow(dd), " with EPA/DHA =", sum(!is.na(dd$n3)), "\n")

COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
mk <- function(x) svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=x)

res <- data.frame()
sub <- dd[!is.na(dd$n3), ]
for(ev in c("event_all","event_cvd")){
  for(m in c("LBXBPB","LBXBCD","LBXTHG")){
    vn <- paste0("z_",m)
    f0 <- as.formula(paste0("Surv(time,",ev,") ~ ",vn," + ",COVS))
    f1 <- as.formula(paste0("Surv(time,",ev,") ~ ",vn," + log_n3 + ",COVS))
    s0 <- summary(svycoxph(f0, design=mk(sub)))$conf.int
    s1 <- summary(svycoxph(f1, design=mk(sub)))$conf.int
    n3c <- summary(svycoxph(f1, design=mk(sub)))$conf.int["log_n3",]
    res <- rbind(res,
      data.frame(outcome=ev, metal=m, model="not adjusted for n-3", HR=s0[vn,1], lo=s0[vn,3], hi=s0[vn,4]),
      data.frame(outcome=ev, metal=m, model="adjusted for n-3",    HR=s1[vn,1], lo=s1[vn,3], hi=s1[vn,4]),
      data.frame(outcome=ev, metal=m, model="n-3 (per log unit)",  HR=n3c[1], lo=n3c[3], hi=n3c[4]))
  }
}
print(res, row.names=FALSE); write.csv(res, file.path(R,"mercury_n3_sensitivity.csv"), row.names=FALSE)

# also: adjust mercury for blood selenium (2011-2018 subset)
sub2 <- d[!is.na(d$LBXBSE), ]; sub2$ln_se <- log(sub2$LBXBSE)
res2 <- data.frame()
for(ev in c("event_all","event_cvd")){
  f0 <- as.formula(paste0("Surv(time,",ev,") ~ z_LBXTHG + ",COVS))
  f1 <- as.formula(paste0("Surv(time,",ev,") ~ z_LBXTHG + ln_se + ",COVS))
  s0 <- summary(svycoxph(f0, design=mk(sub2)))$conf.int
  s1 <- summary(svycoxph(f1, design=mk(sub2)))$conf.int
  res2 <- rbind(res2,
    data.frame(outcome=ev, model="Hg not adj Se", HR=s0["z_LBXTHG",1], lo=s0["z_LBXTHG",3], hi=s0["z_LBXTHG",4]),
    data.frame(outcome=ev, model="Hg adj Se",     HR=s1["z_LBXTHG",1], lo=s1["z_LBXTHG",3], hi=s1["z_LBXTHG",4]))
}
print(res2, row.names=FALSE); write.csv(res2, file.path(R,"mercury_se_sensitivity.csv"), row.names=FALSE)
cat("DONE\n")
