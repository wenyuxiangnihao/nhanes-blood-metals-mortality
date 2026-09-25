#!/usr/bin/env Rscript
# ============================================================
# 血重金属 × 死亡  核心分析
#  Table1 + 单金属加权Cox + 混合暴露 quantile g-computation
# ============================================================
.libPaths("/sandbox/workspace/Rlibs")
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu = "adjust")
D  <- "/sandbox/workspace/heavymetal/data"
R  <- "/sandbox/workspace/heavymetal/results"
dir.create(R, showWarnings=FALSE)
d <- readRDS(file.path(D,"analysis_df.rds"))

MET <- c("LBXBPB","LBXBCD","LBXTHG","LBXBSE","LBXBMN")
LBL <- c(LBXBPB="Lead", LBXBCD="Cadmium", LBXTHG="Mercury", LBXBSE="Selenium", LBXBMN="Manganese")
for(m in MET) d[[paste0("ln_",m)]] <- log(d[[m]])

# z-scores
for(m in MET){ v<-d[[paste0("ln_",m)]]; d[[paste0("z_",m)]] <- (v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE) }

COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=d)

# ---------- Table 1: by lead quartiles (weighted) ----------
q <- quantile(d$LBXBPB, probs=c(0,.25,.5,.75,1), na.rm=TRUE)
d$Qpb <- cut(d$LBXBPB, breaks=q, include.lowest=TRUE, labels=paste0("Q",1:4))
des <- update(des, Qpb=cut(d$LBXBPB, breaks=q, include.lowest=TRUE, labels=paste0("Q",1:4)))

t1 <- data.frame()
wmean <- function(v) as.numeric(svyby(as.formula(paste0("~",v)), ~Qpb, des, svymean, na.rm=TRUE)[,2])
for(v in c("RIDAGEYR","BMXBMI","INDFMPIR","SBP2","sbp","dbp","egfr","tc","hdl","LBXBPB","LBXBCD","LBXTHG")){
  if(v %in% names(d)){ r<-tryCatch(wmean(v), error=function(e) rep(NA,4))
    t1 <- rbind(t1, data.frame(var=v, Q1=r[1],Q2=r[2],Q3=r[3],Q4=r[4])) }
}
wpct <- function(v, lvl) as.numeric(svyby(as.formula(paste0("~I(as.numeric(",v,"==",lvl,"))")), ~Qpb, des, svymean, na.rm=TRUE)[,2])
for(vl in list(c("RIAGENDR",2),c("smoke",2),c("dm",1),c("htn",1),c("cvd",1))){
  r <- tryCatch(wpct(vl[1],vl[2]), error=function(e) rep(NA,4))
  t1 <- rbind(t1, data.frame(var=paste0(vl[1],"=",vl[2],"(%)"), Q1=r[1],Q2=r[2],Q3=r[3],Q4=r[4]))
}
# mortality by quartile
for(ev in c("event_all","event_cvd")){
  r <- as.numeric(svyby(as.formula(paste0("~",ev)), ~Qpb, des, svymean, na.rm=TRUE)[,2])
  t1 <- rbind(t1, data.frame(var=paste0(ev,"(%)"), Q1=100*r[1],Q2=100*r[2],Q3=100*r[3],Q4=100*r[4]))
}
write.csv(t1, file.path(R,"table1_by_leadQ.csv"), row.names=FALSE)
cat("== Table 1 (by lead quartile) ==\n"); print(t1, row.names=FALSE)

# ---------- single-metal weighted Cox ----------
res <- data.frame()
for(m in c("LBXBPB","LBXBCD","LBXTHG")){
  for(ev in c("event_all","event_cvd")){
    f <- as.formula(paste0("Surv(time,",ev,") ~ z_",m," + ",COVS))
    fit <- svycoxph(f, design=des)
    s <- summary(fit)$conf.int; p <- summary(fit)$coefficients
    vn <- paste0("z_",m)
    res <- rbind(res, data.frame(metal=LBL[[m]], outcome=ev, term="per-SD(ln)",
      HR=s[vn,"exp(coef)"], lo=s[vn,"lower .95"], hi=s[vn,"upper .95"], p=p[vn,"Pr(>|z|)"]))
  }
}
write.csv(res, file.path(R,"single_metal_cox.csv"), row.names=FALSE)
cat("\n== single-metal Cox ==\n"); print(res, row.names=FALSE)
