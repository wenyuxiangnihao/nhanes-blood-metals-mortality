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
# 硒/锰 及 5 金属补充分析（2011-2018 子集）
# ============================================================
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu="adjust")
D <- file.path(ROOT, "data"); R <- file.path(ROOT, "results")
d <- readRDS(file.path(D,"analysis_df.rds"))
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
mk <- function(x) svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=x)
for(m in c("LBXBPB","LBXBCD","LBXTHG","LBXBSE","LBXBMN")) d[[paste0("z_",m)]] <- as.numeric(scale(log(d[[m]])))

sub <- d[!is.na(d$LBXBSE) & !is.na(d$LBXBMN), ]
cat("2011-2018 Se/Mn subset n =", nrow(sub), " deaths(all) =", sum(sub$event_all), " cvd =", sum(sub$event_cvd), "\n")

res <- data.frame()
for(m in c("LBXBSE","LBXBMN")){
  for(ev in c("event_all","event_cvd")){
    f <- as.formula(paste0("Surv(time,",ev,") ~ z_",m," + ",COVS))
    s <- summary(svycoxph(f, design=mk(sub)))$conf.int; vn <- paste0("z_",m)
    res <- rbind(res, data.frame(subset="2011-2018", metal=m, outcome=ev, term="per-SD(ln)",
      HR=s[vn,1], lo=s[vn,3], hi=s[vn,4]))
    # nonlinearity (quadratic)
    sq <- sub; sq$z2 <- sq[[paste0("z_",m)]]^2
    f2 <- as.formula(paste0("Surv(time,",ev,") ~ z_",m," + z2 + ",COVS))
    fit2 <- svycoxph(f2, design=mk(sq)); b<-coef(fit2); V<-vcov(fit2)
    p <- 1-pchisq(unname(b["z2"]^2/V["z2","z2"]), 1)
    res <- rbind(res, data.frame(subset="2011-2018", metal=m, outcome=ev, term="p_nonlinear",
      HR=NA, lo=NA, hi=p))
  }
}
print(res, row.names=FALSE); write.csv(res, file.path(R,"selenium_manganese.csv"), row.names=FALSE)

# per-cycle stability of lead/cadmium/mercury (2011-2014 vs 2015-2018)
stab <- data.frame()
periods <- list("2011-2014"=c(2011,2013), "2015-2018"=c(2015,2017), "1999-2010"=c(1999,2001,2005,2007,2009))
for(pn in names(periods)){
  dp <- d[d$cycle %in% periods[[pn]], ]
  for(m in c("LBXBPB","LBXBCD","LBXTHG")){
    fe <- as.formula(paste0("Surv(time,event_all) ~ z_",m," + ",COVS))
    s <- summary(svycoxph(fe, design=mk(dp)))$conf.int; vn <- paste0("z_",m)
    stab <- rbind(stab, data.frame(period=pn, metal=m, n=nrow(dp), HR=s[vn,1], lo=s[vn,3], hi=s[vn,4]))
  }
}
print(stab, row.names=FALSE); write.csv(stab, file.path(R,"period_stability.csv"), row.names=FALSE)
cat("DONE\n")
