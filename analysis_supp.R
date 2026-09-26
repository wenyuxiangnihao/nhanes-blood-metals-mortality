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
# 补充分析：Cox 比例风险假定(暴露xlog时间交互) + 交互检验 FDR + 子集样本量
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu="adjust")
D <- file.path(ROOT, "data"); R <- file.path(ROOT, "results")
d <- readRDS(file.path(D,"analysis_df.rds"))
MET <- c("LBXBPB","LBXBCD","LBXTHG")
for(m in MET){ v <- log(d[[m]]); d[[paste0("z_",m)]] <- (v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE) }
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
mk <- function(x) svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=x)
d$lnt <- log(d$time)

cat("=== PH assumption (cox.zph on unweighted analogue) ===\n")
ph <- data.frame()
for(m in MET){
  fit2 <- coxph(as.formula(paste0("Surv(time,event_all) ~ z_",m," + ",COVS)), data=d)
  z <- cox.zph(fit2)
  pv <- z$table[paste0("z_",m), "p"]
  cat(sprintf("  %-8s exposure PH p = %.4f | global p = %.4f\n", m, pv, z$table["GLOBAL","p"]))
  ph <- rbind(ph, data.frame(metal=m, p_PH_exposure=pv, p_PH_global=z$table["GLOBAL","p"]))
}
write.csv(ph, file.path(R,"ph_test.csv"), row.names=FALSE)

cat("\n=== interaction tests with BH-FDR ===\n")
inr <- read.csv(file.path(R,"subgroup_interaction.csv"))
inr$p_fdr <- p.adjust(inr$p_interaction, method="BH")
inr <- inr[order(inr$p_interaction), ]
print(inr, row.names=FALSE); write.csv(inr, file.path(R,"subgroup_interaction_fdr.csv"), row.names=FALSE)

cat("\n=== subset sizes ===\n")
cat("main  n=", nrow(d), " deaths all=", sum(d$event_all), " cvd=", sum(d$event_cvd), "\n")
s5 <- d[!is.na(d$LBXBSE) & !is.na(d$LBXBMN), ]
cat("Se/Mn subset n=", nrow(s5), " deaths all=", sum(s5$event_all), " cvd=", sum(s5$event_cvd), "\n")

cat("DONE\n")
