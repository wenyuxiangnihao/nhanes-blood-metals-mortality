#!/usr/bin/env Rscript
# ============================================================
# Task 4: BKMR (binomial/probit) on a random subsample N=800, iter=2000
# ============================================================
.libPaths("/sandbox/workspace/Rlibs")
suppressMessages(library(bkmr))
D <- "/sandbox/workspace/heavymetal/data"; R <- "/sandbox/workspace/heavymetal/results"; F <- "/sandbox/workspace/heavymetal/figures"
set.seed(20260917)
ITER <- 2000; NMAX <- 800
d <- readRDS(file.path(D,"analysis_df.rds"))
cat("Full N =", nrow(d), " deaths =", sum(d$event_all), "\n")
if(nrow(d) > NMAX){ d <- d[sample(nrow(d), NMAX), ]; cat("Subsampled N =", nrow(d), " deaths =", sum(d$event_all), "\n") }
metals <- c("LBXBPB","LBXBCD","LBXTHG")
n <- nrow(d)
Z <- sapply(metals, function(m) qnorm((rank(d[[m]], ties.method="first")-0.5)/n))
colnames(Z) <- c("Lead","Cadmium","Mercury")
FORM <- ~ RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd + log(pmax(time,1))
X <- model.matrix(FORM, data=d)[, -1, drop=FALSE]
X <- X[, apply(X,2,function(v) sd(v,na.rm=TRUE))>1e-8, drop=FALSE]; X <- scale(X)
y <- d$event_all
cat("Estimated runtime ~", round(2.2e-9*n^3*ITER/60,1), "min\n"); t0 <- Sys.time()
fit <- kmbayes(y=y, Z=Z, X=X, iter=ITER, family="binomial", verbose=TRUE)
cat("fit done in", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "min\n")
saveRDS(fit, file.path(R,"bkmr_fit.rds"))
qs <- seq(0.1,0.9,by=0.1)
ov <- tryCatch(OverallRiskSummaries(fit, qs=qs, q.fixed=0.5, method="approx"), error=function(e) OverallRiskSummaries(fit, qs=qs, q.fixed=0.5))
write.csv(ov, file.path(R,"bkmr_overall_summaries.csv"), row.names=FALSE)
pr <- tryCatch(PredictorResponseUnivar(fit, qs=qs, method="approx"), error=function(e) PredictorResponseUnivar(fit, qs=qs))
write.csv(pr, file.path(R,"bkmr_single_summaries.csv"), row.names=FALSE)
png(file.path(F,"bkmr_overall.png"), width=1400, height=1100, res=180)
est <- ov[[grep("est",names(ov),ignore.case=TRUE)[1]]]; qy <- ov[[grep("^q$|quantile",names(ov),ignore.case=TRUE)[1]]]
plot(est, qy, type="b", pch=19, xlab="effect (probit scale)", ylab="Quantile of all metals (vs median)", main="BKMR overall mixture effect")
abline(v=0, lty=2, col="grey"); dev.off()
png(file.path(F,"bkmr_single.png"), width=1600, height=1000, res=180)
par(mfrow=c(1,3))
for(nm in c("Lead","Cadmium","Mercury")){
  r <- pr[pr[[grep("variable|^z$",names(pr),ignore.case=TRUE)[1]]]==nm,]
  if(!nrow(r)) next
  e <- r[[grep("est",names(r),ignore.case=TRUE)[1]]]; qv <- r[[grep("^q$|quantile",names(r),ignore.case=TRUE)[1]]]
  plot(e, qv, type="l", lwd=2, col="#c0392b", xlab="effect", ylab="quantile", main=nm); abline(v=0,lty=2,col="grey")
}
dev.off()
cat("BKMR DONE\n")
