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
# Redraw BKMR summary figures from results/bkmr_*_summaries.csv (no BKMR fit needed).
# Content is identical to the original figures/bkmr_*.png; only export formats/resolution added.
R   <- file.path(ROOT, "results")
F   <- file.path(ROOT, "figures")
SUB <- file.path(ROOT, "figures_submission")
dir.create(F, showWarnings = FALSE); dir.create(SUB, showWarnings = FALSE)
ov <- read.csv(file.path(R,"bkmr_overall_summaries.csv"))
pr <- read.csv(file.path(R,"bkmr_single_summaries.csv"))

plot_overall <- function(){
  lo <- ov$est-1.96*ov$sd; hi <- ov$est+1.96*ov$sd
  plot(ov$est, ov$quantile, type="b", pch=19, xlim=range(c(lo,hi,0)), ylim=c(0,1),
       xlab="effect (probit scale, vs all metals at median)", ylab="Quantile of all metals",
       main="BKMR overall mixture effect (all-cause)")
  segments(lo, ov$quantile, hi, ov$quantile); abline(v=0, lty=2, col="grey")
}

plot_single <- function(){
  par(mfrow=c(1,3), mar=c(4,4,3,1))
  for(nm in c("Lead","Cadmium","Mercury")){
    r <- pr[pr$variable==nm,]; if(!nrow(r)) next
    r <- r[order(r$z),]
    plot(r$z, r$est, type="l", lwd=2, col="#c0392b", xlab="metal (quantile-scaled z)", ylab="effect (probit)",
         main=nm, ylim=range(c(r$est-1.96*r$se, r$est+1.96*r$se, 0), finite=TRUE))
    polygon(c(r$z, rev(r$z)), c(r$est-1.96*r$se, rev(r$est+1.96*r$se)), col=rgb(0.75,0.22,0.17,0.18), border=NA)
    abline(h=0, lty=2, col="grey")
  }
}

# ---- legacy 300-dpi-class PNG preview (unchanged content) ----
png(file.path(F,"bkmr_overall.png"), width=1200, height=1000, res=170); plot_overall(); dev.off()
png(file.path(F,"bkmr_single.png"),  width=1600, height=560,  res=150); plot_single();  dev.off()

# ---- submission exports: PNG(600 dpi) + TIFF(600 dpi) + vector EPS ----
# physical size held constant (= legacy preview) while pixels are scaled to a true 600 dpi
png(file.path(SUB,"bkmr_overall.png"),  width=4235, height=3529, res=600); plot_overall(); dev.off()
tiff(file.path(SUB,"bkmr_overall.tiff"), width=4235, height=3529, res=600, compression="lzw"); plot_overall(); dev.off()
setEPS(); postscript(file.path(SUB,"bkmr_overall.eps"), width=1200/170, height=1000/170); plot_overall(); dev.off()

png(file.path(SUB,"bkmr_single.png"),  width=6400, height=2240, res=600); plot_single(); dev.off()
tiff(file.path(SUB,"bkmr_single.tiff"), width=6400, height=2240, res=600, compression="lzw"); plot_single(); dev.off()
setEPS(); postscript(file.path(SUB,"bkmr_single.eps"), width=1600/150, height=560/150); plot_single(); dev.off()

cat("bkmr plots regenerated\n")
