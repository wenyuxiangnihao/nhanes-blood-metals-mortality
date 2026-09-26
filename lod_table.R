# --- project root: $NHANES_ROOT, or the folder containing this script --------
if (!nzchar(Sys.getenv("NHANES_ROOT"))) {
  .ca <- commandArgs(trailingOnly = FALSE)
  .cf <- sub("^--file=", "", .ca[grep("^--file=", .ca)])
  Sys.setenv(NHANES_ROOT = if (length(.cf)) dirname(normalizePath(.cf)) else getwd())
}
ROOT <- Sys.getenv("NHANES_ROOT")
if (dir.exists(file.path(dirname(ROOT), "Rlibs")))
  .libPaths(c(file.path(dirname(ROOT), "Rlibs"), .libPaths()))
suppressMessages(library(foreign))
RAW <- file.path(ROOT, "data_raw"); OUT <- file.path(ROOT, "results")
CYC <- list(list(yr="1999-2000",f="LAB06"),list(yr="2001-2002",f="L06_B"),list(yr="2003-2004",f="L06BMT_C"),
            list(yr="2005-2006",f="PBCD_D"),list(yr="2007-2008",f="PBCD_E"),list(yr="2009-2010",f="PBCD_F"),
            list(yr="2011-2012",f="PBCD_G"),list(yr="2013-2014",f="PBCD_H"),list(yr="2015-2016",f="PBCD_I"),
            list(yr="2017-2018",f="PBCD_J"))
MET <- c("LBXBPB","LBXBCD","LBXTHG")
out <- data.frame()
for(c in CYC){
  p <- file.path(RAW, paste0(c$f,".XPT")); if(!file.exists(p)) { cat("MISSING",c$f,"\n"); next }
  d <- read.xport(p)
  for(m in MET){
    if(!m %in% names(d)) next
    v <- d[[m]]; v <- v[!is.na(v)]; n <- length(v)
    # round to reported precision (2 decimals for these assays)
    key <- sprintf("%.2f", v)
    tb <- sort(table(key), decreasing=TRUE)
    vmin <- min(v); fmin <- sum(v==vmin)
    out <- rbind(out, data.frame(cycle=c$yr, file=c$f, metal=m, n=n,
      n_at_min=fmin, pct_below_LOD=round(100*fmin/n,2),
      min_obs=vmin, inferred_LOD=round(vmin*sqrt(2),3)))
  }
}
write.csv(out, file.path(OUT,"lod_by_cycle_10cyc.csv"), row.names=FALSE)
print(out, row.names=FALSE)
