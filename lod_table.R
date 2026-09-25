.libPaths("/sandbox/workspace/Rlibs")
suppressMessages(library(foreign))
RAW <- "/sandbox/workspace/heavymetal/data_raw"; OUT <- "/sandbox/workspace/heavymetal/results"
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
