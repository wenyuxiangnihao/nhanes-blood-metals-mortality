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
# 血重金属 × 死亡  数据清洗与合并
# 输出: data/analysis_df.rds   data/exclusion_counts.csv
# ============================================================
suppressMessages({library(foreign); library(survey); library(survival)})
RAW <- file.path(ROOT, "data_raw")
OUT <- file.path(ROOT, "data")
dir.create(OUT, showWarnings=FALSE)

CYC <- list(
  list(yr=1999,suf="",  pbcd="LAB06",  mort="NHANES_1999_2000_MORT_2019_PUBLIC.dat", bio="LAB18",   chol="LAB13", hdlf=NULL),
  list(yr=2001,suf="_B",pbcd="L06_B",  mort="NHANES_2001_2002_MORT_2019_PUBLIC.dat", bio="L40_B",   chol="L13_B", hdlf=NULL),
  list(yr=2003,suf="_C",pbcd="L06BMT_C", mort="NHANES_2003_2004_MORT_2019_PUBLIC.dat", bio="L40_C", chol="L13_C", hdlf=NULL),
  list(yr=2005,suf="_D",pbcd="PBCD_D", mort="NHANES_2005_2006_MORT_2019_PUBLIC.dat", bio="BIOPRO_D",chol="TCHOL_D",hdlf="HDL_D"),
  list(yr=2007,suf="_E",pbcd="PBCD_E", mort="NHANES_2007_2008_MORT_2019_PUBLIC.dat", bio="BIOPRO_E",chol="TCHOL_E",hdlf="HDL_E"),
  list(yr=2009,suf="_F",pbcd="PBCD_F", mort="NHANES_2009_2010_MORT_2019_PUBLIC.dat", bio="BIOPRO_F",chol="TCHOL_F",hdlf="HDL_F"),
  list(yr=2011,suf="_G",pbcd="PBCD_G", mort="NHANES_2011_2012_MORT_2019_PUBLIC.dat", bio="BIOPRO_G",chol="TCHOL_G",hdlf="HDL_G"),
  list(yr=2013,suf="_H",pbcd="PBCD_H", mort="NHANES_2013_2014_MORT_2019_PUBLIC.dat", bio="BIOPRO_H",chol="TCHOL_H",hdlf="HDL_H"),
  list(yr=2015,suf="_I",pbcd="PBCD_I", mort="NHANES_2015_2016_MORT_2019_PUBLIC.dat", bio="BIOPRO_I",chol="TCHOL_I",hdlf="HDL_I"),
  list(yr=2017,suf="_J",pbcd="PBCD_J", mort="NHANES_2017_2018_MORT_2019_PUBLIC.dat", bio="BIOPRO_J",chol="TCHOL_J",hdlf="HDL_J"))
NCYC <- length(CYC)

rd <- function(fn){ p<-file.path(RAW,paste0(fn,".XPT")); if(!file.exists(p)) return(NULL)
  tryCatch(read.xport(p), error=function(e) NULL) }
pick <- function(df, cands){ if(is.null(df)) return(rep(NA_real_,0))
  for(v in cands) if(v %in% names(df)) return(df[[v]]); rep(NA_real_, nrow(df)) }
byname <- function(df, cols){ if(is.null(df)) return(NULL)
  keep <- intersect(cols, names(df)); df[, keep, drop=FALSE] }

read_mort <- function(f){
  lines <- readLines(file.path(RAW,f), warn=FALSE)
  lines <- lines[nchar(lines)>0]
  g <- function(a,b){ x<-trimws(substr(lines,a,b)); x[x=="."|x==""]<-NA; as.integer(x) }
  data.frame(SEQN=g(1,6), eligstat=g(15,15), mortstat=g(16,16),
             ucod_leading=g(17,19), permth_int=g(43,45), permth_exm=g(46,48))
}

frames <- list()
for(c in CYC){
  pb <- rd(c$pbcd); if(is.null(pb)) { message("MISSING pbcd ", c$pbcd); next }
  mcols <- intersect(c("SEQN","LBXBPB","LBXBCD","LBXTHG","LBXBSE","LBXBMN"), names(pb))
  pb <- pb[, mcols, drop=FALSE]
  for(v in c("LBXBPB","LBXBCD","LBXTHG","LBXBSE","LBXBMN")) if(!v %in% names(pb)) pb[[v]] <- NA_real_
  pb <- pb[, c("SEQN","LBXBPB","LBXBCD","LBXTHG","LBXBSE","LBXBMN")]

  demo <- rd(paste0("DEMO",c$suf))
  demo <- byname(demo, c("SEQN","RIDAGEYR","RIAGENDR","RIDRETH1","RIDRETH3","DMDEDUC2",
                         "INDFMPIR","WTMEC2YR","SDMVPSU","SDMVSTRA","RIDEXPRG"))
  bmx <- byname(rd(paste0("BMX",c$suf)), c("SEQN","BMXBMI","BMXWAIST"))
  bpx <- byname(rd(paste0("BPX",c$suf)), c("SEQN","BPXSY1","BPXSY2","BPXSY3","BPXSY4","BPXDI1","BPXDI2","BPXDI3","BPXDI4"))
  smq <- byname(rd(paste0("SMQ",c$suf)), c("SEQN","SMQ020","SMQ040","SMD030","SMD055","SMD070","SMD641","SMD650"))
  diq <- byname(rd(paste0("DIQ",c$suf)), c("SEQN","DIQ010"))
  bpq <- byname(rd(paste0("BPQ",c$suf)), c("SEQN","BPQ020","BPQ050A"))
  mcq <- byname(rd(paste0("MCQ",c$suf)), c("SEQN","MCQ160B","MCQ160C","MCQ160E","MCQ160F","MCQ220"))
  alq <- byname(rd(paste0("ALQ",c$suf)), c("SEQN","ALQ101","ALQ110","ALQ120Q","ALQ130"))
  bio <- byname(rd(c$bio), c("SEQN","LBXSCR","LBDSCR"))
  chol <- byname(rd(c$chol), c("SEQN","LBXTC","LBDTCSI","LBDHDL","LBXHDD"))
  hdl <- if(!is.null(c$hdlf)) byname(rd(c$hdlf), c("SEQN","LBDHDD","LBXHDD")) else NULL
  cot <- byname(rd(paste0("COT",c$suf)), c("SEQN","LBXCOT","LBDCOTLC"))
  if(is.null(cot) && c$yr %in% c(1999,2001)) cot <- byname(pb0 <- pb, c())  # handled below

  d <- pb
  for(part in list(demo,bmx,bpx,smq,diq,bpq,mcq,alq,bio,chol,hdl)) if(!is.null(part)) d <- merge(d, part, by="SEQN", all.x=TRUE)
  # cotinine from blood metal file if present in original pbcd (LAB06/L06_B)
  if(c$yr %in% c(1999,2001)) { c2 <- byname(rd(c$pbcd), c("SEQN","LBXCOT","LBDCOTLC")); if(!is.null(c2)) d <- merge(d, c2, by="SEQN", all.x=TRUE) }
  else if(!is.null(cot)) d <- merge(d, cot, by="SEQN", all.x=TRUE)

  m <- read_mort(c$mort)
  d <- merge(d, m, by="SEQN", all.x=TRUE)
  d$cycle <- c$yr
  frames[[as.character(c$yr)]] <- d
  cat(sprintf("cycle %d: pbcd n=%d -> merged n=%d\n", c$yr, nrow(pb), nrow(d)))
}

df <- do.call(rbind, lapply(frames, function(x){ for(cc in unique(unlist(lapply(frames,names)))) if(!cc %in% names(x)) x[[cc]]<-NA; x }))
cat("pooled rows:", nrow(df), "\n")

# ---------- derived variables ----------
mean_pos <- function(d, pat){
  cols <- grep(pat, names(d), value=TRUE)
  if(!length(cols)) return(rep(NA_real_, nrow(d)))
  M <- as.matrix(d[,cols,drop=FALSE]); M[M<=0] <- NA
  rowMeans(M, na.rm=TRUE)
}
df$sbp <- mean_pos(df, "^BPXSY[1-4]$")
df$dbp <- mean_pos(df, "^BPXDI[1-4]$")
df$sbp[is.nan(df$sbp)] <- NA; df$dbp[is.nan(df$dbp)] <- NA

# creatinine coalesce
df$creat <- ifelse(!is.na(df$LBXSCR), df$LBXSCR, df$LBDSCR)
# HDL coalesce (chol file or HDL file)
df$hdl <- ifelse(!is.na(df$LBDHDD), df$LBDHDD, ifelse(!is.na(df$LBXHDD), df$LBXHDD, df$LBDHDL))
df$tc  <- df$LBXTC
# cotinine coalesce
df$cot <- ifelse(!is.na(df$LBXCOT), df$LBXCOT, NA)

# education / missing codes
df$educ <- df$DMDEDUC2; df$educ[df$educ %in% c(7,9)] <- NA
df$smq020 <- df$SMQ020; df$smq020[df$smq020 %in% c(7,9)] <- NA
df$smq040 <- df$SMQ040; df$smq040[df$smq040 %in% c(7,9)] <- NA
df$diq010 <- df$DIQ010; df$diq010[df$diq010 %in% c(7,9)] <- NA
df$bpq020 <- df$BPQ020; df$bpq020[df$bpq020 %in% c(7,9)] <- NA
df$bpq050a <- df$BPQ050A; df$bpq050a[df$bpq050a %in% c(7,9)] <- NA

# eGFR CKD-EPI 2021 (race-free)
scr <- df$creat
k <- ifelse(df$RIAGENDR==2, 0.7, 0.9)
a <- ifelse(df$RIAGENDR==2, -0.241, -0.302)
df$egfr <- 142 * pmin(scr/k,1)^a * pmax(scr/k,1)^(-1.200) * 0.9938^df$RIDAGEYR * ifelse(df$RIAGENDR==2, 1.012, 1)

# smoking 3-level
df$smoke <- NA_integer_
df$smoke[df$smq020==2] <- 0                      # <100 cigs -> never
df$smoke[df$smq020==1 & df$smq040==3] <- 1       # former
df$smoke[df$smq020==1 & df$smq040 %in% c(1,2)] <- 2  # current
# diabetes / hypertension / cvd
# diabetes / hypertension / cvd  (robust: any TRUE -> 1; NA only if ALL sources NA)
df$dm  <- as.integer(df$diq010==1); df$dm[is.na(df$diq010)] <- NA
src <- cbind(df$bpq020==1, df$bpq050a==1,
             !is.na(df$sbp) & df$sbp>=140, !is.na(df$dbp) & df$dbp>=90)
df$htn <- as.integer(rowSums(src, na.rm=TRUE) > 0)
df$htn[is.na(df$bpq020) & is.na(df$sbp) & is.na(df$dbp)] <- NA
mcq_raw <- cbind(df$MCQ160B, df$MCQ160C, df$MCQ160E, df$MCQ160F)
mcq_raw[mcq_raw %in% c(7, 9)] <- NA          # refused / don't know -> missing
df$cvd <- ifelse(rowSums(mcq_raw == 1, na.rm = TRUE) > 0, 1L,
                 ifelse(rowSums(mcq_raw == 2, na.rm = TRUE) > 0, 0L, NA))

# ---- mortality outcome ----
df$cvd_death <- as.integer(!is.na(df$mortstat) & df$mortstat==1 & !is.na(df$ucod_leading) & df$ucod_leading %in% c(1,5))
df$time <- df$permth_exm
df$time[is.na(df$time)] <- df$permth_int[is.na(df$time)]
df$event_all <- as.integer(df$mortstat==1)
df$event_cvd <- df$cvd_death

# ---------- exclusion chain ----------
# v7 (Table S24): keep the full merged frame, so that participants who were
# excluded before the analytic sample (in particular those with no usable blood
# specimen, n = 8,640) can be compared with those included.
saveRDS(df, file.path(OUT, "all_merged_df.rds"))
ec <- data.frame(step=character(), n=integer())
allcyc <- sort(unique(df$cycle))
add <- function(lbl, d){
  ec <<- rbind(ec, data.frame(step=lbl, n=nrow(d)))
  cat(sprintf("  [%-38s] %s\n", lbl,
      if (nrow(d) && "cycle" %in% names(d)) paste(sprintf("%6d", table(factor(d$cycle, levels=allcyc))), collapse="") else ""))
  d }
d0 <- df
d1 <- add("0. 合并血金属+死亡+协变量", d0)
d2 <- add("1. 年龄>=20岁", d1[!is.na(d1$RIDAGEYR) & d1$RIDAGEYR>=20, ])
d3 <- add("2. 非孕妇", d2[is.na(d2$RIDEXPRG) | d2$RIDEXPRG!=1, ])
d4 <- add("3. 死亡随访合格(eligstat=1)", d3[!is.na(d3$eligstat) & d3$eligstat==1, ])
d5 <- add("4. 三种核心金属(Pb/Cd/Hg)齐全", d4[!(is.na(d4$LBXBPB)|is.na(d4$LBXBCD)|is.na(d4$LBXTHG)), ])
d6 <- add("5. 有生存时间(time)", d5[!is.na(d5$time) & d5$time>0, ])
cat("\n>> per-cycle counts / non-missing at step 5:\n")
cs <- sort(unique(d6$cycle))
cat(sprintf("%-10s %s\n", "N", paste(sprintf("%6d", table(factor(d6$cycle, levels=cs))), collapse="")))
for (v in c("RIDAGEYR","educ","INDFMPIR","BMXBMI","smoke","dm","htn","cvd","MCQ160B")) {
  tb <- tapply(!is.na(d6[[v]]), factor(d6$cycle, levels=cs), sum)
  cat(sprintf("%-10s %s\n", v, paste(sprintf("%6d", tb), collapse="")))
}
cat("\n")
cat("\n>> Missingness at step 5 (n=", nrow(d6), "):\n", sep="")
for(v in c("RIDAGEYR","RIAGENDR","RIDRETH1","educ","INDFMPIR","BMXBMI","smoke","dm","htn","cvd"))
  cat(sprintf("   %-10s NA=%6d (%.1f%%)\n", v, sum(is.na(d6[[v]])), 100*mean(is.na(d6[[v]]))))
cat("   complete-case n =", sum(complete.cases(d6[, c("RIDAGEYR","RIAGENDR","RIDRETH1","educ","INDFMPIR","BMXBMI","smoke","dm","htn","cvd")])), "\n")
for(v in c("bpq020","bpq050a","sbp","dbp","MCQ160B","MCQ160C","MCQ160E","MCQ160F","ALQ110","ALQ101","cot")) 
  cat(sprintf("   [raw] %-9s NA=%6d (%.1f%%)\n", v, sum(is.na(d6[[v]])), 100*mean(is.na(d6[[v]]))))
d7 <- add("6. 协变量完整", d6[complete.cases(d6[, c("RIDAGEYR","RIAGENDR","RIDRETH1","educ","INDFMPIR","BMXBMI","smoke","dm","htn","cvd")]), ])
write.csv(ec, file.path(OUT,"exclusion_counts.csv"), row.names=FALSE)
print(ec)

# ---- v6 additions: baseline cancer + smoking pack-years ----
d7$cancer_ever <- ifelse(d7$MCQ220==1, 1L, ifelse(d7$MCQ220==2, 0L, NA_integer_))
d7$cancer_ever[d7$MCQ220 %in% c(7,9)] <- NA
qd <- function(x, bad=c(7,77,777,7777,9,99,999,9999)){ x[x %in% bad] <- NA; x }
d7$cigs_day <- ifelse(!is.na(qd(d7$SMD070)), qd(d7$SMD070), qd(d7$SMD650))
d7$age_start <- qd(d7$SMD030); d7$age_quit <- qd(d7$SMD055)
end_age <- ifelse(d7$smoke==2, d7$RIDAGEYR, ifelse(d7$smoke==1, d7$age_quit, 0))
yrs <- end_age - d7$age_start
d7$packyears <- (d7$cigs_day/20) * yrs
d7$packyears[d7$smoke==0] <- 0
d7$packyears[!is.na(d7$packyears) & (d7$packyears<0 | d7$packyears>200)] <- NA
cat("\n[v6] cancer_ever: ", paste(names(table(d7$cancer_ever,useNA="ifany")), table(d7$cancer_ever,useNA="ifany"), collapse=" "), "\n")
cat("[v6] packyears non-missing:", sum(!is.na(d7$packyears)), " among ever-smokers:", sum(!is.na(d7$packyears[d7$smoke %in% c(1,2)])), "/", sum(d7$smoke %in% c(1,2)), "\n")

# weights
d7$wt <- d7$WTMEC2YR / NCYC

saveRDS(d7, file.path(OUT,"analysis_df_v6.rds"))
cc <- c("RIDAGEYR","RIAGENDR","RIDRETH1","educ","INDFMPIR","BMXBMI","smoke","dm","htn","cvd")
d6$included <- complete.cases(d6[, cc])
d6$wt <- d6$WTMEC2YR / NCYC
saveRDS(d6, file.path(OUT,"pre_exclusion_df.rds"))
cat("pre-exclusion frame (step 5) n =", nrow(d6), " included =", sum(d6$included), "\n")
saveRDS(d7, file.path(OUT,"analysis_df.rds"))
cat("\nFINAL n =", nrow(d7), " all-deaths =", sum(d7$event_all), " cvd-deaths =", sum(d7$event_cvd), "\n")
cat("cycles:", paste(sort(unique(d7$cycle)), collapse=","), "\n")
cat("Se present:", sum(!is.na(d7$LBXBSE)), " Mn present:", sum(!is.na(d7$LBXBMN)), "\n")
