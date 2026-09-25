#!/usr/bin/env Rscript
# 探查各周期协变量文件的关键变量名（为清洗 coalesce 做准备）
.libPaths("/sandbox/workspace/Rlibs")
suppressMessages(library(foreign))
RAW <- "/sandbox/workspace/heavymetal/data_raw"

CYC <- list(
  list(yr=1999,suf="",  bio="LAB18",   chol="LAB13"),
  list(yr=2001,suf="_B",bio="L40_B",   chol="L13_B"),
  list(yr=2003,suf="_C",bio="L40_C",   chol="L13_C"),
  list(yr=2005,suf="_D",bio="BIOPRO_D",chol="TCHOL_D"),
  list(yr=2007,suf="_E",bio="BIOPRO_E",chol="TCHOL_E"),
  list(yr=2009,suf="_F",bio="BIOPRO_F",chol="TCHOL_F"),
  list(yr=2011,suf="_G",bio="BIOPRO_G",chol="TCHOL_G"),
  list(yr=2013,suf="_H",bio="BIOPRO_H",chol="TCHOL_H"),
  list(yr=2015,suf="_I",bio="BIOPRO_I",chol="TCHOL_I"),
  list(yr=2017,suf="_J",bio="BIOPRO_J",chol="TCHOL_J"))

rd <- function(fn) {
  p <- file.path(RAW, paste0(fn, ".XPT"))
  if (!file.exists(p)) return(NULL)
  tryCatch(read.xport(p, ) , error=function(e) NULL)
}
has <- function(x, v) if (is.null(x)) "NOFILE" else if (v %in% names(x)) "Y" else "-"

DEMO_V <- c("RIDAGEYR","RIAGENDR","RIDRETH1","RIDRETH3","DMDEDUC2","INDFMPIR",
            "WTMEC2YR","SDMVPSU","SDMVSTRA","RIDEXPRG")
for (c in CYC) {
  yr <- c$yr; s <- c$suf
  cat(sprintf("\n===== %s =====\n", yr))
  demo <- rd(paste0("DEMO", s))
  cat("DEMO:", paste(sprintf("%s=%s", DEMO_V, sapply(DEMO_V, function(v) has(demo,v))), collapse="  "), "\n")
  bmx <- rd(paste0("BMX",s)); bpx <- rd(paste0("BPX",s))
  cat("BMI:", has(bmx,"BMXBMI"), " WAIST:", has(bmx,"BMXWAIST"),
      " BPXSY1:", has(bpx,"BPXSY1"), " BPXDI1:", has(bpx,"BPXDI1"), "\n")
  smq <- rd(paste0("SMQ",s))
  cat("SMQ020:", has(smq,"SMQ020"), " SMQ040:", has(smq,"SMQ040"), "\n")
  diq <- rd(paste0("DIQ",s)); bpq <- rd(paste0("BPQ",s)); mcq <- rd(paste0("MCQ",s))
  cat("DIQ010:", has(diq,"DIQ010"), " BPQ020:", has(bpq,"BPQ020"), " BPQ050A:", has(bpq,"BPQ050A"),
      " MCQ160B:", has(mcq,"MCQ160B"), " MCQ160E:", has(mcq,"MCQ160E"), " MCQ160F:", has(mcq,"MCQ160F"), "\n")
  alq <- rd(paste0("ALQ",s))
  cat("ALQ101:", has(alq,"ALQ101"), " ALQ110:", has(alq,"ALQ110"), " ALQ120Q:", has(alq,"ALQ120Q"), "\n")
  bio <- rd(c$bio); chol <- rd(c$chol); hdl <- rd(paste0("HDL", s))
  cat("BIO[",c$bio,"] LBXSCR:", has(bio,"LBXSCR"), " LBDSCR:", has(bio,"LBDSCR"),
      " | CHOL[",c$chol,"] LBXTC:", has(chol,"LBXTC"), " LBDTCSI:", has(chol,"LBDTCSI"), "\n")
  cat("HDL[HDL",s,"]: LBDHDD:", has(hdl,"LBDHDD"), " LBXHDD:", has(hdl,"LBXHDD"),
      " LBDHDL:", if(!is.null(chol)) has(chol,"LBDHDL") else "-", "\n")
}
