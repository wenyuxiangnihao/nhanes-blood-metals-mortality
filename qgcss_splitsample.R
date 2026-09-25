.libPaths("/sandbox/workspace/Rlibs")
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu="adjust")
d <- readRDS("/sandbox/workspace/heavymetal/data/analysis_df.rds")
COVS <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
mkdes <- function(dd) svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=dd)
MET <- c("LBXBPB","LBXBCD","LBXTHG")
qg_scores <- function(dat, exps){
  for(m in exps){ br<-unique(quantile(dat[[m]],probs=seq(0,1,.25),na.rm=TRUE)); br[1]<--Inf; br[length(br)]<-Inf
    dat[[paste0("qc_",m)]] <- as.integer(cut(dat[[m]],breaks=br,labels=FALSE,include.lowest=TRUE))-1 }
  dat
}
fit_coef <- function(dat, exps, event, qc){
  f <- as.formula(paste0("Surv(time,",event,") ~ ",paste(qc,collapse="+")," + ",COVS))
  fit <- svycoxph(f, design=mkdes(dat))
  list(b=coef(fit)[qc], V=vcov(fit)[qc,qc], n=nrow(dat))
}
set.seed(20260925); B <- 200
out <- data.frame()
for(ev in c("event_all","event_cvd")){
  for(k in 1:B){
    idx <- sample(rep(1:2, length.out=nrow(d)))
    A <- d[idx==1,]; Bs <- d[idx==2,]
    A <- qg_scores(A, MET); Bs <- qg_scores(Bs, MET)
    qcA <- paste0("qc_",MET)
    rA <- tryCatch(fit_coef(A, MET, ev, qcA), error=function(e) NULL); if(is.null(rA)) next
    bA <- rA$b; names(bA) <- MET           # coefficients from half A
    signs <- sign(bA)                      # direction determined in half A (data-driven)
    rB <- tryCatch(fit_coef(Bs, MET, ev, qcA), error=function(e) NULL); if(is.null(rB)) next
    bB <- rB$b; VB <- rB$V
    pos <- MET[signs>0]; neg <- MET[signs<0]
    posq <- paste0("qc_",pos); negq <- paste0("qc_",neg)
    sump <- if(length(posq)) sum(bB[posq]) else 0
    sumn <- if(length(negq)) sum(bB[negq]) else 0
    sepp <- if(length(posq)) sqrt(sum(VB[posq,posq,drop=FALSE])) else NA
    senn <- if(length(negq)) sqrt(sum(VB[negq,negq,drop=FALSE])) else NA
    out <- rbind(out, data.frame(outcome=ev, rep=k, n_A=rA$n, n_B=rB$n,
      Hg_negative="LBXTHG"%in%neg, sign_pattern=paste(signs,collapse=""),
      HR_pos=exp(sump), pos_lo=exp(sump-1.96*sepp), pos_hi=exp(sump+1.96*sepp),
      HR_neg=exp(sumn), neg_lo=exp(sumn-1.96*senn), neg_hi=exp(sumn+1.96*senn)))
  }
}
write.csv(out, "/sandbox/workspace/heavymetal/results/rev7_qgcss.csv", row.names=FALSE)
cat("\n== sample-splitting (QGCSS, B=200): sign pattern frequency ==\n")
print(table(out$outcome, out$sign_pattern))
cat("\n== distribution of the partial effects across splits ==\n")
for(ev in c("event_all","event_cvd")){
  s <- out[out$outcome==ev,]
  q <- function(x) sprintf("%.2f (%.2f-%.2f)", median(x), quantile(x,.025), quantile(x,.975))
  cat(sprintf("  %-10s positive direction: %s   negative direction: %s\n", ev, q(s$HR_pos), q(s$HR_neg)))
}
