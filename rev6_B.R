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
# REV6 稳健性分析 B：qgcomp 界限值(assay floor)剔除 + 包-年敏感性
#                    + 缺失协变量比较(STROBE) + 多重比较(FDR)
# 输出：results/rev6_*.csv
# ============================================================
suppressMessages({library(survey); library(survival); library(foreign)})
options(survey.lonely.psu="adjust")
D <- file.path(ROOT, "data"); R <- file.path(ROOT, "results")

d <- readRDS(file.path(D,"analysis_df_v6.rds"))       # 分析集 n=33104
p <- readRDS(file.path(D,"pre_exclusion_df.rds"))     # 排除前 n=36743

MET3 <- c("LBXBPB","LBXBCD","LBXTHG")
# ln 与 z：与 A 脚本一致，在全分析集(33104)上标准化
for(m in MET3){ d[[paste0("ln_",m)]] <- log(d[[m]]) }
for(m in MET3){ v <- d[[paste0("ln_",m)]]; d[[paste0("z_",m)]] <- (v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE) }

COVS_M3    <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
COVS_M3_py <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + packyears + dm + htn + cvd"
mkdes <- function(dd) svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=dd)
safe <- function(expr,label){ tryCatch(expr, error=function(e){ cat("  !!! ERROR in",label,":",conditionMessage(e),"\n"); NULL }) }

# ---- qgcomp（严格照抄 A 脚本 (e) 段）----
qg <- function(dat, exps, event, covs=COVS_M3){
  dq <- dat
  for(m in exps){ br<-unique(quantile(dat[[m]],probs=seq(0,1,0.25),na.rm=TRUE)); br[1]<- -Inf; br[length(br)]<-Inf
    dq[[paste0("qc_",m)]] <- as.integer(cut(dat[[m]],breaks=br,labels=FALSE,include.lowest=TRUE))-1 }
  qc <- paste0("qc_",exps)
  desq <- svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~wt,nest=TRUE,data=dq)
  fF <- as.formula(paste0("Surv(time,",event,") ~ ",paste(qc,collapse="+")," + ",covs))
  ff <- svycoxph(fF,design=desq); b <- coef(ff)[qc]; V <- vcov(ff)[qc,qc]
  psi <- sum(b); se <- sqrt(sum(V))
  bp <- b[b>0]; bn <- b[b<0]
  sump <- if(length(bp)) sum(bp) else 0; sumn <- if(length(bn)) sum(bn) else 0
  Vpp <- if(length(bp)) sum(V[names(bp),names(bp),drop=FALSE]) else 0
  Vnn <- if(length(bn)) sum(V[names(bn),names(bn),drop=FALSE]) else 0
  sepp <- if(Vpp>0) sqrt(Vpp) else NA; senn <- if(Vnn>0) sqrt(Vnn) else NA
  list(b=b, psi=psi, se=se, sump=sump, sepp=sepp, sumn=sumn, senn=senn, n=nrow(dq))
}
p2 <- function(est,se) 2*pnorm(-abs(est/se))

cat("分析集 n =", nrow(d), " 全因死亡 =", sum(d$event_all), " CVD死亡 =", sum(d$event_cvd), "\n\n")

# ============================================================
# 任务 1：assay floor 剔除后的 3 金属 qgcomp（两种口径）
# ============================================================
cat("== 任务1: assay floor 剔除 + qgcomp ==\n")
lod <- read.csv(file.path(R,"lod_by_cycle_10cyc.csv"))
lod$yr <- as.integer(substr(lod$cycle,1,4))
gf <- function(yr,metal){ v<-lod$min_obs[lod$yr==yr & lod$metal==metal]; if(length(v)==0) NA_real_ else v[1] }
d$fl_PB <- sapply(d$cycle, gf, metal="LBXBPB")
d$fl_CD <- sapply(d$cycle, gf, metal="LBXBCD")
d$fl_HG <- sapply(d$cycle, gf, metal="LBXTHG")
d$at_PB <- abs(d$LBXBPB - d$fl_PB) < 1e-6
d$at_CD <- abs(d$LBXBCD - d$fl_CD) < 1e-6
d$at_HG <- abs(d$LBXTHG - d$fl_HG) < 1e-6
cat("  floor 命中数: Pb=",sum(d$at_PB)," Cd=",sum(d$at_CD)," Hg=",sum(d$at_HG),"\n",sep="")

scenarios <- list(
  "main(all obs)"        = rep(TRUE,nrow(d)),
  "excl Cd/Hg floor"     = !(d$at_CD | d$at_HG),
  "excl any metal floor" = !(d$at_PB | d$at_CD | d$at_HG)
)
s1 <- data.frame()
for(sn in names(scenarios)){
  dd <- d[scenarios[[sn]],]
  for(ev in c("event_all","event_cvd")){
    r <- safe(qg(dd, MET3, ev), paste0("floor ",sn,ev)); if(is.null(r)) next
    s1 <- rbind(s1, data.frame(
      scenario=sn, outcome=ev, n=r$n, n_deaths=sum(dd[[ev]]),
      psi=r$psi, HR_net=exp(r$psi), net_lo=exp(r$psi-1.96*r$se), net_hi=exp(r$psi+1.96*r$se),
      net_p=p2(r$psi,r$se),
      sum_beta_pos=r$sump, HR_pos=exp(r$sump), pos_lo=exp(r$sump-1.96*r$sepp), pos_hi=exp(r$sump+1.96*r$sepp), pos_p=p2(r$sump,r$sepp),
      sum_beta_neg=r$sumn, HR_neg=exp(r$sumn), neg_lo=exp(r$sumn-1.96*r$senn), neg_hi=exp(r$sumn+1.96*r$senn), neg_p=p2(r$sumn,r$senn)))
  }
}
write.csv(s1, file.path(R,"rev6_qgcomp_floor.csv"), row.names=FALSE)
print(s1, row.names=FALSE, digits=4)

# ============================================================
# 任务 2：包-年敏感性（仅 smoke∈{0,2}；不含前吸烟者）
# ============================================================
cat("\n== 任务2: 包-年敏感性 (smoke in {0,2}) ==\n")
dpy  <- d[d$smoke %in% c(0,2),]
dpy2 <- dpy[!is.na(dpy$packyears),]
cat("  smoke in {0,2}: n=",nrow(dpy)," | packyears 可算: n=",nrow(dpy2),"\n",sep="")

models <- list(
  "M3(smoke)"             = list(data=dpy,  covs=COVS_M3),
  "M3(smoke)[py-subset]"  = list(data=dpy2, covs=COVS_M3),
  "M3(packyears)"         = list(data=dpy2, covs=COVS_M3_py)
)
s2 <- data.frame()
for(mn in names(models)){
  dd <- models[[mn]]$data; cv <- models[[mn]]$covs
  for(ev in c("event_all","event_cvd")){
    for(m in MET3){
      vn <- paste0("z_",m)
      fit <- safe(svycoxph(as.formula(paste0("Surv(time,",ev,") ~ ",vn," + ",cv)), design=mkdes(dd)), paste0("py ",mn,m,ev))
      if(!is.null(fit)){ co<-summary(fit)$conf.int; cc<-summary(fit)$coefficients
        s2 <- rbind(s2, data.frame(model=mn, outcome=ev, term=vn, n=nrow(dd), n_deaths=sum(dd[[ev]]),
          HR=co[vn,1], lo=co[vn,3], hi=co[vn,4], p=cc[vn,ncol(cc)])) }
    }
    r <- safe(qg(dd, MET3, ev, covs=cv), paste0("py qg ",mn,ev)); if(is.null(r)) next
    s2 <- rbind(s2, data.frame(model=mn, outcome=ev, term="qgcomp_net", n=r$n, n_deaths=sum(dd[[ev]]),
        HR=exp(r$psi), lo=exp(r$psi-1.96*r$se), hi=exp(r$psi+1.96*r$se), p=p2(r$psi,r$se)))
    s2 <- rbind(s2, data.frame(model=mn, outcome=ev, term="qgcomp_pos", n=r$n, n_deaths=sum(dd[[ev]]),
        HR=exp(r$sump), lo=exp(r$sump-1.96*r$sepp), hi=exp(r$sump+1.96*r$sepp), p=p2(r$sump,r$sepp)))
    s2 <- rbind(s2, data.frame(model=mn, outcome=ev, term="qgcomp_neg", n=r$n, n_deaths=sum(dd[[ev]]),
        HR=exp(r$sumn), lo=exp(r$sumn-1.96*r$senn), hi=exp(r$sumn+1.96*r$senn), p=p2(r$sumn,r$senn)))
  }
}
write.csv(s2, file.path(R,"rev6_packyears.csv"), row.names=FALSE)
print(s2, row.names=FALSE, digits=4)

# ============================================================
# 任务 3：缺失协变量者 vs 纳入者 基线比较（STROBE）
# ============================================================
cat("\n== 任务3: 缺失 vs 纳入 基线比较 ==\n")
p$grp <- ifelse(p$included==TRUE,"included","excluded")
# 变量清单
cont_vars  <- c("RIDAGEYR","INDFMPIR","BMXBMI","LBXBPB","LBXBCD","LBXTHG","time")
cont_names <- c("年龄(岁)","家庭收入比(PIR)","BMI(kg/m2)","血铅(ug/dL)","血镉(ug/L)","血汞(ug/L)","随访(月)")
cat_vars   <- c("RIAGENDR","RIDRETH1","educ","smoke","dm","htn","cvd")
cat_names  <- c("性别","种族/族裔","教育程度","吸烟状态","糖尿病","高血压","基线CVD")

mkdes_p <- mkdes(p)   # 用于加权统计(全样本)
cnt <- function(x,grp){ x<-x[grp]; c(n=sum(!is.na(x)), m=mean(x,na.rm=TRUE), sd=sd(x,na.rm=TRUE)) }
# 未加权 & 加权
uw <- function(x,g) if(all(is.na(x[g]))) NA else mean(x[g],na.rm=TRUE)
uwsd <- function(x,g) if(length(x[g][!is.na(x[g])])<2) NA else sd(x[g],na.rm=TRUE)
wmean <- function(x,g){ sub<-p[g,]; ds<-mkdes(sub)
                        v<-svymean(as.formula(paste0("~",x)),ds,na.rm=TRUE)
                        s<-sqrt(as.numeric(svyvar(as.formula(paste0("~",x)),ds,na.rm=TRUE)))
                        c(m=as.numeric(v), sd=s) }
wprop <- function(x,lvl,g){ sub<-p[g,]; sub$zz <- as.numeric(as.character(sub[[x]])==as.character(lvl))
                        ds<-mkdes(sub); as.numeric(coef(svymean(~zz,ds,na.rm=TRUE))) }

gI <- p$grp=="included"; gE <- p$grp=="excluded"
out <- data.frame()
# 连续变量
for(i in seq_along(cont_vars)){
  v <- cont_vars[i]
  wi <- tryCatch(wmean(v,gI),error=function(e)c(m=NA,sd=NA)); we <- tryCatch(wmean(v,gE),error=function(e)c(m=NA,sd=NA))
  pv_uw <- tryCatch(t.test(p[[v]][gI],p[[v]][gE])$p.value, error=function(e) NA)
  pv_w  <- tryCatch({sub<-p; sub$g<-factor(gI); ds<-mkdes(sub)
                     svyttest(as.formula(paste0(v,"~g")), ds)$p.value}, error=function(e) NA)
  out <- rbind(out, data.frame(
    variable=v, varname=cont_names[i], level="", type="continuous",
    uw_incl=sprintf("%.2f (%.2f)",uw(p[[v]],gI),uwsd(p[[v]],gI)),
    uw_excl=sprintf("%.2f (%.2f)",uw(p[[v]],gE),uwsd(p[[v]],gE)),
    w_incl=sprintf("%.2f (%.2f)",wi["m"],wi["sd"]),
    w_excl=sprintf("%.2f (%.2f)",we["m"],we["sd"]), p_uw=pv_uw, p_w=pv_w))
}
# 分类变量
for(i in seq_along(cat_vars)){
  v <- cat_vars[i]
  levs <- sort(unique(na.omit(as.character(p[[v]]))))
  tab <- table(p[[v]], p$grp)
  pv_uw <- tryCatch(chisq.test(tab)$p.value, error=function(e) NA)
  pv_w  <- tryCatch({sub<-p; sub$g<-factor(gI); sub$.f<-factor(as.character(sub[[v]])); ds<-mkdes(sub)
                     svychisq(as.formula("~.f+g"), ds)$p.value}, error=function(e) NA)
  # header row
  out <- rbind(out, data.frame(variable=v, varname=cat_names[i], level="(分布%)", type="categorical",
    uw_incl="", uw_excl="", w_incl="", w_excl="", p_uw=pv_uw, p_w=pv_w))
  for(l in levs){
    nI <- sum(p[[v]][gI]==l,na.rm=TRUE); nE <- sum(p[[v]][gE]==l,na.rm=TRUE)
    uI <- 100*nI/sum(!is.na(p[[v]][gI])); uE <- 100*nE/sum(!is.na(p[[v]][gE]))
    wI <- 100*wprop(v,l,gI); wE <- 100*wprop(v,l,gE)
    out <- rbind(out, data.frame(variable=paste0("  ",v), varname="", level=l, type="cat-level",
      uw_incl=sprintf("%.1f%%",uI), uw_excl=sprintf("%.1f%%",uE),
      w_incl=sprintf("%.1f%%",wI), w_excl=sprintf("%.1f%%",wE), p_uw=NA, p_w=NA))
  }
}
# 全因死亡比例
for(ev in c("event_all")){
  nI<-sum(p[[ev]][gI]); nE<-sum(p[[ev]][gE])
  wI<-100*wprop(ev,1,gI); wE<-100*wprop(ev,1,gE)
  out <- rbind(out, data.frame(variable=ev, varname="全因死亡比例", level="death%", type="event",
    uw_incl=sprintf("%.1f%% (%d)",100*nI/sum(gI),nI), uw_excl=sprintf("%.1f%% (%d)",100*nE/sum(gE),nE),
    w_incl=sprintf("%.1f%%",wI), w_excl=sprintf("%.1f%%",wE),
    p_uw=tryCatch(chisq.test(table(p[[ev]],p$grp))$p.value,error=function(e)NA), p_w=NA))
}
# 组规模行
out <- rbind(data.frame(variable="N", varname="样本量", level="", type="size",
  uw_incl=as.character(sum(gI)), uw_excl=as.character(sum(gE)), w_incl="", w_excl="", p_uw=NA, p_w=NA), out)
write.csv(out, file.path(R,"rev6_missingness_compare.csv"), row.names=FALSE)
cat("  rev6_missingness_compare.csv 写出 (",nrow(out),"行 )\n",sep="")

# ---- 缺失频数表 ----
allcov <- c("RIDAGEYR","RIAGENDR","RIDRETH1","educ","INDFMPIR","BMXBMI","smoke","dm","htn","cvd")
freq <- data.frame()
wmiss <- function(col, idx=NULL){ dd <- if(is.null(idx)) p else p[idx,]
  dd$.na <- as.numeric(is.na(if(is.null(idx)) col else col[idx])); ds <- mkdes(dd)
  100*as.numeric(coef(svymean(~.na, ds))) }
for(v in c(allcov,"LBXBPB","LBXBCD","LBXTHG","time")){
  col <- p[[v]]
  nm_all <- sum(is.na(col)); nm_exc <- sum(is.na(col[gE])); nm_inc <- sum(is.na(col[gI]))
  wpct_all <- tryCatch(wmiss(col),error=function(e)NA)
  wpct_exc <- tryCatch(wmiss(col, which(gE)),error=function(e)NA)
  freq <- rbind(freq, data.frame(variable=v, n_missing_all=nm_all, pct_all=100*nm_all/nrow(p),
    n_missing_excl=nm_exc, pct_excl=100*nm_exc/sum(gE), n_missing_incl=nm_inc, pct_incl=100*nm_inc/sum(gI),
    w_pct_all=wpct_all, w_pct_excl=wpct_exc))
}
freq <- rbind(freq, data.frame(variable="ANY covariate missing", n_missing_all=sum(gE), pct_all=100*sum(gE)/nrow(p),
  n_missing_excl=sum(gE), pct_excl=100, n_missing_incl=0, pct_incl=0, w_pct_all=NA, w_pct_excl=100))
write.csv(freq, file.path(R,"rev6_missingness_freq.csv"), row.names=FALSE)
print(freq, row.names=FALSE, digits=3)

# ============================================================
# 任务 4：多重比较 FDR（BH）
# ============================================================
cat("\n== 任务4: FDR (BH) ==\n")
bh <- function(pv){ o<-order(pv); q<-pv[o]*length(pv)/seq_along(pv)
  q<-rev(cummin(rev(q))); out<-numeric(length(pv)); out[o]<-q; pmin(out,1) }
mc <- read.csv(file.path(R,"single_metal_cox_M1M2M3_10cyc.csv"))
m3 <- mc[mc$model=="M3" & mc$metal %in% MET3,]
m3$test <- paste0(m3$metal,"_",m3$outcome)
# 主检验 6 项
mainset <- m3[,c("test","HR","lo","hi","p")]
q6 <- bh(mainset$p)
main6 <- data.frame(set="main6", test=mainset$test, HR=mainset$HR, lo=mainset$lo, hi=mainset$hi, p=mainset$p, q_BH=q6)
# 扩展集 = 6 主检验 + qgcomp 净效应(2) + 方向分解(4) = 12 项
qg_ov3  <- subset(read.csv(file.path(R,"qgcomp_overall_10cyc.csv")), set=="3-metal")
qg_dir3 <- subset(read.csv(file.path(R,"qgcomp_directional_10cyc.csv")), set=="3-metal")
ext <- data.frame(test=character(),HR=numeric(),lo=numeric(),hi=numeric(),p=numeric())
for(i in 1:nrow(qg_ov3)) ext <- rbind(ext, data.frame(test=paste0("qgcomp_net_",qg_ov3$outcome[i]),
  HR=qg_ov3$HR_overall[i],lo=qg_ov3$lo[i],hi=qg_ov3$hi[i],p=qg_ov3$p[i]))
for(i in 1:nrow(qg_dir3)){
  sepos <- (log(qg_dir3$HR_pos[i])-log(qg_dir3$pos_lo[i]))/1.96
  seneg <- (log(qg_dir3$HR_neg[i])-log(qg_dir3$neg_lo[i]))/1.96
  ext <- rbind(ext, data.frame(test=paste0("qgcomp_pos_",qg_dir3$outcome[i]),HR=qg_dir3$HR_pos[i],lo=qg_dir3$pos_lo[i],hi=qg_dir3$pos_hi[i],p=p2(log(qg_dir3$HR_pos[i]),sepos)))
  ext <- rbind(ext, data.frame(test=paste0("qgcomp_neg_",qg_dir3$outcome[i]),HR=qg_dir3$HR_neg[i],lo=qg_dir3$neg_lo[i],hi=qg_dir3$neg_hi[i],p=p2(log(qg_dir3$HR_neg[i]),seneg)))
}
extall <- rbind(mainset, ext)
extq  <- bh(extall$p)
extended12 <- data.frame(set="extended12", test=extall$test, HR=extall$HR, lo=extall$lo, hi=extall$hi, p=extall$p, q_BH=extq)
fdr <- rbind(main6, extended12)
write.csv(fdr, file.path(R,"rev6_fdr_main.csv"), row.names=FALSE)
print(fdr, row.names=FALSE, digits=4)
cat("\n主检验 6 项 q<0.05 个数:", sum(q6<0.05), "/6\n")
cat("扩展集", nrow(extended12), "项 q<0.05 个数:", sum(extq<0.05), "/", nrow(extended12), "\n")
cat("\n== REV6_B 完成 ==\n")
