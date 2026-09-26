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
# REV6 Part A -- 稳健性分析A：分解/方向稳定性、Fine-Gray、排除基线癌症、
#               逐步调整、铅四分位人年
# 数据: data/analysis_df_v6.rds  (10 周期, n=33104, 全因 4260, CVD 1332)
# 输出: results/rev6_*.csv
# 严禁修改 data/ 与 results/ 下任何既有文件
# ============================================================
suppressMessages({library(survey); library(survival)})
options(survey.lonely.psu = "adjust")
D <- file.path(ROOT, "data")
R <- file.path(ROOT, "results")

d <- readRDS(file.path(D, "analysis_df_v6.rds"))
MET3 <- c("LBXBPB", "LBXBCD", "LBXTHG")
MET5 <- c(MET3, "LBXBSE", "LBXBMN")
COVS_M3 <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"

# 与 analysis_reanalysis_A.R 一致的 ln/z 转换
for (m in MET5) d[[paste0("ln_", m)]] <- log(d[[m]])
for (m in MET5) { v <- d[[paste0("ln_", m)]]; d[[paste0("z_", m)]] <- (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE) }

mkdes <- function(dd) svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = dd)
getc <- function(fit, term) { capture.output(sm <- summary(fit))  # survey 会在 summary 时自动打印设计, 静默之
  co <- sm$coefficients
  c(HR = sm$conf.int[term, 1], lo = sm$conf.int[term, 3], hi = sm$conf.int[term, 4], p = co[term, ncol(co)]) }
safe <- function(e, l) tryCatch(e, error = function(x) { cat("  !!! ERROR", l, ":", conditionMessage(x), "\n"); NULL })

cat("=== REV6-A START ===\n")
cat("analytic n =", nrow(d), " all-deaths =", sum(d$event_all, na.rm = TRUE),
    " cvd =", sum(d$event_cvd, na.rm = TRUE), "\n")

# 核对主模型基准 --------------------------------------------------------------
cat("\n=== 主模型基准核对 (svycoxph M3) ===\n")
base_chk <- data.frame()
for (ev in c("event_all", "event_cvd")) {
  des <- mkdes(d)
  for (m in MET3) {
    fit <- safe(svycoxph(as.formula(paste0("Surv(time,", ev, ") ~ z_", m, " + ", COVS_M3)), design = des), paste(ev, m))
    if (is.null(fit)) next
    r <- getc(fit, paste0("z_", m))
    base_chk <- rbind(base_chk, data.frame(outcome = ev, metal = m, HR = r["HR"], lo = r["lo"], hi = r["hi"], p = r["p"]))
  }
}
print(base_chk, row.names = FALSE)

# ##### 任务1：方向归属稳定性 ###################################################
cat("\n\n=== 任务1: 方向归属稳定性 (复权重拟合) ===\n")

# 四分位得分 qc_*，与 analysis_reanalysis_A.R 的 qg() 完全一致
score_q <- function(dat, exps) {
  dq <- dat
  for (m in exps) {
    br <- unique(quantile(dat[[m]], probs = seq(0, 1, 0.25), na.rm = TRUE))
    br[1] <- -Inf; br[length(br)] <- Inf
    dq[[paste0("qc_", m)]] <- as.integer(cut(dat[[m]], breaks = br, labels = FALSE, include.lowest = TRUE)) - 1
  }
  dq
}
dq <- score_q(d, MET5)

des_main <- mkdes(dq)
# 分层删除-one-PSU 刀切 (JKn -> 301 个复权点 = 301 PSU)。JK1 不支持分层设计。
repdes <- as.svrepdesign(des_main, type = "JKn")
RW <- weights(repdes, type = "analysis")   # n x 301
NREP <- ncol(RW)
cat("复权方案: as.svrepdesign(type='JKn')(=分层删除一个PSU刀切, JK1不支持分层设计) 复权点数 =", NREP, "\n")

# PSU 整群自助 (备用方案, 按 stratum 内有放回抽 PSU)
psu_tab <- unique(dq[, c("SDMVSTRA", "SDMVPSU")])
boot_coefs <- function(f, terms, B = 200) {
  out <- matrix(NA_real_, B, length(terms), dimnames = list(NULL, terms))
  strat <- as.character(psu_tab$SDMVSTRA)
  for (b in 1:B) {
    sel <- unlist(lapply(split(psu_tab$SDMVPSU, strat), function(p) sample(p, length(p), replace = TRUE)))
    # 重建抽样后的数据 (按 stratum-psu 组合)
    key_sel <- paste(rep(names(split(psu_tab$SDMVPSU, strat)), lengths(split(psu_tab$SDMVPSU, strat))), sel, sep = "_")
    dq$.key <- paste(dq$SDMVSTRA, dq$SDMVPSU, sep = "_")
    idx <- match(key_sel, dq$.key)
    idx <- idx[!is.na(idx)]
    dd <- dq[idx, ]
    fit <- safe(coxph(f, data = dd, weights = dd$wt), paste0("boot", b))
    if (!is.null(fit)) out[b, ] <- coef(fit)[terms]
  }
  out
}

run_stability <- function(exps, event) {
  qc <- paste0("qc_", exps)
  f <- as.formula(paste0("Surv(time,", event, ") ~ ", paste(qc, collapse = "+"), " + ", COVS_M3))
  # 全样本参考系数 (qg 使用的 svycoxph 口径)
  ffull <- safe(svycoxph(f, design = mkdes(dq)), "full")
  bfull <- coef(ffull)[qc]
  # 计时前 5 个复权点 (删除的 PSU 复权权重=0 -> coxph 不接受0权重, 需剔除该 PSU 行)
  t0 <- Sys.time(); ntest <- min(5, NREP)
  for (r in 1:ntest) { wr <- RW[, r]; keep <- wr > 0; safe(coxph(f, data = dq[keep, ], weights = wr[keep]), "test") }
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / ntest
  est <- dt * NREP / 60
  cat(sprintf("  [%s/%s] 单次 %.3fs, 外推 %d 复权点 = %.1f 分钟\n", paste(exps, collapse = "+"), event, dt, NREP, est))
  if (est <= 40) {
    method <- sprintf("svyrep_JKn(%d)", NREP)
    M <- matrix(NA_real_, NREP, length(qc), dimnames = list(NULL, qc))
    for (r in 1:NREP) {
      wr <- RW[, r]; keep <- wr > 0
      fit <- safe(coxph(f, data = dq[keep, ], weights = wr[keep]), paste0("rep", r))
      if (!is.null(fit)) M[r, ] <- coef(fit)[qc]
    }
  } else {
    B <- 200; method <- sprintf("psu_bootstrap(B=%d)", B)
    cat("  -> 超过40分钟, 降级为 PSU 整群自助 B=", B, "\n")
    M <- boot_coefs(f, qc, B)
  }
  res <- data.frame()
  for (j in seq_along(qc)) {
    b <- M[, j]; b <- b[is.finite(b)]
    res <- rbind(res, data.frame(
      mixture = if (length(exps) == 3) "3-metal" else "5-metal",
      metal = sub("qc_", "", qc[j]), outcome = event, method = method,
      n_reps = length(b), b_full = bfull[j],
      prop_b_pos = mean(b > 0), prop_sign_consistent = mean(sign(b) == sign(bfull[j])),
      mean_b = mean(b), sd_b = sd(b),
      q025 = as.numeric(quantile(b, .025)), q975 = as.numeric(quantile(b, .975))))
  }
  res
}

stab <- data.frame()
for (exps in list(MET3, MET5)) {
  for (ev in c("event_all", "event_cvd")) {
    r <- safe(run_stability(exps, ev), "stab")
    if (!is.null(r)) stab <- rbind(stab, r)
  }
}
write.csv(stab, file.path(R, "rev6_direction_stability.csv"), row.names = FALSE)
cat("\n-- 方向稳定性结果 --\n"); print(stab, row.names = FALSE)

# ##### 任务2：Fine-Gray 竞争风险 ################################################
cat("\n\n=== 任务2: Fine-Gray 竞争风险 (未加权) ===\n")
if (!requireNamespace("cmprsk", quietly = TRUE)) {
  try(install.packages("cmprsk", repos = "https://mirrors.aliyun.com/CRAN/", quiet = TRUE), silent = TRUE)
}
cat("cmprsk 可用:", requireNamespace("cmprsk", quietly = TRUE), "\n")

d$event_cancer <- as.integer(d$mortstat == 1 & !is.na(d$ucod_leading) & d$ucod_leading == 2)
d$event_cerebro <- as.integer(d$mortstat == 1 & !is.na(d$ucod_leading) & d$ucod_leading == 5)

causes <- list(CVD_death = c(1, 5), cancer_death = 2, cerebrovascular_death = 5, CLRD_death = 3)
# 已加权的 cause-specific Cox (用于对比)
cs <- read.csv(file.path(R, "cause_specific_10cyc.csv"), stringsAsFactors = FALSE)
csmap <- c("CVD_death" = "heart_disease", "cancer_death" = "cancer",
           "cerebrovascular_death" = "cerebrovascular", "CLRD_death" = "CLRD")

# 预先命名协变量列(finegray 会字面保留因子项名, 需用简单列名)
d$age <- d$RIDAGEYR; d$sex <- factor(d$RIAGENDR); d$race <- factor(d$RIDRETH1)
d$edu_f <- factor(d$educ); d$income <- d$INDFMPIR; d$bmi <- d$BMXBMI
d$smokef <- factor(d$smoke); d$dmz <- d$dm; d$htnz <- d$htn; d$cvdz <- d$cvd
COVS_FG <- "age + sex + race + edu_f + income + bmi + smokef + dmz + htnz + cvdz"

fg_out <- data.frame()
for (cn in names(causes)) {
  dd <- d
  dd$status <- factor(ifelse(dd$mortstat == 0, "censor",
                             ifelse(!is.na(dd$ucod_leading) & dd$ucod_leading %in% causes[[cn]], "event", "competing")),
                      levels = c("censor", "event", "competing"))
  dd <- dd[!is.na(dd$time) & !is.na(dd$status), ]
  n_event <- sum(dd$status == "event")
  for (m in MET3) {
    dd$zm <- dd[[paste0("z_", m)]]
    fgdat <- safe(finegray(as.formula(paste0("Surv(time,status) ~ zm + ", COVS_FG)),
                           data = dd, etype = "event"), paste("fg", cn, m))
    if (is.null(fgdat)) next
    fit <- safe(coxph(as.formula(paste0("Surv(fgstart,fgstop,fgstatus) ~ zm + ", COVS_FG)),
                      data = fgdat, weights = fgdat$fgwt), paste("fgfit", cn, m))
    if (is.null(fit)) next
    r <- getc(fit, "zm")
    # 对应加权 cause-specific Cox
    wrow <- cs[cs$cause == csmap[[cn]] & cs$metal == m, ]
    fg_out <- rbind(fg_out, data.frame(
      cause = cn, metal = m, n = nrow(dd), n_event = n_event,
      sHR = r["HR"], sHR_lo = r["lo"], sHR_hi = r["hi"], p = r["p"], weighted = FALSE,
      wCS_HR = if (nrow(wrow)) wrow$HR else NA, wCS_lo = if (nrow(wrow)) wrow$lo else NA,
      wCS_hi = if (nrow(wrow)) wrow$hi else NA, wCS_p = if (nrow(wrow)) wrow$p else NA))
  }
}
write.csv(fg_out, file.path(R, "rev6_finegray.csv"), row.names = FALSE)
cat("\n-- Fine-Gray 结果 --\n"); print(fg_out, row.names = FALSE)

# ##### 任务3：排除基线癌症的敏感性分析 #########################################
cat("\n\n=== 任务3: 排除基线癌症 敏感性分析 ===\n")
run_sub <- function(dd, label, outcomes = c("event_all", "event_cvd"), metals = MET3, covs = COVS_M3) {
  out <- data.frame()
  for (m in metals) for (ev in outcomes) {
    n_ev <- sum(dd[[ev]], na.rm = TRUE)
    fit <- safe(svycoxph(as.formula(paste0("Surv(time,", ev, ") ~ z_", m, " + ", covs)), design = mkdes(dd)), paste(label, m, ev))
    if (is.null(fit)) next
    r <- getc(fit, paste0("z_", m))
    out <- rbind(out, data.frame(subset = label, outcome = ev, metal = m, n = nrow(dd), n_deaths = n_ev,
                                 HR = r["HR"], lo = r["lo"], hi = r["hi"], p = r["p"]))
  }
  out
}
ex <- data.frame()
# 剔除基线自报癌症者(cancer_ever==1); 未报告(NA)并入非癌症组
subA <- d[is.na(d$cancer_ever) | d$cancer_ever == 0, ]
subB <- d[(is.na(d$cancer_ever) | d$cancer_ever == 0) & d$cvd == 0, ]
COVS_nocvd <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn"
ex <- rbind(ex, run_sub(subA, "excl_cancer"))
ex <- rbind(ex, run_sub(subB, "excl_cancer_and_cvd", covs = COVS_nocvd))  # 该样本 cvd 恒定, 需去掉 cvd 项
# 汞 cause-specific (全因/癌症死亡/脑血管死亡), 在排除基线癌症样本
ex <- rbind(ex, run_sub(subA, "excl_cancer_Hg_cause", outcomes = c("event_all", "event_cancer", "event_cerebro"), metals = "LBXTHG"))
write.csv(ex, file.path(R, "rev6_exclude_cancer.csv"), row.names = FALSE)
cat("\n-- 排除基线癌症结果 --\n"); print(ex, row.names = FALSE)

# ##### 任务4：逐步调整表 (stepwise) ############################################
cat("\n\n=== 任务4: 逐步调整表 ===\n")
blocks <- c(
  S1 = "",
  S2 = "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1)",
  S3 = "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR",
  S4 = "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI",
  S5 = "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke)",
  S6 = "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn",
  S7 = COVS_M3)
labels <- c(S1 = "1_crude", S2 = "2_+age/sex/race", S3 = "3_+educ/income", S4 = "4_+BMI",
            S5 = "5_+smoking", S6 = "6_+diabetes/HTN", S7 = "7_+baseline CVD(=M3)")
sw <- data.frame()
for (m in MET3) for (ev in c("event_all", "event_cvd")) {
  for (s in names(blocks)) {
    rhs <- paste0("z_", m, if (nzchar(blocks[[s]])) paste0(" + ", blocks[[s]]) else "")
    fit <- safe(svycoxph(as.formula(paste0("Surv(time,", ev, ") ~ ", rhs)), design = mkdes(d)), paste("sw", m, ev, s))
    if (is.null(fit)) next
    r <- getc(fit, paste0("z_", m))
    sw <- rbind(sw, data.frame(metal = m, outcome = ev, step = s, model = labels[[s]],
                               n = nrow(d), HR = r["HR"], lo = r["lo"], hi = r["hi"], p = r["p"]))
  }
}
write.csv(sw, file.path(R, "rev6_stepwise_mercury.csv"), row.names = FALSE)
cat("\n-- 逐步调整结果 --\n"); print(sw, row.names = FALSE)

# ##### 任务5：铅四分位的人年 ####################################################
cat("\n\n=== 任务5: 金属四分位人年 ===\n")
qpy <- data.frame()
for (m in MET3) {
  br <- quantile(d[[m]], probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
  qv <- cut(d[[m]], breaks = unique(br), include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))
  for (lv in c("Q1", "Q2", "Q3", "Q4")) {
    idx <- which(qv == lv)
    py <- sum(d$time[idx] / 12, na.rm = TRUE)  # 随访月 -> 年
    qpy <- rbind(qpy, data.frame(
      metal = m, Q = lv, n = length(idx), person_years = py,
      all_deaths = sum(d$event_all[idx], na.rm = TRUE),
      cvd_deaths = sum(d$event_cvd[idx], na.rm = TRUE),
      all_rate_per1000py = 1000 * sum(d$event_all[idx], na.rm = TRUE) / py,
      cvd_rate_per1000py = 1000 * sum(d$event_cvd[idx], na.rm = TRUE) / py,
      cut_lo = br[c(1, 2, 3, 4)][match(lv, c("Q1", "Q2", "Q3", "Q4"))],
      cut_hi = br[c(2, 3, 4, 5)][match(lv, c("Q1", "Q2", "Q3", "Q4"))]))
  }
}
write.csv(qpy, file.path(R, "rev6_quartile_personyears.csv"), row.names = FALSE)
cat("\n-- 四分位人年结果 --\n"); print(qpy, row.names = FALSE)

# ##### 任务2b：cmprsk::crr 交叉验证 (未加权, 竞争风险) ##########################
cat("\n\n=== 任务2b: cmprsk::crr 交叉验证 ===\n")
if (requireNamespace("cmprsk", quietly = TRUE)) {
  suppressMessages(library(cmprsk))
  crr_out <- data.frame()
  for (cn in names(causes)) {
    dd <- d
    dd$status <- ifelse(dd$mortstat == 0, 0L,
                        ifelse(!is.na(dd$ucod_leading) & dd$ucod_leading %in% causes[[cn]], 1L, 2L))
    cc <- dd[complete.cases(dd[, c("time", "status", paste0("z_", MET3), "RIDAGEYR", "RIAGENDR",
                                   "RIDRETH1", "educ", "INDFMPIR", "BMXBMI", "smoke", "dm", "htn", "cvd")]), ]
    for (m in MET3) {
      X <- model.matrix(as.formula(paste0("~ z_", m, " + ", COVS_M3)), data = cc)[, -1]
      fit <- safe(crr(cc$time, cc$status, cov1 = X, failcode = 1, cencode = 0), paste("crr", cn, m))
      if (is.null(fit)) next
      i <- match(paste0("z_", m), names(fit$coef)); b <- fit$coef[i]; se <- sqrt(as.numeric(fit$var[i, i]))
      crr_out <- rbind(crr_out, data.frame(cause = cn, metal = m, n = nrow(cc), n_event = sum(cc$status == 1),
                                           crr_sHR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
                                           p = 2 * pnorm(-abs(b / se)), weighted = FALSE))
    }
  }
  write.csv(crr_out, file.path(R, "rev6_finegray_crrcheck.csv"), row.names = FALSE)
  cat("-- crr 结果 --\n"); print(crr_out, row.names = FALSE)
} else cat("cmprsk 不可用, 跳过 crr 交叉验证\n")

cat("\n=== REV6-A DONE ===\n")
