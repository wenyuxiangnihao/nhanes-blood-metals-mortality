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
# REV7 Part C -- 缺失数据敏感性分析：多重插补 (MI) 与逆概率加权 (IPW)
# 回应审稿意见 M4: "Covariates analysed complete-case ... MNAR cannot be excluded"
#
# 数据:
#   data/pre_exclusion_df.rds  eligible n=36,743 (33,104 完整 + 3,639 协变量缺失; 含 `included`)
#   data/analysis_df_v6.rds    complete-case 基准 n=33,104
# 口径:
#   - ln/z 金属变换使用与主分析完全相同的常数 (由 33,104 计算)
#   - 协变量公式与主模型 M3 完全一致
#   - MI: M=20 个插补集, 链式方程 10 轮; 每个插补集在全部 36,743 上拟合
#   - IPW: 入选概率 logistic(在 36,743 上, 用同一插补集), w = wt / p̂ 归一化;
#          在 33,104 上加权拟合。含 1%/99% 截尾版本
#   - 合并: 手动 Rubin 规则 (点估计均值, 总方差 = W + (1+1/M)B)
# 输出: results/rev7_mi_ipw_summary.csv  (Table S20 数据源)
#       results/rev7_mi_ipw_detail.csv   (每插补集明细)
#       results/rev7_ipw_diagnostics.txt
# ============================================================
suppressMessages({library(survey); library(survival); library(MASS)})
options(survey.lonely.psu = "adjust")

ARGS <- commandArgs(trailingOnly = TRUE)
M    <- if (length(ARGS) >= 1) as.integer(ARGS[1]) else 20L
ITER <- if (length(ARGS) >= 2) as.integer(ARGS[2]) else 10L
SEED <- 20260925L
set.seed(SEED)

D <- file.path(ROOT, "data")
R <- file.path(ROOT, "results")
dir.create(R, showWarnings = FALSE)

MET3   <- c("LBXBPB", "LBXBCD", "LBXTHG")
COVS   <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) + INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd"
OUTS   <- c("event_all", "event_cvd")
IMPV   <- c("educ", "INDFMPIR", "BMXBMI", "smoke", "dm", "htn", "cvd")   # 需插补
PREDV  <- c("RIDAGEYR", "RIAGENDR", "RIDRETH1", "cycle", "time", "event_all", "event_cvd",
            "z_LBXBPB", "z_LBXBCD", "z_LBXTHG")                          # 完全观测预测子

cc <- readRDS(file.path(D, "analysis_df_v6.rds"))
el <- readRDS(file.path(D, "pre_exclusion_df.rds"))
stopifnot(nrow(el) == 36743L, nrow(cc) == 33104L, sum(el$included) == 33104L)

# --- ln/z 常数来自 complete-case 基准 (与主分析一致) ---
CONST <- sapply(MET3, function(m) { v <- log(cc[[m]]); c(mean = mean(v), sd = sd(v)) })
for (m in MET3) el[[paste0("z_", m)]] <- (log(el[[m]]) - CONST["mean", m]) / CONST["sd", m]
for (m in MET3) cc[[paste0("z_", m)]] <- (log(cc[[m]]) - CONST["mean", m]) / CONST["sd", m]   # complete-case 用同一常数

# --- 四分位切点也取自 complete-case 基准 ---
QCUT <- lapply(MET3, function(m) {
  b <- unique(quantile(cc[[m]], probs = seq(0, 1, 0.25), na.rm = TRUE))
  b[1] <- -Inf; b[length(b)] <- Inf; b
})
names(QCUT) <- MET3

cat("=== REV7-C MI/IPW START ===\n")
cat(sprintf("eligible=%d  included=%d  excluded=%d  M=%d  ITER=%d\n",
            nrow(el), sum(el$included), sum(el$included == 0), M, ITER))

# ============================================================
# 1. 链式方程多重插补 (MICE, 自行实现)
# ============================================================
draw_cat <- function(P, levs) {
  P <- P / rowSums(P)
  vapply(seq_len(nrow(P)), function(i) levs[sample.int(length(levs), 1L, prob = P[i, ])],
         integer(1))
}

impute_once <- function(dat) {
  d <- dat
  for (v in IMPV) {
    if (v %in% c("INDFMPIR", "BMXBMI")) d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)
    else d[[v]][is.na(d[[v]])] <- as.integer(names(which.max(table(d[[v]]))))
  }
  preds <- c(PREDV, IMPV)
  for (it in seq_len(ITER)) {
    for (v in IMPV) {
      rhs <- setdiff(preds, v)
      obs <- !is.na(dat[[v]]); mis <- is.na(dat[[v]])
      if (!any(mis)) next
      f <- as.formula(paste(v, "~", paste(rhs, collapse = " + ")))
      # 连续变量: 线性回归 + 残差正态抽样
      if (v %in% c("INDFMPIR", "BMXBMI")) {
        fit <- lm(f, data = d[obs, , drop = FALSE])
        s2  <- summary(fit)$sigma
        pr  <- predict(fit, newdata = d[mis, , drop = FALSE])
        val <- pr + rnorm(sum(mis), 0, s2)
        lo <- min(dat[[v]], na.rm = TRUE); hi <- max(dat[[v]], na.rm = TRUE)
        d[[v]][mis] <- pmin(pmax(val, lo), hi)
      } else if (v == "educ") {            # 1-5 有序 -> 比例优势模型
        fit <- try(polr(factor(educ) ~ ., data = d[obs, c("educ", rhs), drop = FALSE], Hess = FALSE), silent = TRUE)
        if (inherits(fit, "try-error")) {
          d[[v]][mis] <- draw_cat(tab <- prop.table(table(dat[[v]][obs])), as.integer(names(tab)))
        } else {
          P <- predict(fit, newdata = d[mis, , drop = FALSE], type = "probs")
          if (is.null(dim(P))) P <- matrix(P, nrow = sum(mis))
          d[[v]][mis] <- draw_cat(P, as.integer(colnames(P)))
        }
      } else if (v == "smoke") {           # 0/1/2 无序多分类
        fit <- try(nnet::multinom(f, data = d[obs, , drop = FALSE], trace = FALSE), silent = TRUE)
        if (inherits(fit, "try-error")) {
          d[[v]][mis] <- draw_cat(prop.table(table(dat[[v]][obs])), as.integer(names(table(dat[[v]][obs]))))
        } else {
          P <- predict(fit, newdata = d[mis, , drop = FALSE], type = "probs")
          if (is.null(dim(P))) P <- matrix(P, nrow = sum(mis))
          d[[v]][mis] <- draw_cat(P, as.integer(colnames(P)))
        }
      } else {                             # dm / htn / cvd : logistic
        fit <- try(glm(f, data = d[obs, , drop = FALSE], family = binomial()), silent = TRUE)
        if (inherits(fit, "try-error")) {
          d[[v]][mis] <- as.integer(runif(sum(mis)) < mean(dat[[v]][obs]))
        } else {
          p <- predict(fit, newdata = d[mis, , drop = FALSE], type = "response")
          d[[v]][mis] <- as.integer(runif(sum(mis)) < pmin(pmax(p, 1e-6), 1 - 1e-6))
        }
      }
    }
  }
  d
}

cat("\n[1] 生成 M =", M, "个插补数据集 ...\n")
imps <- vector("list", M)
for (i in seq_len(M)) {
  imps[[i]] <- impute_once(el)
  nchg <- sum(sapply(IMPV, function(v) sum(is.na(el[[v]]) != is.na(imps[[i]][[v]]))))
  cat(sprintf("    imputation %2d/%d  填补 %d 个格子\n", i, M, nchg))
}
saveRDS(imps, file.path(D, "rev7_imps.rds"))

# ============================================================
# 2. 估计器: 从 svycoxph 拟合中取单个暴露项的 beta / var / 及方向分解
# ============================================================
SE_CHK <- FALSE
get1 <- function(fit, term) {
  sm <- summary(fit); co <- sm$coefficients
  b  <- co[term, 1]
  # design-based (robust) 方差: vcov.svycoxph 返回带行名的设计基矩阵 (fit$var 无行名)
  V <- vcov(fit)
  se_rob <- sqrt(V[term, term])
  if (!SE_CHK && ncol(co) >= 4) {
    SE_CHK <<- TRUE
    cat(sprintf("    [check] robust se: vcov=%.6f  summary col4=%.6f  (col3 model-based=%.6f) -> %s\n",
                se_rob, co[term, 4], co[term, 3],
                ifelse(abs(se_rob - co[term, 4]) < 1e-8, "OK", "MISMATCH")))
  }
  list(b = b, v = se_rob^2, se_rob = se_rob, se_mod = co[term, 3])
}

rubin <- function(bs, vs) {                 # 标量 Rubin 规则
  M <- length(bs); bbar <- mean(bs)
  W <- mean(vs); B <- var(bs)
  Tv <- W + (1 + 1 / M) * B
  df <- if (B > 0) (M - 1) * (1 + W / ((1 + 1 / M) * B))^2 else Inf
  se <- sqrt(Tv)
  tc <- if (is.finite(df)) qt(0.975, df) else 1.959964
  c(est = bbar, se = se, lo = bbar - tc * se, hi = bbar + tc * se, df = df,
    fmi = (1 + 1 / M) * B / Tv)
}

# 方向分解: psi+ = sum(b>0), psi- = sum(b<0), 同一切点
dirsum <- function(fit, terms) {
  sm <- summary(fit); co <- sm$coefficients
  b <- co[terms, 1]
  V <- vcov(fit)                      # fit$var 无 dimnames, 必须用 vcov()
  V <- V[terms, terms, drop = FALSE]
  idx <- function(s) which(if (s > 0) b > 0 else b < 0)
  qf <- function(S) { if (!length(S)) return(c(est = 0, v = 0))
    w <- rep(1, length(S)); e <- sum(w * b[S]); v <- as.numeric(t(w) %*% V[S, S, drop = FALSE] %*% w)
    c(est = e, v = v) }
  list(pos = qf(idx(1)), neg = qf(idx(-1)))
}

# ============================================================
# 3. MI 分析 (在全部 36,743 上)
# ============================================================
cat("\n[2] MI 分析 (n = 36,743, 每个插补集) ...\n")
mi_b  <- list()   # key -> 数值向量
mi_v  <- list()
add <- function(k, b, v) { mi_b[[k]] <<- c(mi_b[[k]], b); mi_v[[k]] <<- c(mi_v[[k]], v) }

for (i in seq_len(M)) {
  di <- imps[[i]]
  des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = di)
  for (ev in OUTS) for (m in MET3) {
    f <- as.formula(paste0("Surv(time,", ev, ") ~ z_", m, " + ", COVS))
    r <- get1(svycoxph(f, design = des), paste0("z_", m))
    add(paste("single", ev, m, sep = "|"), r$b, r$v)
  }
  for (ev in OUTS) {                       # 方向分解 (四分位得分)
    dd <- di
    for (m in MET3) dd[[paste0("qc_", m)]] <- as.integer(cut(dd[[m]], breaks = QCUT[[m]],
                                                             labels = FALSE, include.lowest = TRUE)) - 1
    des2 <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = dd)
    f <- as.formula(paste0("Surv(time,", ev, ") ~ qc_LBXBPB + qc_LBXBCD + qc_LBXTHG + ", COVS))
    ds <- dirsum(svycoxph(f, design = des2), c("qc_LBXBPB", "qc_LBXBCD", "qc_LBXTHG"))
    add(paste("dirpos", ev, sep = "|"), ds$pos["est"], ds$pos["v"])
    add(paste("dirneg", ev, sep = "|"), ds$neg["est"], ds$neg["v"])
  }
  cat(sprintf("    MI %2d/%d done\n", i, M))
}
saveRDS(list(b = mi_b, v = mi_v), file.path(R, "rev7_mi_raw.rds"))

# ============================================================
# 4. IPW 分析
# ============================================================
cat("\n[3] IPW 分析 (入选概率模型 + 加权 Cox) ...\n")
ipw_b <- list(); ipw_v <- list(); ipw_t <- list()
add2 <- function(k, b, v) { ipw_b[[k]] <<- c(ipw_b[[k]], b); ipw_v[[k]] <<- c(ipw_v[[k]], v) }
diag_txt <- character()

for (i in seq_len(M)) {
  di <- imps[[i]]
  fs <- as.formula(paste("included ~", COVS, "+ z_LBXBPB + z_LBXBCD + z_LBXTHG"))
  pfit <- glm(fs, data = di, family = binomial())
  ph <- pmin(pmax(fitted(pfit), 0.02), 0.98)
  di$ipw_raw <- di$wt / ph
  inc <- di[di$included == 1, , drop = FALSE]
  inc$ipw <- inc$ipw_raw * (nrow(inc) / sum(inc$ipw_raw))     # 归一化到 n
  inc$ipw_t <- pmin(inc$ipw, quantile(inc$ipw, 0.99))
  inc$ipw_t <- inc$ipw_t * (nrow(inc) / sum(inc$ipw_t))
  if (i == 1) diag_txt <- c(diag_txt,
    sprintf("p(include): min=%.4f p1=%.4f median=%.4f max=%.4f",
            min(fitted(pfit)), quantile(fitted(pfit), .01), median(fitted(pfit)), max(fitted(pfit))),
    sprintf("ipw: min=%.2f median=%.2f max=%.2f  sum=%.0f  (min wtsq=%s)",
            min(inc$ipw), median(inc$ipw), max(inc$ipw), sum(inc$ipw),
            paste(sprintf("%.1f", min(inc$ipw)^2 * nrow(inc)), collapse = "")))
  for (tag in c("ipw", "ipw_t")) {
    des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = as.formula(paste0("~", tag)),
                     nest = TRUE, data = inc)
    for (ev in OUTS) for (m in MET3) {
      f <- as.formula(paste0("Surv(time,", ev, ") ~ z_", m, " + ", COVS))
      r <- get1(svycoxph(f, design = des), paste0("z_", m))
      add2(paste(tag, ev, m, sep = "|"), r$b, r$v)
    }
    for (ev in OUTS) {
      dd <- inc
      for (m in MET3) dd[[paste0("qc_", m)]] <- as.integer(cut(dd[[m]], breaks = QCUT[[m]],
                                                               labels = FALSE, include.lowest = TRUE)) - 1
      des2 <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA,
                        weights = as.formula(paste0("~", tag)), nest = TRUE, data = dd)
      f <- as.formula(paste0("Surv(time,", ev, ") ~ qc_LBXBPB + qc_LBXBCD + qc_LBXTHG + ", COVS))
      ds <- dirsum(svycoxph(f, design = des2), c("qc_LBXBPB", "qc_LBXBCD", "qc_LBXTHG"))
      add2(paste(tag, "dirpos", ev, sep = "|"), ds$pos["est"], ds$pos["v"])
      add2(paste(tag, "dirneg", ev, sep = "|"), ds$neg["est"], ds$neg["v"])
    }
  }
  cat(sprintf("    IPW %2d/%d done\n", i, M))
}
saveRDS(list(b = ipw_b, v = ipw_v), file.path(R, "rev7_ipw_raw.rds"))
writeLines(diag_txt, file.path(R, "rev7_ipw_diagnostics.txt"))

# ============================================================
# 5. complete-case 基准 (与主分析同口径, 用于并列比较)
# ============================================================
cat("\n[4] complete-case 基准 ...\n")
cc_b <- list(); cc_v <- list()
descc <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = cc)
add3 <- function(k, b, v) { cc_b[[k]] <<- b; cc_v[[k]] <<- v }
for (ev in OUTS) for (m in MET3) {
  f <- as.formula(paste0("Surv(time,", ev, ") ~ z_", m, " + ", COVS))
  r <- get1(svycoxph(f, design = descc), paste0("z_", m))
  add3(paste("single", ev, m, sep = "|"), r$b, r$v)
}
for (ev in OUTS) {
  dd <- cc
  for (m in MET3) dd[[paste0("qc_", m)]] <- as.integer(cut(dd[[m]], breaks = QCUT[[m]],
                                                           labels = FALSE, include.lowest = TRUE)) - 1
  des2 <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = dd)
  f <- as.formula(paste0("Surv(time,", ev, ") ~ qc_LBXBPB + qc_LBXBCD + qc_LBXTHG + ", COVS))
  ds <- dirsum(svycoxph(f, design = des2), c("qc_LBXBPB", "qc_LBXBCD", "qc_LBXTHG"))
  add3(paste("dirpos", ev, sep = "|"), ds$pos["est"], ds$pos["v"])
  add3(paste("dirneg", ev, sep = "|"), ds$neg["est"], ds$neg["v"])
}

# ============================================================
# 6. 汇总
# ============================================================
mk <- function(k) {
  if (is.null(mi_b[[k]])) return(NULL)
  out <- data.frame(key = k)
  kbase <- sub("^single\\|", "", k)          # IPW 侧 key 无 "single|" 前缀
  for (nm in c("MI", "IPW", "IPW_trunc")) {
    tg <- if (nm == "IPW_trunc") "ipw_t" else "ipw"
    src_b <- if (nm == "MI") mi_b[[k]] else ipw_b[[paste(tg, kbase, sep = "|")]]
    src_v <- if (nm == "MI") mi_v[[k]] else ipw_v[[paste(tg, kbase, sep = "|")]]
    if (is.null(src_b)) { out[[nm]] <- NA; next }
    r <- rubin(src_b, src_v)
    out[[nm]] <- sprintf("%.3f (%.3f-%.3f)", exp(r["est"]), exp(r["lo"]), exp(r["hi"]))
    if (nm == "MI") { out$MI_fmi <- round(r["fmi"], 3) }
  }
  if (!is.null(cc_b[[k]])) {
    b <- cc_b[[k]]; v <- cc_v[[k]]
    out$CC <- sprintf("%.3f (%.3f-%.3f)", exp(b - 1.96 * sqrt(v)), exp(b), exp(b + 1.96 * sqrt(v)))
  }
  out
}
keys <- c(paste("single", rep(OUTS, each = 3), MET3, sep = "|"),
          paste("dirpos", OUTS, sep = "|"), paste("dirneg", OUTS, sep = "|"))
res <- do.call(rbind, lapply(keys, mk))
write.csv(res, file.path(R, "rev7_mi_ipw_summary.csv"), row.names = FALSE)

det <- data.frame()
for (k in names(mi_b)) det <- rbind(det, data.frame(key = k, method = "MI", i = seq_along(mi_b[[k]]), est = mi_b[[k]]))
for (k in names(ipw_b)) det <- rbind(det, data.frame(key = k, method = "IPW", i = seq_along(ipw_b[[k]]), est = ipw_b[[k]]))
write.csv(det, file.path(R, "rev7_mi_ipw_detail.csv"), row.names = FALSE)

cat("\n=== 汇总 (Table S20) ===\n"); print(res, row.names = FALSE)
cat("\n=== IPW 诊断 ===\n"); cat(paste(diag_txt, collapse = "\n"), "\n")
cat("\n=== REV7-C DONE ===\n")
