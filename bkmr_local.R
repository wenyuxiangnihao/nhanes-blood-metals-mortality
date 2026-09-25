# =============================================================================
# Blood metal mixture x mortality  --  BKMR auto-runner   (v10)
#
#   v10:
#     * REUSES an existing bkmr_fit.rds if present (skips the 20-40 min MCMC)
#     * plotting no longer assumes a column called 'sd' (column names of
#       OverallRiskSummaries / PredictorResponseUnivar differ across bkmr
#       versions) -- it auto-detects them and prints them for the record
#
#   Rscript bkmr_local.R                 # default: N = 800, iter = 2000
#   Rscript bkmr_local.R quick           # N = 500, iter = 1000
#   Rscript bkmr_local.R n=600 iter=3000
#   (delete bkmr_fit.rds if you want a fresh fit)
# =============================================================================
a <- commandArgs(trailingOnly = FALSE)
fa <- grep("^--file=", a, value = TRUE)
here <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa))) else getwd()
setwd(here); cat("Working dir:", getwd(), "\n")

argv <- commandArgs(trailingOnly = TRUE)
getnum <- function(key, default) {
  hit <- grep(paste0("^", key, "="), argv, value = TRUE)
  if (length(hit)) as.numeric(sub(paste0("^", key, "="), "", hit[1])) else default
}
QUICK <- ("quick" %in% argv)
ITER  <- if (QUICK) 1000 else getnum("iter", 2000)
NMAX  <- if (QUICK) 500  else getnum("n", 800)
FITFILE <- "bkmr_fit.rds"

DATA <- "heavymetal_analysis_for_local.csv"
if (!file.exists(DATA)) {
  cand <- list.files(pattern = "\\.csv$")
  if (length(cand) == 1) DATA <- cand[1] else stop("Data file not found in: ", here)
}
suppressMessages(library(bkmr))
cat("bkmr version:", as.character(packageVersion("bkmr")), "\n")
set.seed(20260917)
metals <- c("LBXBPB", "LBXBCD", "LBXTHG")

## ================= fit (or reuse) =================
if (file.exists(FITFILE)) {
  cat("\n>>> Found", FITFILE, "- reusing it (skipping MCMC).\n")
  fit <- readRDS(FITFILE)
} else {
  d <- read.csv(DATA, stringsAsFactors = FALSE)
  cat("Full N =", nrow(d), " | all-cause deaths =", sum(d$event_all), "\n")
  if (nrow(d) > NMAX) {
    d <- d[sample(nrow(d), NMAX), ]
    cat("Subsampled (random) to N =", nrow(d), " | deaths =", sum(d$event_all), "\n")
  }
  mkZ <- function(dat) {
    n <- nrow(dat)
    Z <- sapply(metals, function(m) qnorm((rank(dat[[m]], ties.method = "first") - 0.5) / n))
    matrix(Z, nrow = n, dimnames = list(NULL, c("Lead", "Cadmium", "Mercury")))
  }
  mkX <- function(form, dat) {
    X <- model.matrix(form, data = dat)[, -1, drop = FALSE]
    s <- apply(X, 2, function(v) sd(v, na.rm = TRUE))
    X <- X[, is.finite(s) & s > 1e-8, drop = FALSE]
    q <- qr(X); if (q$rank < ncol(X)) X <- X[, q$pivot[seq_len(q$rank)], drop = FALSE]
    scale(X)
  }
  FORM_FULL <- ~ RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(educ) +
    INDFMPIR + BMXBMI + factor(smoke) + dm + htn + cvd + log(pmax(time, 1))
  FORM_SLIM <- ~ RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + BMXBMI + factor(smoke)
  Z <- mkZ(d); y <- d$event_all
  cat("Estimated run time ~", round(2.2e-9 * nrow(d)^3 * ITER / 60, 1), "min\n")
  fit <- NULL; err <- NULL
  for (at in list(list(X = mkX(FORM_FULL, d), fam = "binomial", tag = "binomial | fullX"),
                  list(X = mkX(FORM_SLIM, d), fam = "binomial", tag = "binomial | slimX"),
                  list(X = mkX(FORM_SLIM, d), fam = "gaussian", tag = "gaussian | slimX"))) {
    cat(">> trying:", at$tag, "...\n")
    fit <- tryCatch(kmbayes(y = y, Z = Z, X = at$X, iter = ITER, family = at$fam, verbose = TRUE),
                    error = function(e) { err <<- conditionMessage(e); NULL })
    if (!is.null(fit)) { cat(">> SUCCESS:", at$tag, "\n"); break }
    cat("!! failed:", at$tag, "->", err, "\n\n")
  }
  if (is.null(fit)) stop("All attempts failed. Last error: ", err)
  saveRDS(fit, FITFILE); gc()
}

## ================= helpers =================
pick <- function(df, cand) { for (cc in cand) if (cc %in% names(df)) return(df[[cc]]); NULL }

## ================= 1. overall mixture effect =================
qs <- seq(0.1, 0.9, by = 0.1)
ov <- tryCatch(OverallRiskSummaries(fit, qs = qs, q.fixed = 0.5, method = "approx"),
               error = function(e) OverallRiskSummaries(fit, qs = qs, q.fixed = 0.5))
cat("\n--- names(OverallRiskSummaries) ---\n"); print(names(ov)); print(head(ov, 10))
write.csv(ov, "bkmr_overall_summaries.csv", row.names = FALSE)

qy  <- pick(ov, c("q", "quantile", "qs", "qval"))
est <- pick(ov, c("est", "estimate", "mean", "beta", "risk"))
sdv <- pick(ov, c("sd", "se", "std", "sd.est", "est_sd", "sd_est"))
if (is.null(est)) est <- ov[[ncol(ov)]]
if (is.null(qy))  qy  <- seq_along(est)

png("bkmr_overall.png", width = 1400, height = 1100, res = 180)
if (!is.null(sdv) && length(sdv) == length(est)) {
  lo <- est - 1.96 * sdv; hi <- est + 1.96 * sdv
  plot(est, qy, type = "b", pch = 19, xlim = range(c(lo, hi, 0), finite = TRUE),
       xlab = "effect on mortality scale", ylab = "Quantile of all metals (vs median)",
       main = "Overall mixture effect (BKMR)")
  segments(lo, qy, hi, qy)
} else {
  plot(est, qy, type = "b", pch = 19,
       xlab = "effect on mortality scale", ylab = "Quantile of all metals (vs median)",
       main = "Overall mixture effect (BKMR)")
}
abline(v = 0, lty = 2, col = "grey"); dev.off()
cat("wrote bkmr_overall.png\n")

## ================= 2. single-metal dose-response =================
pr <- tryCatch(PredictorResponseUnivar(fit, q.fixed = 0.5, method = "approx"),
               error = function(e) PredictorResponseUnivar(fit, q.fixed = 0.5))
cat("\n--- names(PredictorResponseUnivar) ---\n"); print(names(pr))
write.csv(pr, "bkmr_single_data.csv", row.names = FALSE)

vvar <- pick(pr, c("variable", "var", "exposure", "varname"))
vz   <- pick(pr, c("z", "zval", "x", "value"))
vest <- pick(pr, c("est", "estimate", "mean"))
vsd  <- pick(pr, c("sd", "se", "std", "sd.est", "est_sd", "sd_est"))

if (!is.null(vvar) && !is.null(vz) && !is.null(vest)) {
  labs <- unique(as.character(vvar))
  png("bkmr_single.png", width = 1900, height = 650, res = 170)
  par(mfrow = c(1, length(labs)))
  for (m in labs) {
    idx <- as.character(vvar) == m
    zz <- vz[idx]; ee <- vest[idx]
    o <- order(zz); zz <- zz[o]; ee <- ee[o]
    if (!is.null(vsd)) { ss <- vsd[idx][o] } else ss <- NULL
    ylim <- if (!is.null(ss)) range(c(ee - 1.96 * ss, ee + 1.96 * ss), finite = TRUE) else range(ee, finite = TRUE)
    plot(zz, ee, type = "l", lwd = 2, xlab = paste0(m, " (score)"), ylab = "h(m)",
         main = m, ylim = ylim)
    if (!is.null(ss)) polygon(c(zz, rev(zz)), c(ee - 1.96 * ss, rev(ee + 1.96 * ss)),
                              col = rgb(.8, .2, .2, .2), border = NA)
    abline(h = 0, lty = 2, col = "grey")
  }
  dev.off(); cat("wrote bkmr_single.png\n")
} else {
  cat("!! could not identify columns for the single-metal plot; see names() above.\n")
}

cat("\nDONE. Files in", here, ":\n")
cat("  bkmr_overall.png, bkmr_single.png, bkmr_overall_summaries.csv, bkmr_single_data.csv\n")
try(shell.exec(file.path(here, "bkmr_overall.png")), silent = TRUE)
