# =============================================================================
# bkmr diagnostic: is the binomial (probit) sampler usable in this setup?
#   Double-click  test_bkmr.bat  to run.  Takes about 1 minute.
#   It tests simulated data and a small slice of your real data, with and
#   without an intercept column, binomial vs gaussian.
# =============================================================================
a <- commandArgs(trailingOnly = FALSE)
fa <- grep("^--file=", a, value = TRUE)
here <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa))) else getwd()
setwd(here)
suppressMessages(library(bkmr))
cat("R   :", R.version.string, "\n")
cat("bkmr:", as.character(packageVersion("bkmr")), "\n")
cat("dir :", here, "\n\n")

res <- character(0)
run <- function(lab, expr) {
  cat("--", lab, " ... ")
  out <- tryCatch({ force(expr); "OK" }, error = function(e) paste("FAIL:", conditionMessage(e)))
  cat(out, "\n")
  res <<- c(res, paste0(lab, "  =>  ", out))
  invisible(out)
}

set.seed(1)
n <- 300
Z <- matrix(rnorm(n * 3), n, 3); colnames(Z) <- c("a", "b", "c")
X <- cbind(age = rnorm(n), sex = rbinom(n, 1, 0.5))
y <- rbinom(n, 1, plogis(0.6 * Z[, 1] - 0.4 * Z[, 2] + 0.1 * X[, 1]))

cat("=== A. simulated data (n=300, iter=300) ===\n")
run("A1 binomial | X no intercept",
    kmbayes(y = y, Z = Z, X = X, iter = 300, family = "binomial", verbose = FALSE))
run("A2 binomial | X WITH intercept",
    kmbayes(y = y, Z = Z, X = cbind(intercept = 1, X), iter = 300, family = "binomial", verbose = FALSE))
run("A3 gaussian | X no intercept",
    kmbayes(y = y, Z = Z, X = X, iter = 300, family = "gaussian", verbose = FALSE))

cat("\n=== B. your real data (n=500, iter=300) ===\n")
if (file.exists("heavymetal_analysis_for_local.csv")) {
  d <- read.csv("heavymetal_analysis_for_local.csv", stringsAsFactors = FALSE)
  set.seed(2); d <- d[sample(nrow(d), 500), ]
  mets <- c("LBXBPB", "LBXBCD", "LBXTHG")
  Z2 <- as.matrix(scale(log(as.matrix(d[, mets]))))
  colnames(Z2) <- c("Lead", "Cadmium", "Mercury")
  X2 <- model.matrix(~ RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + BMXBMI + factor(smoke),
                     data = d)[, -1, drop = FALSE]
  y2 <- d$event_all
  run("B1 real | binomial | X no intercept",
      kmbayes(y = y2, Z = Z2, X = X2, iter = 300, family = "binomial", verbose = FALSE))
  run("B2 real | binomial | X WITH intercept",
      kmbayes(y = y2, Z = Z2, X = cbind(intercept = 1, X2), iter = 300, family = "binomial", verbose = FALSE))
  run("B3 real | gaussian | X no intercept",
      kmbayes(y = y2, Z = Z2, X = X2, iter = 300, family = "gaussian", verbose = FALSE))
} else {
  cat("(heavymetal_analysis_for_local.csv not found in this folder - skipped)\n")
}

cat("\n===== SUMMARY =====\n")
cat(paste(res, collapse = "\n"), "\n")
