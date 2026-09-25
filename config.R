# ============================================================
# config.R -- single place that defines every project path.
#
# Every R script sources this file at start-up.  Values can be
# overridden with environment variables (useful for HPC / CI):
#
#   RAW_DIR      -> directory holding the downloaded NHANES *.XPT files
#   DATA_DIR     -> directory for cleaned analysis datasets (*.rds)
#   RESULTS_DIR  -> directory for all result tables (*.csv)
#   FIG_DIR      -> directory for all figures (png / tiff / eps)
#   HM_CONFIG    -> path to this file if it is not next to the script
#
# Example:
#   RAW_DIR=/data/nhanes RESULTS_DIR=/tmp/res Rscript rcs3.R
# ============================================================

## Repo root = directory that contains this config.R
.hm_cfg <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
.hm_root <- if (!is.na(.hm_cfg)) dirname(.hm_cfg) else getwd()

RAW_DIR     <- Sys.getenv("RAW_DIR",     unset = file.path(.hm_root, "data_raw"))
DATA_DIR    <- Sys.getenv("DATA_DIR",    unset = file.path(.hm_root, "data"))
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = file.path(.hm_root, "results"))
FIG_DIR     <- Sys.getenv("FIG_DIR",     unset = file.path(.hm_root, "figures"))

## make sure the writable output directories exist
invisible(lapply(c(DATA_DIR, RESULTS_DIR, FIG_DIR),
                 dir.create, showWarnings = FALSE, recursive = TRUE))

message(sprintf("[config] root=%s\n         RAW_DIR=%s\n         DATA_DIR=%s\n         RESULTS_DIR=%s\n         FIG_DIR=%s",
                .hm_root, RAW_DIR, DATA_DIR, RESULTS_DIR, FIG_DIR))
