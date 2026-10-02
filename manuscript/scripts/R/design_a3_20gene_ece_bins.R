#!/usr/bin/env Rscript
# Unique-writer DESIGN for A3 20-gene per-bin ECE (leftover8-comment-response).
#
# Cluster: A3-LATER-ECE
# Window: after A3 printable clusters, before the remaining seven papers.
# This ralplan session: run_this_session = FALSE. Do not invoke this file now.
#
# Purpose (later, when a human tickets a run):
#   If — and only if — 20-gene reliability-bin rows already exist on disk
#   (columns include mean_pred and obs_rate), emit an auditable per-bin table
#   and stop. Never invent bins. Never re-fit Platt to reprint 0.2443 / 0.0207.
#
# This file is the design. The inventory-only sibling
#   make_f7_20gene_perbin.R
# already recorded VERDICT STOP. Do not treat that inventory as a bin writer.

RUN_THIS_SESSION <- FALSE

HERE <- tryCatch(
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
  error = function(e) "."
)
if (length(HERE) == 0L || identical(HERE, "")) HERE <- "scripts/R"
BASE <- normalizePath(file.path(HERE, "..", ".."))
CAP  <- normalizePath(file.path(BASE, ".."))
PD   <- file.path(BASE, "figures", "plotdata")

# ---- required inputs (capsule-relative; fail-closed if absent) ----
# Pass path (would allow a later write of F7_20gene_reliability.csv):
#   a table whose columns include mean_pred and obs_rate for 20-gene
#   held-out reliability bins, raw and calibrated, already computed from
#   locked per-variant probabilities. That table is not in this capsule.
REQUIRED_BIN_TABLE_CANDIDATES <- c(
  file.path(PD, "F7_20gene_reliability.csv"),
  file.path(PD, "F7_20gene_perbin.csv")
)
# Token-only objects (present; NOT sufficient to recover bins):
TOKEN_ONLY <- c(
  file.path(PD, "F7_20gene_ece_summary.csv"),          # pooled ECE tokens
  file.path(PD, "F7_20gene_coverage.csv"),             # coverage curve, not bins
  file.path(CAP, "results-snapshot", "N6_clinvar", "clinical_calibration_full20.json")
)
# Six-gene bins are a different panel; do not paste them as 20-gene bins:
SIXGENE_BINS_ONLY <- file.path(PD, "F7_reliability.csv")

LOCKED_ECE_RAW <- 0.2443
LOCKED_ECE_CAL <- 0.0207

fail_closed <- function(msg) {
  cat("FAIL-CLOSED: ", msg, "\n", sep = "")
  quit(save = "no", status = 2)
}

if (!isTRUE(RUN_THIS_SESSION)) {
  cat(paste0(
    "DESIGN ONLY \u2014 leftover8-comment-response A3-LATER-ECE.\n",
    "run_this_session=FALSE. Unique writer not executed.\n",
    "Required later: a 20-gene table with columns mean_pred and obs_rate.\n",
    "Absent those rows: STOP. Do not synthesise bins.\n",
    "Forbidden: Platt re-fit to manufacture ", LOCKED_ECE_RAW, " / ",
    LOCKED_ECE_CAL, ".\n",
    "Forbidden: copy six-gene F7_reliability.csv rows onto the 20-gene panel.\n",
    "Token-only inputs (not sufficient): F7_20gene_ece_summary.csv, ",
    "F7_20gene_coverage.csv, clinical_calibration_full20.json (pooled ECE + cov_curve).\n"
  ))
  quit(save = "no", status = 0)
}

# ---- later-run gate (dead this session because RUN_THIS_SESSION is FALSE) ----
has_mean_pred_obs <- function(path) {
  if (!file.exists(path)) return(FALSE)
  d <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(d) || nrow(d) < 1L) return(FALSE)
  all(c("mean_pred", "obs_rate") %in% names(d))
}

ok <- vapply(REQUIRED_BIN_TABLE_CANDIDATES, has_mean_pred_obs, logical(1))
if (!any(ok)) {
  fail_closed(paste0(
    "no 20-gene mean_pred/obs_rate rows on disk. ",
    "Will not synthesise bins and will not re-fit Platt to reprint ",
    LOCKED_ECE_RAW, "/", LOCKED_ECE_CAL, "."
  ))
}

# If a real bin table appears later: copy-through only. No Platt. No matplotlib.
src <- REQUIRED_BIN_TABLE_CANDIDATES[which(ok)[1]]
d <- read.csv(src, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("which", "bin_mid", "mean_pred", "obs_rate", "n")
if (!all(need %in% names(d))) {
  fail_closed(paste("bin table missing columns:", paste(setdiff(need, names(d)), collapse = ",")))
}
out <- file.path(PD, "F7_20gene_reliability.csv")
if (normalizePath(src, mustWork = FALSE) != normalizePath(out, mustWork = FALSE)) {
  write.csv(d[need], out, row.names = FALSE)
}
cat("copy-through 20-gene bins from existing mean_pred/obs_rate rows; Platt not re-fit\n")
