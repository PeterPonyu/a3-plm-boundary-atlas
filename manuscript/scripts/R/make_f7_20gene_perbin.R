#!/usr/bin/env Rscript
# Inventory-only writer: full-20 per-bin reliability coordinates.
#
# Direct-run, capsule-local, R-first. Does NOT re-fit Platt.
# Does NOT emit a substitute mean_pred/obs_rate table.
#
# Pass condition (would write F7_20gene_reliability.csv): a file on disk that
# already stores 20-gene reliability-bin coordinates (bin_mid, mean_pred,
# obs_rate, n) for raw and calibrated scores. That object is not in this
# capsule. This script records the gap and exits non-zero.
suppressMessages({
  if (requireNamespace("jsonlite", quietly = TRUE)) library(jsonlite)
})

HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0L || identical(HERE, "")) HERE <- "scripts/R"
BASE <- normalizePath(file.path(HERE, "..", ".."))
CAP  <- normalizePath(file.path(BASE, ".."))
PD   <- file.path(BASE, "figures", "plotdata")
SNAP <- file.path(CAP, "results-snapshot", "N6_clinvar", "clinical_calibration_full20.json")
STATUS <- file.path(PD, "F7_20gene_perbin_STATUS.txt")

has <- function(p) file.exists(p)
rd  <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)
relpath <- function(p) {
  cap_slash <- paste0(CAP, .Platform$file.sep)
  if (startsWith(p, cap_slash)) sub(cap_slash, "", p, fixed = TRUE) else p
}

bin_like <- function(names) {
  grepl("bin_mid|mean_pred|obs_rate|reliability_point|per_bin", names, ignore.case = TRUE)
}

walk_names <- function(x, prefix = "") {
  out <- character()
  if (is.list(x)) {
    nm <- names(x)
    if (is.null(nm)) nm <- rep("", length(x))
    for (i in seq_along(x)) {
      key <- if (nzchar(nm[i])) paste0(prefix, nm[i]) else paste0(prefix, "[", i, "]")
      out <- c(out, key, walk_names(x[[i]], paste0(key, ".")))
    }
  }
  unique(out)
}

lines <- c(
  "A3 full-20 per-bin ECE inventory",
  "writer: manuscript/scripts/R/make_f7_20gene_perbin.R",
  "policy: no Platt re-fit; no synthesised mean_pred/obs_rate substitutes",
  paste("run_utc:", format(as.POSIXlt(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")),
  ""
)

# --- lock-token CSVs already in plotdata ---
summ_path <- file.path(PD, "F7_20gene_ece_summary.csv")
cov_path  <- file.path(PD, "F7_20gene_coverage.csv")
rel_path  <- file.path(PD, "F7_reliability.csv")
stopifnot(has(summ_path), has(cov_path), has(rel_path))
summ <- rd("F7_20gene_ece_summary.csv")
cov  <- rd("F7_20gene_coverage.csv")
rel_df <- rd("F7_reliability.csv")

lines <- c(lines,
  "PRESENT (capsule plotdata):",
  paste0("  ", relpath(summ_path)),
  paste0("    rows=", nrow(summ), " cols=", paste(names(summ), collapse = ",")),
  paste0("    full20 raw ECE=", summ$ece[summ$panel == "full20" & summ$score == "raw"],
         " calibrated=", summ$ece[summ$panel == "full20" & summ$score == "calibrated"],
         " n=", summ$n[summ$panel == "full20" & summ$score == "raw"],
         " n_heldout=", summ$n_heldout[summ$panel == "full20" & summ$score == "raw"]),
  paste0("  ", relpath(cov_path), "  (coverage curve; not reliability bins; rows=", nrow(cov), ")"),
  paste0("  ", relpath(rel_path), "  (SIX-GENE bins only; rows=", nrow(rel_df),
         "; which=", paste(unique(rel_df$which), collapse = ","), ")"),
  ""
)

# six-gene bin-weighted ECE must still match the lock tokens
ece_from_bins <- function(d) sum(d$n * abs(d$mean_pred - d$obs_rate)) / sum(d$n)
ece_raw_bins <- ece_from_bins(rel_df[rel_df$which == "raw", ])
ece_cal_bins <- ece_from_bins(rel_df[rel_df$which == "calibrated", ])
stopifnot(abs(ece_raw_bins - 0.1574) < 0.02)
stopifnot(abs(ece_cal_bins - 0.0466) < 0.02)
lines <- c(lines,
  sprintf("six-gene bin-weighted ECE from F7_reliability.csv: raw=%.4f cal=%.4f (lock 0.1574 / 0.0466)",
          ece_raw_bins, ece_cal_bins),
  ""
)

# --- pinned full20 JSON (tokens + coverage; no bins) ---
json_ok <- has(SNAP)
lines <- c(lines, paste0("PINNED JSON: ", relpath(SNAP), " exists=", json_ok))
json_has_coords <- FALSE
json_keys_hit <- character()
if (json_ok) {
  rawtxt <- readLines(SNAP, warn = FALSE)
  hit <- grep("bin_mid|mean_pred|obs_rate|reliability", rawtxt, ignore.case = TRUE, value = TRUE)
  json_keys_hit <- hit
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    j <- jsonlite::fromJSON(SNAP, simplifyVector = FALSE)
    alln <- walk_names(j)
    json_has_coords <- any(bin_like(alln))
    top <- names(j)
    lines <- c(lines, paste0("  top_keys: ", paste(top, collapse = ", ")))
    if (!is.null(j$full20)) {
      f20n <- names(j$full20)
      lines <- c(lines, paste0("  full20_keys: ", paste(f20n, collapse = ", ")))
      lines <- c(lines, sprintf("  full20 ece_raw=%s ece_calibrated=%s n=%s n_test_heldout=%s",
                                j$full20$ece_raw, j$full20$ece_calibrated,
                                j$full20$n, j$full20$n_test_heldout))
    }
  } else {
    json_has_coords <- length(hit) > 0L
    lines <- c(lines, "  jsonlite missing; fell back to text scan")
  }
  lines <- c(lines,
    paste0("  text_scan_binlike_lines=", length(json_keys_hit)),
    paste0("  json_has_bin_coordinates=", json_has_coords)
  )
}

# --- objects that would be needed to RECOMPUTE bins (forbidden here) ---
# Named so the missing list is exact. None of these live in the capsule.
missing <- c(
  "20-gene reliability-bin CSV with columns which,bin_mid,mean_pred,obs_rate,n (never written)",
  "clinical_calibration_full20.json reliability_points / per-bin arrays (JSON has pooled ECE tokens + cov_curve only)",
  "the original full20 calibration writer (not in this capsule; a3_clinical_calibration.py writes six-gene F7_reliability.csv only)",
  "capsule-local 20-gene per-variant held-out probabilities after the locked Platt fit (re-fitting Platt is forbidden)"
)
lines <- c(lines, "", "MISSING (cannot emit F7_20gene_reliability.csv):", paste0("  - ", missing), "")

candidate <- file.path(PD, "F7_20gene_reliability.csv")
if (has(candidate)) {
  cand <- rd("F7_20gene_reliability.csv")
  ok_cols <- all(c("which", "bin_mid", "mean_pred", "obs_rate", "n") %in% names(cand))
  lines <- c(lines, paste0("UNEXPECTED: ", relpath(candidate), " exists; ok_cols=", ok_cols,
                           " rows=", nrow(cand)))
  if (ok_cols && json_has_coords) {
    writeLines(lines, STATUS)
    cat("unexpected recoverable table at", candidate, "\n")
    quit(save = "no", status = 0)
  }
}

lines <- c(lines,
  "VERDICT: STOP. Full-20 per-bin coordinates are not recoverable from capsule files.",
  "No substitute CSV written.",
  paste0("status_file: ", relpath(STATUS))
)
writeLines(lines, STATUS)
cat(paste(lines, collapse = "\n"), "\n", sep = "")
quit(save = "no", status = 1)
