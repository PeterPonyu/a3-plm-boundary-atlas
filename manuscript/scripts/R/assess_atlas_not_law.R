#!/usr/bin/env Rscript
# Honest science-bar inventory from capsule CSVs only.
# Does not fit a new slope, does not invent a law, does not write plotdata locks.
# Prints whether a monotone quantitative physical-world law is on disk.
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0L || identical(HERE, "")) HERE <- "scripts/R"
BASE <- normalizePath(file.path(HERE, "..", ".."))
PD   <- file.path(BASE, "figures", "plotdata")
OUT  <- file.path(PD, "atlas_not_law_STATUS.txt")
rel_out <- "manuscript/figures/plotdata/atlas_not_law_STATUS.txt"
rd   <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)

bind <- rd("F4_binding.csv")
disc <- rd("F10_discrimination.csv")
cal  <- rd("F10_calibration.csv")
atl  <- rd("F11_atlas_spearman.csv")
ece  <- rd("F7_20gene_ece_summary.csv")
cov  <- rd("F7_20gene_coverage.csv")

ci_overlap0 <- function(lo, hi) lo <= 0 && hi >= 0
ci_overlap  <- function(a_lo, a_hi, b_lo, b_hi) !(a_hi < b_lo || b_hi < a_lo)

row_of <- function(d, key) d[d$key == key, , drop = FALSE]
ov  <- row_of(bind, "overall")
ab  <- row_of(bind, "antibody_ABAG")
ni  <- row_of(bind, "noninterface")
gi  <- row_of(bind, "general_nonAB")
ifc <- row_of(bind, "interface")

# Held-out MegaScale Spearman / RMSE ranges from committed tables (no new fit)
sp_ho <- range(disc$spearman[disc$split == "homology-held-out"])
rmse_all <- range(cal$calib_rmse)

f20 <- ece[ece$panel == "full20", ]
cov_full <- cov[abs(cov$coverage - 1.0) < 1e-9, ]
cov_half <- cov[abs(cov$coverage - 0.5) < 1e-9, ]

# Coverage-accuracy Spearman on the 10 committed operating points (descriptive only)
cov_rho <- suppressWarnings(cor(cov$coverage, cov$accuracy, method = "spearman"))

lines <- c(
  "A3 science-bar inventory (capsule CSVs only)",
  "option: (b) atlas, not a fitted physical-world law",
  paste("run_utc:", format(as.POSIXlt(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")),
  "",
  "LOCKED PHYSICAL READOUTS (rank correlation / RMSE, not a kcal slope):",
  sprintf("  SKEMPI overall Spearman=%.6f  95%% CI [%.6f, %.6f]  n=%d  CI_includes_0=%s",
          ov$spearman, ov$ci_lo, ov$ci_hi, ov$n, ci_overlap0(ov$ci_lo, ov$ci_hi)),
  sprintf("  general_nonAB Spearman=%.6f  95%% CI [%.6f, %.6f]  n=%d  CI_includes_0=%s",
          gi$spearman, gi$ci_lo, gi$ci_hi, gi$n, ci_overlap0(gi$ci_lo, gi$ci_hi)),
  sprintf("  antibody_ABAG Spearman=%.6f  95%% CI [%.6f, %.6f]  n=%d  CI_includes_0=%s",
          ab$spearman, ab$ci_lo, ab$ci_hi, ab$n, ci_overlap0(ab$ci_lo, ab$ci_hi)),
  sprintf("  interface Spearman=%.6f  95%% CI [%.6f, %.6f]  n=%d",
          ifc$spearman, ifc$ci_lo, ifc$ci_hi, ifc$n),
  sprintf("  noninterface Spearman=%.6f  95%% CI [%.6f, %.6f]  n=%d  CI_includes_0=%s",
          ni$spearman, ni$ci_lo, ni$ci_hi, ni$n, ci_overlap0(ni$ci_lo, ni$ci_hi)),
  sprintf("  interface vs noninterface CI_overlap=%s",
          ci_overlap(ifc$ci_lo, ifc$ci_hi, ni$ci_lo, ni$ci_hi)),
  sprintf("  antibody_ABAG vs general_nonAB CI_overlap=%s",
          ci_overlap(ab$ci_lo, ab$ci_hi, gi$ci_lo, gi$ci_hi)),
  "",
  "F11_atlas_spearman.csv (three-decimal synthesis face):",
  paste0("  ", apply(atl, 1, function(r) paste(r, collapse = " | "))),
  "",
  sprintf("  MegaScale homology-held-out Spearman locked cells from F10_discrimination.csv: v2 WT %.3f, v2 masked %.3f",
          disc$spearman[disc$variant == "v2 WT" & disc$split == "homology-held-out"],
          disc$spearman[disc$variant == "v2 masked" & disc$split == "homology-held-out"]),
  "  stale range [0.336, 0.362] SUPERSEDED (0.336 is not a cell in F10_discrimination.csv)",
  sprintf("  MegaScale calib_rmse range [%.3f, %.3f] from F10_calibration.csv",
          rmse_all[1], rmse_all[2]),
  "",
  "CALIBRATION / ABSTENTION TOKENS (not a physical law):",
  sprintf("  full20 ECE raw=%s calibrated=%s n=%s n_heldout=%s",
          f20$ece[f20$score == "raw"], f20$ece[f20$score == "calibrated"],
          f20$n[f20$score == "raw"], f20$n_heldout[f20$score == "raw"]),
  sprintf("  coverage-accuracy Spearman (10 operating points)=%.4f  [descriptive; n_points=%d]",
          cov_rho, nrow(cov)),
  sprintf("  full coverage accuracy=%s  half coverage accuracy=%s",
          cov_full$accuracy, cov_half$accuracy),
  "",
  "WHY THIS IS NOT OPTION (a):",
  "  - No committed kcal/mol slope, intercept, or residual law vs experimental \u0394\u0394G.",
  "  - Antibody\u2013antigen Spearman CI lower bound is 0.00179 (near 0); not a judged A\u2260B collapse.",
  "  - Non-interface CI includes 0; manuscript already scopes it as corroborative / range-restricted.",
  "  - Coverage vs accuracy is an in-panel operating curve on ClinVar, not a physical measurement.",
  "  - ECE drop is a ClinVar calibration token (in-panel Platt), not transport to a new lab.",
  "",
  "MEASUREMENT THAT WOULD BE REQUIRED FOR A LAW PAPER:",
  "  Paired experimental free-energy or binding measurements (kcal/mol) on a physically",
  "  held-out condition (new partner, temperature, or complex), with a fitted quantitative",
  "  map (slope, intercept, residual distribution) whose uncertainty is estimated on that",
  "  hold-out; and/or an experimental pathogenic/benign outcome measured outside UniRef and",
  "  ClinVar. Another AUROC or ECE table would not clear the bar.",
  "",
  "VERDICT: option (b). Reliability atlas. Do not manufacture a law figure."
)
writeLines(lines, OUT)
cat(paste(lines, collapse = "\n"), "\n", sep = "")
quit(save = "no", status = 0)
