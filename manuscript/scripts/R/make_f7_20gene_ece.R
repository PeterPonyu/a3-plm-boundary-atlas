#!/usr/bin/env Rscript
# A3 dedicated 20-gene ECE / abstention float (A3-8).
# Reads existing F7_reliability.csv (six-gene bins) plus locked full20 summaries.
# Does not re-fit Platt. Does not invent 20-gene per-bin coordinates.
suppressMessages({ library(ggplot2); library(patchwork); library(dplyr) })
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "scripts/R"
source(file.path(HERE, "theme_a3.R"))

BASE <- normalizePath(file.path(HERE, "..", ".."))
PD   <- file.path(BASE, "figures", "plotdata")
FIG  <- Sys.getenv("A3_FIGURE_OUTPUT", unset = "")
if (FIG == "") stop("Set A3_FIGURE_OUTPUT to an explicit staging directory.", call. = FALSE)
if (!dir.exists(FIG)) stop(sprintf("A3_FIGURE_OUTPUT does not exist: %s", FIG), call. = FALSE)
rd   <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)

rel <- rd("F7_reliability.csv")
ece_from_bins <- function(d) {
  sum(d$n * abs(d$mean_pred - d$obs_rate)) / sum(d$n)
}
ece_raw_bins <- ece_from_bins(rel[rel$which == "raw", ])
ece_cal_bins <- ece_from_bins(rel[rel$which == "calibrated", ])
# Six-gene lock tokens are 0.1574 / 0.0466. Bin-weighted ECE must stay near those.
stopifnot(abs(ece_raw_bins - 0.1574) < 0.02)
stopifnot(abs(ece_cal_bins - 0.0466) < 0.02)

summ <- rd("F7_20gene_ece_summary.csv")
summ$panel <- factor(summ$panel, levels = c("full20", "sixgene"),
                     labels = c("20-gene panel", "six-gene subset"))
summ$score <- factor(summ$score, levels = c("raw", "calibrated"),
                     labels = c("raw", "Platt-calibrated"))

pA <- ggplot(summ, aes(panel, ece, fill = score)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68) +
  geom_text(aes(label = sprintf("%.3f", ece)),
            position = position_dodge(width = 0.72),
            vjust = -0.35, size = 2.6, family = A3_FONT) +
  scale_fill_manual(values = c("raw" = oi("verm"), "Platt-calibrated" = oi("blue"))) +
  scale_y_continuous(limits = c(0, 0.30), expand = expansion(c(0, 0.06))) +
  labs(x = NULL, y = "pooled ECE (held-out)") +
  theme_a3() + a3_legend_outside("bottom")

cov <- rd("F7_20gene_coverage.csv")
pB <- ggplot(cov, aes(coverage, accuracy)) +
  geom_hline(yintercept = 0.7631, linewidth = 0.35, linetype = 2, colour = "grey50") +
  geom_line(linewidth = 0.7, colour = oi("blue")) +
  geom_point(size = 2.0, colour = oi("blue")) +
  scale_x_continuous(limits = c(0.05, 1.0), breaks = seq(0.2, 1.0, 0.2)) +
  scale_y_continuous(limits = c(0.70, 0.90)) +
  labs(x = "coverage (most-confident fraction)", y = "overall accuracy") +
  theme_a3()

prec <- cov %>% filter(!is.na(pathogenic_precision))
pC <- ggplot(prec, aes(coverage, pathogenic_precision)) +
  geom_vline(xintercept = 0.6, linewidth = 0.35, linetype = 2, colour = "grey50") +
  geom_line(linewidth = 0.7, colour = oi("orange")) +
  geom_point(size = 2.0, colour = oi("orange")) +
  scale_x_continuous(limits = c(0.05, 1.0), breaks = seq(0.2, 1.0, 0.2)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "coverage", y = "pathogenic-call precision") +
  theme_a3()

rel$which <- factor(rel$which, levels = c("raw", "calibrated"),
                    labels = c("raw score", "Platt-calibrated"))
pD <- ggplot(rel, aes(mean_pred, obs_rate, colour = which)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4, colour = "grey55", linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  geom_point(aes(size = n)) +
  scale_colour_manual(values = c("raw score" = oi("verm"), "Platt-calibrated" = oi("blue"))) +
  scale_size_continuous(range = c(1.0, 3.6), guide = "none") +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "predicted P(pathogenic)", y = "observed pathogenic rate") +
  theme_a3() + a3_legend_outside("bottom")

fig <- (pA | pB) / (pC | pD) + plot_annotation(tag_levels = "A") &
  # Single tag scale across the capsule: plot_annotation tags inherit the
  # theme_a3 margin band; only the size is pinned here to the shared 9pt.
  theme(plot.tag = element_text(size = A3_PANEL_LABEL_PT, face = "bold", family = A3_FONT))
a3_save(fig, "F7_20gene_ece", 9.2, 7.0, FIG)
cat(sprintf("six-gene bin-weighted ECE raw=%.4f cal=%.4f (lock 0.1574 / 0.0466)\n",
            ece_raw_bins, ece_cal_bins))
