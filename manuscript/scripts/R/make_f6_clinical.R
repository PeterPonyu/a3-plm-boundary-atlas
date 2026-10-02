#!/usr/bin/env Rscript
# A3 F7 — clinical (ClinVar) calibration + abstention, ggplot2 + Times New Roman (Arial).
# Two panels: (a) reliability diagram (raw vs Platt-calibrated, pooled held-out ESM-2 P(pathogenic))
#             (b) selective abstention: accuracy + pathogenic-call precision vs coverage.
# Reads tidy CSVs written by a3_clinical_calibration.py (figures/plotdata/), identical numbers.
suppressMessages({ library(ggplot2); library(patchwork); library(dplyr) })
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "scripts/R"
source(file.path(HERE, "theme_a3.R"))

BASE <- normalizePath(file.path(HERE, "..", ".."))
PD   <- file.path(BASE, "figures", "plotdata")
FIG  <- file.path(BASE, "figures")
rd   <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)

# ---------------- (a) reliability diagram ----------------
rel <- rd("F7_reliability.csv") %>%
  mutate(which = factor(which, levels = c("raw", "calibrated"),
                        labels = c("raw score", "Platt-calibrated")))
pa <- ggplot(rel, aes(mean_pred, obs_rate, colour = which)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4, colour = "grey55", linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  geom_point(aes(size = n)) +
  scale_colour_manual(values = c("raw score" = oi("verm"), "Platt-calibrated" = oi("blue"))) +
  scale_size_continuous(range = c(1.0, 3.6), guide = "none") +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "predicted P(pathogenic)", y = "observed pathogenic rate",
       title = "pooled ClinVar, ESM-2 650M (held-out)") +
  theme_a3() + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))

# ---------------- (b) coverage / abstention ----------------
cov <- rd("F7_coverage.csv")
covl <- bind_rows(
  data.frame(coverage = cov$coverage, value = cov$accuracy, metric = "overall accuracy"),
  data.frame(coverage = cov$coverage, value = cov$pathogenic_precision, metric = "pathogenic-call precision")
) %>% filter(!is.na(value))
pb <- ggplot(covl, aes(coverage, value, colour = metric)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
  scale_colour_manual(values = c("overall accuracy" = oi("blue"),
                                 "pathogenic-call precision" = oi("orange"))) +
  scale_x_continuous(limits = c(0.05, 1.0), breaks = seq(0.2, 1.0, 0.2)) +
  scale_y_continuous(limits = c(0.4, 0.9)) +
  labs(x = "coverage (retained most-confident fraction)", y = "accuracy / precision",
       title = "selective abstention (|P-0.5| ranking)") +
  theme_a3() + theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0))

p <- pa + pb + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", family = A3_FONT))
a3_save(p, "F4_clinical_calibration", 9.2, 4.0, FIG)
