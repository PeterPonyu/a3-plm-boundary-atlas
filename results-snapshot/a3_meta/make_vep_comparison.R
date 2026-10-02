#!/usr/bin/env Rscript
# Where ESM-1v zero-shot sits among published variant-effect predictors on DMS.
# Per-protein |Spearman(predictor, DMS)|, median across the Livesey & Marsh set; ESM highlighted.
suppressMessages({
  source("/home/zeyufu/Desktop/labs/active/reliability-program/biochemical/biochem-research/paper/B1-cofolding-confidence/scripts/R/theme_b1.R")
  library(ggplot2)
})
OUT <- "/home/zeyufu/Desktop/labs/active/reliability-program/biochemical/biochem-research/research/results/a3_meta"
d <- read.csv(file.path(OUT, "vep_comparison.csv"), stringsAsFactors = FALSE)
d <- d[order(d$median_abs_spearman), ]
d$predictor <- factor(d$predictor, levels = d$predictor)
d$is_esm <- d$predictor == "ESM-1v"

p <- ggplot(d, aes(x = median_abs_spearman, y = predictor, fill = is_esm)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", median_abs_spearman)), hjust = -0.2, size = 2.8, family = B1_FONT) +
  scale_fill_manual(values = c(`TRUE` = oi("blue"), `FALSE` = "grey70"), guide = "none") +
  scale_x_continuous(limits = c(0, 0.58), breaks = seq(0, 0.5, 0.1)) +
  labs(x = expression("median |Spearman(predictor, DMS)| across proteins"), y = NULL) +
  theme_b1(10) + theme(panel.grid.major.y = element_blank())

b1_save(p, "vep_comparison", w = 5.6, h = 3.0, figdir = OUT)
