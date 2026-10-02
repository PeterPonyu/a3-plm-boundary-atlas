#!/usr/bin/env Rscript
# A3 graphical abstract / TOC graphic: a single conceptual panel of the reliability boundary.
# Honest by construction — qualitative zones with each task's ACTUAL native metric labelled
# (AUROC for discrimination, Spearman for physical), NOT two metrics conflated on one numeric axis.
suppressMessages({
  source("/home/zeyufu/Desktop/labs/active/reliability-program/biochemical/biochem-research/paper/B1-cofolding-confidence/scripts/R/theme_b1.R")
  library(ggplot2)
})
OUT <- "/home/zeyufu/Desktop/labs/active/reliability-program/biochemical/biochem-research/research/results/a3_meta"
oiv <- function(k) oi(k)

p <- ggplot() + coord_fixed(1, xlim = c(0, 100), ylim = c(0, 52), expand = FALSE) +
  theme_void(base_family = B1_FONT) + theme(plot.margin = margin(4, 4, 4, 4))

# reliability spectrum strip (green -> red), boundary in the middle
strip <- data.frame(x = seq(8, 92, length.out = 200))
p <- p +
  geom_tile(data = strip, aes(x = x, y = 26, fill = x), height = 5, width = (92 - 8) / 199) +
  scale_fill_gradient2(low = oiv("green"), mid = "grey85", high = oiv("verm"),
                       midpoint = 50, guide = "none") +
  annotate("segment", x = 8, xend = 92, y = 26, yend = 26, colour = NA) +
  annotate("text", x = 50, y = 21, label = "ESM-2 zero-shot likelihood: one model, two regimes",
           family = B1_FONT, size = 3.2, fontface = "bold", colour = "grey25")

# left zone: RELIABLE
p <- p +
  annotate("rect", xmin = 6, xmax = 48, ymin = 33, ymax = 49, fill = scales::alpha(oiv("green"), 0.16), colour = oiv("green"), linewidth = 0.5) +
  annotate("text", x = 27, y = 46, label = "RELIABLE", family = B1_FONT, size = 4.2, fontface = "bold", colour = oiv("green")) +
  annotate("text", x = 27, y = 42.3, label = "naturalness · conservation", family = B1_FONT, size = 2.8, colour = "grey25") +
  annotate("text", x = 27, y = 38.3, label = "DMS fitness   AUROC 0.72", family = B1_FONT, size = 2.9, colour = "grey15") +
  annotate("text", x = 27, y = 35.0, label = "ClinVar path/benign   AUROC 0.78", family = B1_FONT, size = 2.9, colour = "grey15")

# right zone: BLIND
p <- p +
  annotate("rect", xmin = 52, xmax = 94, ymin = 33, ymax = 49, fill = scales::alpha(oiv("verm"), 0.14), colour = oiv("verm"), linewidth = 0.5) +
  annotate("text", x = 73, y = 46, label = "BLIND", family = B1_FONT, size = 4.2, fontface = "bold", colour = oiv("verm")) +
  annotate("text", x = 73, y = 42.3, label = "physical interaction · liability", family = B1_FONT, size = 2.8, colour = "grey25") +
  annotate("text", x = 73, y = 38.3, label = "binding ΔΔG   Spearman 0.08–0.27", family = B1_FONT, size = 2.9, colour = "grey15") +
  annotate("text", x = 73, y = 35.0, label = "antibody / developability   ≈ null", family = B1_FONT, size = 2.9, colour = "grey15")

# boundary marker + connectors
p <- p +
  annotate("segment", x = 50, xend = 50, y = 23.5, yend = 29.3, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  annotate("text", x = 50, y = 30.6, label = "reliability boundary", family = B1_FONT, size = 2.5, colour = "grey40") +
  annotate("segment", x = 27, xend = 22, y = 33, yend = 28.7, colour = oiv("green"), linewidth = 0.4,
           arrow = grid::arrow(length = unit(0.07, "inches"), type = "closed")) +
  annotate("segment", x = 73, xend = 78, y = 33, yend = 28.7, colour = oiv("verm"), linewidth = 0.4,
           arrow = grid::arrow(length = unit(0.07, "inches"), type = "closed"))

# bottom takeaway
p <- p +
  annotate("text", x = 50, y = 12.5, label = "A uniform, self-de-circularizing audit maps where to trust the score — and where to abstain.",
           family = B1_FONT, size = 2.9, colour = "grey20") +
  annotate("text", x = 50, y = 7.5, label = "6 axes · DMS · ClinVar · SKEMPI · MegaScale · FLAb · raw scores, bootstrap CIs",
           family = B1_FONT, size = 2.5, colour = "grey45")

b1_save(p, "a3_graphical_abstract", w = 7.4, h = 3.9, figdir = OUT)
