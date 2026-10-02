#!/usr/bin/env Rscript
# A3 schematics in B1 style (Times/Tinos, Okabe-Ito): research design, route, architecture.
#   a3_research_design        — study logic / boundary-atlas claim chain
#   a3_research_route         — the N-series roadmap (N1->N5->N6->N4->N10->N9->atlas)
#   a3_research_architecture  — data sources -> ESM-2 scorer -> metrics -> atlas/figures -> manuscript
suppressMessages({
  source("/home/zeyufu/Desktop/labs/active/reliability-program/biochemical/biochem-research/paper/B1-cofolding-confidence/scripts/R/theme_b1.R")
  library(ggplot2)
})
OUT <- "/home/zeyufu/Desktop/labs/active/reliability-program/biochemical/biochem-research/research/results/a3_meta"

boxdf <- function(x, y, hw, hh, label, fill, textcol = "black", size = 3.0, bold = FALSE)
  data.frame(x = x, y = y, hw = hw, hh = hh, label = label, fill = fill,
             textcol = textcol, size = size, bold = bold, stringsAsFactors = FALSE)
draw_boxes <- function(p, B) p +
  geom_rect(data = B, aes(xmin = x - hw, xmax = x + hw, ymin = y - hh, ymax = y + hh),
            fill = B$fill, colour = "grey30", linewidth = 0.4, inherit.aes = FALSE) +
  geom_text(data = B, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            colour = B$textcol, size = B$size, family = B1_FONT, lineheight = 0.92,
            fontface = ifelse(B$bold, "bold", "plain"))
arrow_seg <- function(x1, y1, x2, y2, col = "grey25", lwd = 0.5)
  annotate("segment", x = x1, y = y1, xend = x2, yend = y2, colour = col, linewidth = lwd,
           arrow = grid::arrow(length = unit(0.10, "inches"), type = "closed"))
canvas <- function(xl, yl) ggplot() + coord_fixed(1, xlim = xl, ylim = yl, expand = FALSE) +
  theme_void(base_family = B1_FONT) + theme(plot.margin = margin(2, 2, 2, 2))
C <- list(grey = "#E6E6E6", lblue = "#D6E8F5", lgreen = "#D7EFE4", lorange = "#FBE7C6",
          lpurple = "#EFD9E6", lverm = "#F6D9CC")

############################ 1. A3 RESEARCH DESIGN ############################
B <- rbind(
  boxdf(50, 95, 46, 4.5, "PREMISE: protein language models (ESM-2) are increasingly used as zero-shot variant-effect oracles", C$grey, size = 2.7, bold = TRUE),
  boxdf(50, 83, 34, 5.5, "QUESTION\nWhere does zero-shot LLR reliability HOLD,\nand where does it BREAK?", C$lblue, size = 3.0, bold = TRUE),
  boxdf(26, 66, 21, 5.5, "RELIABLE axis\n(naturalness / conservation)", C$lgreen, size = 2.8, bold = TRUE),
  boxdf(26, 55, 21, 5, "DMS fitness AUROC ~0.72\n(scale ceiling 650M; N1)", C$lgreen, size = 2.6),
  boxdf(26, 45.5, 21, 5, "ClinVar path/benign 0.78;\nencoder-robust ProtBERT (N5,N6c)", C$lgreen, size = 2.6),
  boxdf(74, 66, 21, 5.5, "BOUNDARY axis\n(physical / biophysical)", C$lverm, size = 2.8, bold = TRUE),
  boxdf(74, 55, 21, 5, "Binding ΔΔG weak (rho 0.20),\nantibody ~0.08 (N4)", C$lverm, size = 2.6),
  boxdf(74, 45.5, 21, 5, "Stability ΔΔG weak 0.36 (N10);\ndevelopability null 0/14 (N9)", C$lverm, size = 2.6),
  boxdf(50, 30, 46, 5, "DE-CIRCULARIZATION CONTROLS: cross-encoder applicability gate - homology-held-out split - germline control (withdrawn where inapplicable)", C$lpurple, size = 2.5, bold = TRUE),
  boxdf(50, 16, 46, 5, "BOUNDED ATLAS: ESM-2 ranks naturalness/pathogenicity moderately (never reliably, <0.9)\nand cannot read physical binding. Retrospective, raw-score, no clinical deployment claim.", C$lorange, size = 2.6, bold = TRUE)
)
p1 <- draw_boxes(canvas(c(0, 100), c(8, 100)), B) +
  arrow_seg(50, 89.5, 50, 88.5) + arrow_seg(40, 80, 30, 71.5) + arrow_seg(60, 80, 70, 71.5) +
  arrow_seg(26, 42.9, 40, 35) + arrow_seg(74, 42.9, 60, 35) + arrow_seg(50, 25, 50, 21)
b1_save(p1, "a3_research_design", w = 7.6, h = 7.0, figdir = OUT)

############################ 2. A3 RESEARCH ROUTE ############################
ph <- function(x, y, lab, sub) boxdf(x, y, 13.5, 8.5, paste0(lab, "\n\n", sub), C$lgreen, size = 2.5, bold = TRUE)
B2 <- rbind(
  ph(16, 72, "N1 Scale", "DMS AUROC vs\nmodel size (ceiling)"),
  ph(44, 72, "N5 Architecture", "ProtBERT replicates\n(encoder-robust)"),
  ph(72, 72, "N6c Clinical", "ClinVar transfer\nper gene"),
  ph(16, 44, "N4 Binding", "SKEMPI ΔΔG\n(interface collapse)"),
  ph(44, 44, "N10 Stability", "MegaScale ΔΔG\n+ homology split"),
  ph(72, 44, "N9 Developability", "FLAb nulls\n(polyreact/aggreg)")
)
B2b <- boxdf(50, 18, 30, 5.5, "BOUNDARY ATLAS\n(unifying reliability map)", C$lorange, size = 2.8, bold = TRUE)
p2 <- draw_boxes(draw_boxes(canvas(c(0, 100), c(8, 86)), B2), B2b) +
  arrow_seg(29.5, 72, 30.5, 72) + arrow_seg(57.5, 72, 58.5, 72) +
  arrow_seg(72, 63.5, 30, 52.5, col = "grey45") +
  arrow_seg(29.5, 44, 30.5, 44) + arrow_seg(57.5, 44, 58.5, 44) +
  arrow_seg(72, 35.5, 55, 23) + arrow_seg(44, 35.5, 50, 23.5) + arrow_seg(16, 35.5, 45, 23)
b1_save(p2, "a3_research_route", w = 8.2, h = 5.2, figdir = OUT)

############################ 3. A3 RESEARCH ARCHITECTURE ############################
B3 <- rbind(
  boxdf(12, 90, 11, 6, "ProteinGym\nDMS\n69 prot", C$lblue, size = 2.3),
  boxdf(31, 90, 11, 6, "ClinVar +\ngnomAD\n20 genes", C$lblue, size = 2.3),
  boxdf(50, 90, 11, 6, "SKEMPI v2\nbinding ΔΔG\n3,424", C$lblue, size = 2.3),
  boxdf(69, 90, 11, 6, "MegaScale\nstability ΔΔG\n55k held-out", C$lblue, size = 2.3),
  boxdf(88, 90, 11, 6, "FLAb\ndevelopability\n63 sets", C$lblue, size = 2.3),
  boxdf(50, 72, 46, 5.5, "ESM-2 SCORER  raw LLR (WT- / masked-marginal), no clamping\n+ ProtBERT cross-encoder (N5) - one deterministic re-score discipline", C$lorange, size = 2.5, bold = TRUE),
  boxdf(22, 55, 17, 5, "AUROC\n(discrimination)", C$lgreen, size = 2.5),
  boxdf(50, 55, 17, 5, "Spearman rho\n(ΔΔG tasks)", C$lgreen, size = 2.5),
  boxdf(78, 55, 17, 5, "calibration RMSE\n+ bootstrap CI", C$lgreen, size = 2.5),
  boxdf(30, 38, 20, 5, "Boundary atlas\n(this work)", C$lpurple, size = 2.5),
  boxdf(70, 38, 20, 5, "Per-axis figures\nF1-F5, F10", C$lpurple, size = 2.5),
  boxdf(50, 22, 30, 5, "A3 manuscript (R/Times)", C$grey, size = 2.7, bold = TRUE)
)
p3 <- draw_boxes(canvas(c(0, 100), c(15, 98)), B3) +
  arrow_seg(12, 84, 40, 77.5) + arrow_seg(31, 84, 45, 77.5) + arrow_seg(50, 84, 50, 77.5) +
  arrow_seg(69, 84, 56, 77.5) + arrow_seg(88, 84, 60, 77.5) +
  arrow_seg(40, 66.5, 24, 60) + arrow_seg(50, 66.5, 50, 60) + arrow_seg(60, 66.5, 76, 60) +
  arrow_seg(30, 50, 30, 43) + arrow_seg(70, 50, 70, 43) +
  arrow_seg(30, 35.5, 45, 27) + arrow_seg(70, 35.5, 55, 27)
b1_save(p3, "a3_research_architecture", w = 8.6, h = 6.6, figdir = OUT)
