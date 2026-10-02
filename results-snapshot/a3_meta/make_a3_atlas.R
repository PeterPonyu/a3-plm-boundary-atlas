#!/usr/bin/env Rscript
# A3 reliability BOUNDARY ATLAS: one figure mapping where ESM-2 zero-shot LLR holds vs collapses.
#  A: discrimination tasks (AUROC) — DMS fitness + ClinVar, moderate (0.72-0.78), below reliable.
#  B: physical/biophysical tasks (Spearman) — binding/stability ΔΔG, weak; antibody & non-interface ~0.
# Within-study summaries from A3's committed result JSONs (N1/N5/N6c/N4/N10). NOT an independent-dataset
# pool (A3's tasks are heterogeneous) — it is a task-boundary map. B1 R/Times style.
suppressMessages({
  source("/home/zeyufu/Desktop/labs/active/reliability-program/biochemical/biochem-research/paper/B1-cofolding-confidence/scripts/R/theme_b1.R")
  library(ggplot2); library(patchwork)
})
OUT <- "/home/zeyufu/Desktop/labs/active/reliability-program/biochemical/biochem-research/research/results/a3_meta"
A <- read.csv(file.path(OUT, "a3_atlas_auroc.csv"), check.names = FALSE)
B <- read.csv(file.path(OUT, "a3_atlas_spearman.csv"), check.names = FALSE)
A$task <- factor(A$task, levels = rev(A$task))
B$task <- factor(B$task, levels = rev(B$task))
num <- function(x) suppressWarnings(as.numeric(x))

pA <- ggplot(A, aes(x = num(value), y = task)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0.9, linetype = "dotted", colour = oi("verm"), linewidth = 0.45) +
  geom_errorbar(aes(xmin = num(lo), xmax = num(hi)), orientation = "y", width = 0.25,
                linewidth = 0.5, colour = oi("green"), na.rm = TRUE) +
  geom_point(size = 2.8, colour = oi("green")) +
  annotate("text", x = 0.5, y = 0.55, label = "chance", size = 2.4, colour = "grey45",
           family = B1_FONT, hjust = -0.1) +
  annotate("text", x = 0.9, y = 0.55, label = "reliable", size = 2.4, colour = oi("verm"),
           family = B1_FONT, hjust = 1.1) +
  scale_x_continuous(limits = c(0.45, 1.0), breaks = seq(0.5, 1.0, 0.1)) +
  labs(x = "AUROC", y = NULL, tag = "A",
       title = "Discrimination tasks (naturalness / pathogenicity)") +
  theme_b1(10) + theme(plot.title = element_text(size = 8.5, hjust = 0.5),
                       panel.grid.major.y = element_blank())

pB <- ggplot(B, aes(x = num(value), y = task)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(xmin = num(lo), xmax = num(hi)), orientation = "y", width = 0.25,
                linewidth = 0.5, colour = oi("blue"), na.rm = TRUE) +
  geom_point(size = 2.8, colour = oi("blue")) +
  scale_x_continuous(limits = c(-0.15, 0.45), breaks = seq(0, 0.4, 0.1)) +
  labs(x = expression("Spearman " * rho), y = NULL, tag = "B",
       title = "Physical / biophysical tasks (binding, stability ΔΔG)") +
  theme_b1(10) + theme(plot.title = element_text(size = 8.5, hjust = 0.5),
                       panel.grid.major.y = element_blank())

p <- pA / pB + plot_layout(heights = c(1, 1.4))
b1_save(p, "a3_boundary_atlas", w = 6.8, h = 5.2, figdir = OUT)
