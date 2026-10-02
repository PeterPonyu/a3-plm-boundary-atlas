#!/usr/bin/env Rscript
# A3 Figure F9: protein uncertainty quantification under distribution shift --
# the missing UQ axis. Four UQ methods (BayesianRidge baseline, Deep Ensemble,
# GP-RBF, MC-Dropout) on identical ESM-2 650M embeddings + identical FLIP splits.
# Panel A: rho_unc (predicted-uncertainty vs actual-error rank correlation) by
#   method x split, grouped by shift severity; learned uncertainty (ensemble /
#   MC-dropout) recovers weak signal on MODERATE shift but collapses on the
#   EXTREME split, while the simple baseline is ~0 everywhere.
# Panel B: retrospective Bayesian-optimization acquisition (UCB vs greedy vs
#   random) final-best per split; UCB never beats greedy => uncertainty-driven
#   acquisition is null. ggplot2 + Times New Roman (Arial).
suppressMessages({ library(ggplot2); library(patchwork); library(dplyr); library(tidyr) })
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "scripts/R"
source(file.path(HERE, "theme_a3.R"))

BASE <- normalizePath(file.path(HERE, "..", ".."))
PD   <- file.path(BASE, "figures", "plotdata")
FIG  <- file.path(BASE, "figures")

rd <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)

METHODS <- c("BayesianRidge", "Deep Ensemble", "GP-RBF", "MC-Dropout")
METHCOL <- c("BayesianRidge" = oi("grey"), "Deep Ensemble" = oi("blue"),
             "GP-RBF" = oi("orange"), "MC-Dropout" = oi("green"))

# ---- Panel A: rho_unc by method x split, grouped by shift severity ----
a <- rd("F9_rho_unc.csv")
a$method <- factor(a$method, levels = METHODS)
# short split labels
a$split_lab <- sub("^gb1/", "gb1: ", a$split)
a$split_lab <- sub("^meltome/", "melt: ", a$split_lab)
# order severity: extreme -> moderate -> in-distribution
a$shift_class <- factor(a$shift_class,
                        levels = c("extreme", "moderate-shift", "in-distribution"),
                        labels = c("extreme shift", "moderate shift", "in-distribution"))
# order splits within facets to keep extreme split first
split_ord <- a %>% distinct(split, split_lab, shift_class) %>%
  arrange(shift_class, split) %>% pull(split_lab)
a$split_lab <- factor(a$split_lab, levels = unique(split_ord))

pA <- ggplot(a, aes(split_lab, rho_unc, fill = method)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, colour = NA) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                position = position_dodge(width = 0.8), width = 0.3,
                linewidth = 0.3, colour = "grey30") +
  facet_grid(. ~ shift_class, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = METHCOL) +
  labs(x = NULL,
       y = expression("uncertainty-vs-error " * rho[unc]),
       title = "Learned UQ recovers weak signal on moderate shift, collapses on extreme") +
  theme_a3() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major.x = element_blank())

# ---- Panel B: BO acquisition (UCB vs greedy vs random) final-best per split ----
b <- rd("F9_bo_acquisition.csv")
bl <- b %>%
  pivot_longer(c(greedy, ucb, random), names_to = "acq", values_to = "best") %>%
  mutate(split_lab = sub("^gb1/", "gb1: ", split),
         split_lab = sub("^meltome/", "melt: ", split_lab),
         acq = factor(acq, levels = c("greedy", "ucb", "random"),
                      labels = c("greedy", "UCB", "random")))
bl$split_lab <- factor(bl$split_lab, levels = unique(split_ord))

ACQCOL <- c("greedy" = oi("black"), "UCB" = oi("blue"), "random" = oi("grey"))
pB <- ggplot(bl, aes(split_lab, best, colour = acq, group = acq)) +
  geom_line(linewidth = 0.4, alpha = 0.5,
            position = position_dodge(width = 0.35)) +
  geom_point(aes(shape = acq), size = 2.4,
             position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = ACQCOL) +
  scale_shape_manual(values = c("greedy" = 16, "UCB" = 17, "random" = 4)) +
  labs(x = NULL, y = "best label found (retrospective BO)",
       title = "Uncertainty-driven acquisition (UCB) never beats greedy") +
  theme_a3() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = c(0.5, 0.16),
        legend.direction = "horizontal")

fig <- pA / pB + plot_annotation(tag_levels = "A",
  theme = theme(plot.tag = element_text(size = 13, face = "bold", family = A3_FONT)))

a3_save(fig, "F6_protein_uq", 9.2, 7.4, FIG)
