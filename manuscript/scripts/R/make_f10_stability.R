#!/usr/bin/env Rscript
# A3 Figure F10: protein stability ddG zero-shot -- a pre-registered controlled
# negative. ESM-2 zero-shot ddG (Megascale, homology-clustered hold-out).
# Panel A: discrimination Spearman, train (in-distribution) vs homology-held-out
#   -- they transfer ~1:1 (~0.38 vs ~0.34), i.e. NO reliability gap. A dashed
#   reference line at the supervised ThermoMPNN literature level (~0.65) shows the
#   ESM zero-shot signal is WEAK though stable.
# Panel B: calibration RMSE, train vs homology-held-out (~equal, inflation ~+0.04)
#   = no homology calibration gap. ggplot2 + Times New Roman (Arial).
# Honest: signal is weak everywhere; the "no gap" is partly because a foundation
# model's pretraining is not truly out-of-distribution for these domains.
suppressMessages({ library(ggplot2); library(patchwork); library(dplyr); library(tidyr) })
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "scripts/R"
source(file.path(HERE, "theme_a3.R"))

BASE <- normalizePath(file.path(HERE, "..", ".."))
PD   <- file.path(BASE, "figures", "plotdata")
FIG  <- file.path(BASE, "figures")

rd <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)

SPLITS <- c("train", "homology-held-out")
SPLITCOL <- c("train" = oi("blue"), "homology-held-out" = oi("verm"))

# ---- Panel A: discrimination Spearman, train vs homology-held-out ----
a <- rd("F10_discrimination.csv")
a$split <- factor(a$split, levels = SPLITS)
a$variant <- factor(a$variant, levels = c("v1 (orig)", "v2 WT", "v2 masked"))
ref <- rd("F10_reference.csv")
sref <- ref$spearman[1]

pA <- ggplot(a, aes(variant, spearman, fill = split)) +
  geom_hline(yintercept = sref, linewidth = 0.4, linetype = 2, colour = "grey45") +
  geom_col(position = position_dodge(width = 0.7), width = 0.62, colour = NA) +
  geom_text(aes(label = sprintf("%.2f", spearman)),
            position = position_dodge(width = 0.7), vjust = -0.5,
            size = 2.6, family = A3_FONT, colour = "grey20") +
  annotate("text", x = 0.55, y = sref + 0.025, hjust = 0,
           label = "supervised reference (ThermoMPNN, literature)",
           size = 2.7, colour = "grey35", family = A3_FONT) +
  scale_fill_manual(values = SPLITCOL) +
  scale_y_continuous(limits = c(0, 0.78), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "discrimination Spearman",
       title = "Discrimination transfers ~1:1 (no reliability gap), but weak vs supervised") +
  theme_a3() +
  theme(legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0),
        panel.grid.major.x = element_blank())

# ---- Panel B: calibration RMSE, train vs homology-held-out ----
b <- rd("F10_calibration.csv")
b$split <- factor(b$split, levels = SPLITS)
b$variant <- factor(b$variant, levels = c("v1 (orig)", "v2 WT", "v2 masked"))
# one inflation value per variant (heldout - heldin)
infl <- b %>% distinct(variant, inflation)
infl_lab <- paste0(infl$variant, ": ", sprintf("%+.3f", infl$inflation))

pB <- ggplot(b, aes(variant, calib_rmse, fill = split)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62, colour = NA) +
  geom_text(aes(label = sprintf("%.3f", calib_rmse)),
            position = position_dodge(width = 0.7), vjust = -0.5,
            size = 2.6, family = A3_FONT, colour = "grey20") +
  scale_fill_manual(values = SPLITCOL, guide = "none") +
  scale_y_continuous(limits = c(0, 1.12), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "RMSE (MegaScale ddG_ML units)",
       title = "Calibration RMSE ~equal across the homology split (inflation ~+0.04)") +
  theme_a3() +
  theme(panel.grid.major.x = element_blank())

fig <- (pA | pB) + plot_annotation(tag_levels = "A",
  theme = theme(plot.tag = element_text(size = 13, face = "bold", family = A3_FONT)))

a3_save(fig, "F5_stability_ddg", 9.2, 4.2, FIG)
