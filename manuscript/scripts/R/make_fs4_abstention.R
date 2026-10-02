#!/usr/bin/env Rscript
# A3 Supplementary Fig FS4 — ESM-2's clamped confidence is anti-helpful for selective
# abstention, across model tiers. Panel A: per tier (35M/150M/650M/ESM-1v) protein-mean
# AUROC at full (100%) vs keep-highest-confidence-50% vs shuffle-50% control — keep-hi
# falls below BOTH full and shuffle in 4/4 tiers (confidence does NOT select reliable
# predictions). Panel B: 650M coverage-risk curve, confidence-ranked abstention vs the
# shuffled-confidence control (flat). Honest negative. ggplot2 + Times New Roman.
suppressMessages({ library(ggplot2); library(patchwork); library(dplyr); library(tidyr) })
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "scripts/R"
source(file.path(HERE, "theme_a3.R"))
BASE <- normalizePath(file.path(HERE, "..", ".."))
PD <- file.path(BASE, "figures", "plotdata"); FIG <- file.path(BASE, "figures")
rd <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)

# ---- Panel A: per-tier full / keep-hi-50% / shuffle-50% ----
a <- rd("FS4_abstention_tiers.csv")
tier_lv <- a$tier[order(a$order)] |> unique()
cond_lv <- c("full (100%)", "keep-hi 50%", "shuffle 50%")
a <- a %>% mutate(tier = factor(tier, levels = tier_lv),
                  condition = factor(condition, levels = cond_lv))
condcol <- c("full (100%)" = oi("grey"),
             "keep-hi 50%" = oi("verm"),
             "shuffle 50%" = oi("blue"))
pA <- ggplot(a, aes(tier, auroc, fill = condition)) +
  geom_col(position = position_dodge(0.74), width = 0.68) +
  geom_text(aes(label = sprintf("%.3f", auroc)), position = position_dodge(0.74),
            vjust = -0.4, size = 2.3, family = A3_FONT, colour = "grey20") +
  scale_fill_manual(values = condcol, name = NULL) +
  scale_y_continuous(limits = c(0, 0.86), breaks = seq(0, 0.8, 0.2),
                     expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "protein-mean AUROC (raw score)", tag = "A",
       title = "keep-highest-confidence 50% < full and < shuffle (4/4 tiers)") +
  theme_a3() +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "top", legend.direction = "horizontal")

# ---- Panel B: 650M coverage-risk curve ----
b <- rd("FS4_abstention_coverage.csv") %>%
  pivot_longer(c(conf_auroc, shuffle_auroc), names_to = "series", values_to = "auroc") %>%
  mutate(series = recode(series,
                         conf_auroc = "confidence-ranked",
                         shuffle_auroc = "shuffled control"))
sercol <- c("confidence-ranked" = oi("verm"), "shuffled control" = oi("blue"))
pB <- ggplot(b, aes(coverage, auroc, colour = series)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_colour_manual(values = sercol, name = NULL) +
  scale_x_reverse(breaks = seq(1.0, 0.5, -0.1)) +
  scale_y_continuous(limits = c(0.62, 0.72), breaks = seq(0.62, 0.72, 0.02)) +
  labs(x = "coverage (fraction of highest-confidence variants kept)",
       y = "protein-mean AUROC (raw score)", tag = "B",
       title = "ESM-2 650M: trusting confidence degrades discrimination") +
  theme_a3() +
  theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0))

a3_save(pA / pB, "FS4_abstention_tiers", 6.8, 6.6, FIG)
cat("A3 make_fs4_abstention.R: done\n")
