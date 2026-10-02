#!/usr/bin/env Rscript
# A3 Supplementary Fig FS3 — antibody binding-boundary robustness. Grouped forest:
# per subset (overall / general / antibody / interface / non-interface) the
# Spearman(LLR, -ddG) for ESM-2 650M (WT-marginal) vs ESM-2 3B (WT-marginal) vs
# ESM-2 650M (masked-marginal), with bootstrap CIs. Key messages: 3B does NOT
# rescue (overall 0.148 < 650M 0.197) and masked approx wt (scoring-robust), so the
# binding boundary is not fixed by scaling or scoring choice. ggplot2 + Times New Roman.
suppressMessages({ library(ggplot2); library(dplyr) })
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "scripts/R"
source(file.path(HERE, "theme_a3.R"))
BASE <- normalizePath(file.path(HERE, "..", ".."))
PD <- file.path(BASE, "figures", "plotdata"); FIG <- file.path(BASE, "figures")
rd <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)

d <- rd("FS3_antibody_robustness.csv")
sub_lv <- d$subset[order(d$order)] |> unique()
meth_lv <- c("ESM-2 650M (WT)", "ESM-2 3B (WT)", "ESM-2 650M (masked)")
d <- d %>% mutate(
  subset = factor(subset, levels = rev(sub_lv)),   # rev: top subset highest on y axis
  method = factor(method, levels = meth_lv))

dodge <- position_dodge(width = 0.62)
cols <- c("ESM-2 650M (WT)" = oi("blue"),
          "ESM-2 3B (WT)" = oi("verm"),
          "ESM-2 650M (masked)" = oi("green"))

p <- ggplot(d, aes(rho, subset, colour = method, group = method)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "black") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.32, linewidth = 0.6,
                 position = dodge) +
  geom_point(size = 2.2, position = dodge) +
  scale_colour_manual(values = cols, name = NULL,
                      guide = guide_legend(reverse = FALSE)) +
  scale_x_continuous(limits = c(-0.15, 0.36), breaks = seq(-0.1, 0.3, 0.1)) +
  labs(x = expression(paste("Spearman(LLR, ", -Delta*Delta*"G binding)")),
       y = NULL,
       title = "3B does not rescue the binding boundary; masked \u2248 WT scoring") +
  theme_a3() +
  theme(panel.grid.major.y = element_blank(),
        legend.position = c(0.98, 0.02), legend.justification = c(1, 0))

a3_save(p, "FS3_antibody_robustness", 6.6, 3.8, FIG)
cat("A3 make_fs3_antibody.R: done\n")
