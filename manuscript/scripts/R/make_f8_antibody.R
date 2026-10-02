#!/usr/bin/env Rscript
# A3 Fig 8 — antibody-specific LM (AntiBERTy) head-to-head vs ESM-2 on the antibody
# axes. Two panels: (a) N9 developability liability nulls (BH-significant dataset
# count); (b) N4 antibody-antigen binding ddG correlation. ggplot2 + Times New Roman.
suppressMessages({ library(ggplot2); library(patchwork); library(dplyr); library(tidyr) })
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "scripts/R"
source(file.path(HERE, "theme_a3.R"))
BASE <- normalizePath(file.path(HERE, "..", ".."))
PD <- file.path(BASE, "figures", "plotdata"); FIG <- file.path(BASE, "figures")
rd <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)

# Panel A: N9 liability nulls — BH-significant dataset fraction, AntiBERTy vs ESM-2
n9 <- rd("F8_antibody_n9.csv") %>%
  pivot_longer(c(antiberty_bh_sig, esm2_bh_sig), names_to = "model", values_to = "sig") %>%
  mutate(model = recode(model, antiberty_bh_sig = "AntiBERTy", esm2_bh_sig = "ESM-2 650M"),
         frac = sig / n_datasets,
         lab = paste0(sig, "/", n_datasets),
         property = factor(property, c("polyreactivity", "aggregation")))
pA <- ggplot(n9, aes(property, frac, fill = model)) +
  geom_col(position = position_dodge(0.72), width = 0.66) +
  geom_text(aes(label = lab), position = position_dodge(0.72), vjust = -0.4, size = 2.8, family = A3_FONT) +
  scale_fill_manual(values = c("AntiBERTy" = oi("verm"), "ESM-2 650M" = oi("blue")), name = NULL) +
  scale_y_continuous(limits = c(0, 0.32), expand = expansion(c(0, 0.04))) +
  labs(x = NULL, y = "fraction of datasets significant (BH q<0.05)", tag = "A",
       title = "developability liabilities") +
  theme_a3() + theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
                     axis.title.y = element_text(size = 8.5),
                     plot.tag.position = c(0.015, 0.985),
                     plot.margin = margin(5, 6, 4, 16),
                     panel.grid.major.x = element_blank())

# Panel B: N4 antibody-antigen binding ddG correlation, AntiBERTy vs ESM-2
n4 <- rd("F8_antibody_n4.csv") %>% mutate(model = factor(model, levels = rev(model)))
pB <- ggplot(n4, aes(rho, model, colour = model)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "black") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.16, linewidth = 0.7) +
  geom_point(size = 2.8) +
  geom_text(aes(label = sprintf("%.3f", rho)), vjust = -1.1, size = 2.8, family = A3_FONT, colour = "grey20") +
  scale_colour_manual(values = c("AntiBERTy" = oi("verm"), "ESM-2 650M" = oi("blue")), guide = "none") +
  scale_x_continuous(limits = c(-0.1, 0.2)) +
  labs(x = "Spearman(LLR, \u2212\u0394\u0394G binding), n=627 antibody muts", y = NULL, tag = "B",
       title = "antibody\u2013antigen binding") +
  theme_a3() + theme(panel.grid.major.y = element_blank(),
                     axis.title.x = element_text(size = 8.5),
                     plot.margin = margin(5, 10, 4, 5))

a3_save(pA + pB, "F9_antibody_lm", 8.4, 3.6, FIG)
cat("A3 make_f8_antibody.R: done\n")
