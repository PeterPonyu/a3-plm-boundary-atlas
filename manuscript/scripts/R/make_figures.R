#!/usr/bin/env Rscript
# A3 main figures F1-F5 (PLM mutation-effect boundary atlas), ggplot2 + Times New
# Roman (Arial). Reads tidy CSVs from export_plot_data.py (figures/plotdata/),
# which reshape the committed result artifacts — identical numbers.
suppressMessages({ library(ggplot2); library(patchwork); library(dplyr) })
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "scripts/R"
source(file.path(HERE, "theme_a3.R"))

BASE <- normalizePath(file.path(HERE, "..", ".."))
PD   <- file.path(BASE, "figures", "plotdata")
FIG  <- file.path(BASE, "figures")
rd   <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)

# ---------------- F1 scale plateau (N1) ----------------
f1 <- rd("F1_scale.csv") %>% mutate(model_size = factor(model_size, c("35M", "150M", "650M", "3B")))
p1 <- ggplot(f1, aes(model_size, protein_mean_auroc, group = 1)) +
  geom_line(linewidth = 0.7, colour = oi("blue")) +
  geom_point(size = 2.6, colour = oi("blue")) +
  geom_text(aes(label = sprintf("%.3f", protein_mean_auroc)), vjust = -1.0, size = 3, family = A3_FONT) +
  annotate("rect", xmin = 2.6, xmax = 4.4, ymin = 0.718, ymax = 0.732,
           fill = oi("orange"), alpha = 0.12) +
  scale_y_continuous(limits = c(0.64, 0.745)) +
  labs(x = "ESM-2 parameters",
       y = "protein-mean AUROC (raw-score, variant effect, n = 69)") +
  theme_a3()
a3_save(p1, "F1_scale_plateau", 5.6, 4.0, FIG)

# ---------------- F2 cross-MLM-encoder agreement (N5) ----------------
f2 <- rd("F2_crossmodel.csv") %>% mutate(model = factor(model, model))
ag <- rd("F2_agreement.csv")
sp <- ag$spearman[ag$comparison == "ProtBERT vs ESM-2 (masked)"]
p2 <- ggplot(f2, aes(model, protein_mean_auroc, fill = model)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.3f", protein_mean_auroc)), vjust = -0.5, size = 3, family = A3_FONT) +
  scale_fill_manual(values = c(oi("orange"), oi("blue"), oi("sky")), guide = "none") +
  scale_y_continuous(limits = c(0, 0.8), expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "protein-mean AUROC") +
  theme_a3() + theme(axis.text.x = element_text(angle = 12, hjust = 1),
                     panel.grid.major.x = element_blank())
a3_save(p2, "F2_cross_model", 5.6, 4.0, FIG)

# ---------------- F3 clinical per-gene forest (N6c) ----------------
f3 <- rd("F3_clinvar.csv") %>% mutate(gene = factor(gene, levels = gene))
p3 <- ggplot(f3, aes(auroc_clinvar, gene)) +
  geom_vline(xintercept = 0.5, linewidth = 0.4, colour = "grey50", linetype = "dashed") +
  geom_vline(xintercept = 0.7, linewidth = 0.4, colour = oi("green"), linetype = "dotted") +
  geom_segment(aes(x = 0.5, xend = auroc_clinvar, yend = gene), colour = "grey80", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.22, colour = oi("blue"), linewidth = 0.45) +
  geom_point(aes(size = n_path), colour = oi("blue")) +
  scale_size_continuous(range = c(1.2, 4), name = "n pathogenic") +
  scale_x_continuous(limits = c(0.45, 1.0)) +
  labs(x = "per-gene ClinVar AUROC (pathogenic vs benign)", y = NULL) +
  theme_a3() + theme(panel.grid.major.y = element_blank(),
                     legend.position = c(0.98, 0.03), legend.justification = c(1, 0))
a3_save(p3, "F3_clinvar_forest", 5.8, 5.2, FIG)

# ---------------- F4 binding boundary by subset (N4) ----------------
f4 <- rd("F4_binding.csv") %>% mutate(subset = factor(subset, levels = rev(subset)),
                                      hl = key %in% c("antibody_ABAG", "noninterface"))
p4 <- ggplot(f4, aes(spearman, subset, colour = hl)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "black") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18, linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.2f", spearman)), vjust = -0.9, size = 2.7, family = A3_FONT, colour = "grey20") +
  scale_colour_manual(values = c(`FALSE` = oi("blue"), `TRUE` = oi("verm")), guide = "none") +
  labs(x = "Spearman(ESM-2 LLR, \u2212\u0394\u0394G binding)", y = NULL) +
  theme_a3() + theme(panel.grid.major.y = element_blank())
a3_save(p4, "F7_binding_boundary", 6.0, 3.6, FIG)

# ---------------- F5 developability map (N9) ----------------
f5 <- rd("F5_developability.csv") %>%
  mutate(property = factor(property, levels = property),
         lab = paste0(n_sig_raw_bh, "/", n_datasets),
         null_liability = property %in% c("polyreactivity", "aggregation"))
# Zero-height nulls must still show colour: plot a short verm stub + count label.
f5 <- f5 %>% mutate(frac_plot = ifelse(frac_sig <= 0, 0.035, frac_sig),
                    lab_vjust = ifelse(frac_sig <= 0, -0.2, -0.5))
p5 <- ggplot(f5, aes(property, frac_plot, fill = null_liability)) +
  geom_col(width = 0.66, colour = NA) +
  geom_text(aes(label = lab, y = frac_plot, vjust = lab_vjust), size = 2.8, family = A3_FONT) +
  scale_fill_manual(values = c(`FALSE` = oi("blue"), `TRUE` = oi("verm")),
                    labels = c("other", "liability (no BH-sig)"), name = NULL) +
  scale_y_continuous(limits = c(0, 1.08), expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "fraction of datasets with raw signal (BH q<0.05)") +
  theme_a3() + theme(axis.text.x = element_text(angle = 25, hjust = 1),
                     panel.grid.major.x = element_blank(),
                     legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
a3_save(p5, "F8_developability_map", 6.4, 4.0, FIG)

cat("A3 make_figures.R: done\n")
