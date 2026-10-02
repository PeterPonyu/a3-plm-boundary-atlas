#!/usr/bin/env Rscript
# A3 Fig S2 — clinical-axis population control (N6c canonical, fp32 consistent scoring):
# per-gene ClinVar AUROC holds when the benign set is swapped to population-tolerant
# gnomAD (AC>=2). Shows the per-gene transfer is not a ClinVar-benign-ascertainment
# artifact. Uses N6c only (NOT the withdrawn N6/N6b long-protein-artifact numbers).
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) "."); if (length(HERE) == 0 || HERE == "") HERE <- "scripts/R"
source(file.path(HERE, "theme_a3.R"))
BASE <- normalizePath(file.path(HERE, "..", "..")); PD <- file.path(BASE, "figures", "plotdata"); FIG <- file.path(BASE, "figures")
rd <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)

g <- rd("FS2_clinvar_vs_gnomad.csv") %>% arrange(auroc_clinvar_benign) %>%
  mutate(gene = factor(gene, levels = unique(gene)))
gl <- g %>% select(gene, auroc_clinvar_benign, auroc_gnomad_benign) %>%
  pivot_longer(-gene, names_to = "benign_source", values_to = "auroc") %>%
  mutate(benign_source = recode(benign_source, auroc_clinvar_benign = "ClinVar benign",
                                auroc_gnomad_benign = "gnomAD AC\u22652 benign"))

p <- ggplot() +
  geom_vline(xintercept = 0.5, linewidth = 0.4, colour = "grey50", linetype = "dashed") +
  geom_vline(xintercept = 0.7, linewidth = 0.4, colour = oi("green"), linetype = "dotted") +
  geom_segment(data = g, aes(x = auroc_clinvar_benign, xend = auroc_gnomad_benign, y = gene, yend = gene),
               colour = "grey75", linewidth = 0.5) +
  geom_point(data = gl, aes(auroc, gene, colour = benign_source), size = 2.1) +
  scale_colour_manual(values = c("ClinVar benign" = oi("blue"), "gnomAD AC\u22652 benign" = oi("orange")), name = NULL) +
  scale_x_continuous(limits = c(0.5, 1.0)) +
  labs(x = "per-gene AUROC (ESM-2 650M, fp32 consistent; pathogenic vs benign)", y = NULL,
       title = "ClinVar transfer is robust to the benign source (gnomAD gives a similar-to-harder picture)") +
  theme_a3() + theme(panel.grid.major.y = element_blank(),
                     legend.position = c(0.98, 0.03), legend.justification = c(1, 0))
a3_save(p, "FS2_gnomad_control", 6.2, 5.0, FIG)
cat("FS2 done\n")
