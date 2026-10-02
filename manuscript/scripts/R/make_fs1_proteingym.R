#!/usr/bin/env Rscript
# A3 Fig S1 — ProteinGym variant-effect reliability atlas (B05): per-assay zero-shot
# discrimination AUROC across the DMS panel, with the pooled estimate. ggplot2 + Arial.
suppressMessages({ library(ggplot2); library(dplyr) })
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) "."); if (length(HERE) == 0 || HERE == "") HERE <- "scripts/R"
source(file.path(HERE, "theme_a3.R"))
BASE <- normalizePath(file.path(HERE, "..", "..")); PD <- file.path(BASE, "figures", "plotdata"); FIG <- file.path(BASE, "figures")
rd <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)

a <- rd("FS1_proteingym_atlas.csv")
a$assay <- sub("_[0-9]{4}$", "", a$assay)                 # drop trailing year
a$assay <- sub("_[A-Za-z]+$", "", a$assay)                # drop author token
a$assay <- substr(a$assay, 1, 22)
a <- a %>% arrange(auroc) %>% mutate(assay = factor(assay, levels = unique(assay)),
                                     disc = ifelse(ci_excludes_chance == 1, "discriminates", "at chance"))
pool <- rd("FS1_proteingym_pooled.csv")

p <- ggplot(a, aes(auroc, assay, colour = disc)) +
  annotate("rect", xmin = pool$ci_lo[1], xmax = pool$ci_hi[1], ymin = -Inf, ymax = Inf,
           fill = oi("blue"), alpha = 0.10) +
  geom_vline(xintercept = 0.5, linewidth = 0.4, colour = "grey50", linetype = "dashed") +
  geom_vline(xintercept = 0.7, linewidth = 0.4, colour = oi("green"), linetype = "dotted") +
  geom_vline(xintercept = pool$value[1], linewidth = 0.5, colour = oi("blue")) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.0, linewidth = 0.5) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = c("discriminates" = oi("blue"), "at chance" = oi("verm")), name = NULL) +
  scale_x_continuous(limits = c(0.3, 1.0)) +
  labs(x = "per-assay AUROC (ESM-2 650M zero-shot, position-held-out)", y = NULL,
       title = sprintf("pooled AUROC %.3f [%.3f, %.3f], %d/%d assays discriminate",
                       pool$value[1], pool$ci_lo[1], pool$ci_hi[1], pool$n_discriminating[1], pool$n_assays[1])) +
  theme_a3() + theme(panel.grid.major.y = element_blank(),
                     legend.position = "top", legend.justification = "left",
                     legend.direction = "horizontal")
a3_save(p, "FS1_proteingym_atlas", 6.2, 5.0, FIG)
cat("FS1 done\n")
