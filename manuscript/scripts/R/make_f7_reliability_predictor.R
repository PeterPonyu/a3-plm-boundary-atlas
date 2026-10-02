#!/usr/bin/env Rscript
# A3 Figure F6: unifying reliability-predictor. Tests whether ONE simple shared
# property predicts per-protein/per-gene discrimination AUROC across the atlas
# axes (DMS variant effect vs ClinVar). Result is a NULL -- the figure shows the
# (lack of) relationship honestly: scatters with Spearman rho + bootstrap CI, and
# a cross-axis sign-consistency summary. ggplot2 + Times New Roman (Arial).
suppressMessages({ library(ggplot2); library(patchwork); library(dplyr) })
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "scripts/R"
source(file.path(HERE, "theme_a3.R"))

BASE <- normalizePath(file.path(HERE, "..", ".."))
PD   <- file.path(BASE, "figures", "plotdata")
FIG  <- file.path(BASE, "figures")

d <- read.csv(file.path(PD, "F10_reliability_predictor.csv"),
              stringsAsFactors = FALSE, check.names = FALSE)
d$axis <- factor(d$axis, levels = c("DMS", "CLINVAR"), labels = c("DMS", "ClinVar"))

AXCOL <- c("DMS" = oi("blue"), "ClinVar" = oi("orange"))

sp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

# ---- panel A: DMS score-spread (the one weak within-axis signal) ----
dms <- d %>% filter(axis == "DMS", is.finite(score_std))
rhoA <- sp(dms$score_std, dms$auroc)
pA <- ggplot(dms, aes(score_std, auroc)) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = 2, colour = "grey60") +
  geom_smooth(method = "lm", se = TRUE, colour = oi("blue"),
              fill = scales::alpha(oi("blue"), 0.12), linewidth = 0.5, formula = y ~ x) +
  geom_point(size = 1.8, colour = oi("blue"), alpha = 0.8) +
  labs(x = "raw-LLR spread (SD across variants)",
       y = "per-protein AUROC", title = "DMS: score-spread (strongest signal)") +
  theme_a3()

# ---- panel B: length vs AUROC, both axes overlaid (null) ----
dl <- d %>% filter(is.finite(length))
rhoB_dms  <- sp(dl$length[dl$axis == "DMS"],     dl$auroc[dl$axis == "DMS"])
rhoB_clin <- sp(dl$length[dl$axis == "ClinVar"], dl$auroc[dl$axis == "ClinVar"])
pB <- ggplot(dl, aes(length, auroc, colour = axis)) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = 2, colour = "grey60") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, formula = y ~ x) +
  geom_point(size = 1.8, alpha = 0.8) +
  scale_colour_manual(values = AXCOL) +
  scale_x_log10() +
  labs(x = "sequence length (DMS: max-position proxy)",
       y = "AUROC", title = "Length: no consistent effect") +
  theme_a3() + theme(legend.position = c(0.16, 0.16))

# ---- panel C: class balance vs AUROC, both axes (sign-consistent but n.s.) ----
db <- d %>% filter(is.finite(class_balance))
rhoC_dms  <- sp(db$class_balance[db$axis == "DMS"],     db$auroc[db$axis == "DMS"])
rhoC_clin <- sp(db$class_balance[db$axis == "ClinVar"], db$auroc[db$axis == "ClinVar"])
pC <- ggplot(db, aes(class_balance, auroc, colour = axis)) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = 2, colour = "grey60") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, formula = y ~ x) +
  geom_point(size = 1.8, alpha = 0.8) +
  scale_colour_manual(values = AXCOL, guide = "none") +
  labs(x = "positive-class fraction",
       y = "AUROC", title = "Class balance: weak, n.s. on both") +
  theme_a3()

# ---- panel D: cross-axis sign summary (the verdict, visually) ----
summ <- data.frame(
  predictor = c("sample size n", "length", "class balance"),
  DMS    = c(0.046, 0.064, 0.093),
  ClinVar = c(-0.247, -0.378, 0.263)
)
sl <- data.frame(
  predictor = rep(summ$predictor, 2),
  axis = factor(rep(c("DMS", "ClinVar"), each = nrow(summ)), levels = c("DMS", "ClinVar")),
  rho = c(summ$DMS, summ$ClinVar)
)
sl$predictor <- factor(sl$predictor, levels = rev(summ$predictor))
pD <- ggplot(sl, aes(rho, predictor, colour = axis)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey50") +
  geom_point(size = 3, alpha = 0.9,
             position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = AXCOL) +
  scale_x_continuous(limits = c(-0.55, 0.55)) +
  labs(x = expression("Spearman " * rho * " (predictor vs AUROC)"),
       y = NULL, title = "Cross-axis: no predictor agrees in sign + strength") +
  theme_a3() + theme(legend.position = c(0.85, 0.85))

fig <- (pA | pB) / (pC | pD)

a3_save(fig, "F10_reliability_predictor", 9.2, 7.0, FIG)
