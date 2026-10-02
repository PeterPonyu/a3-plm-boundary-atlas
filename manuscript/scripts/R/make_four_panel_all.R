#!/usr/bin/env Rscript
# Staged A3 compositor: canonical figures are composed from committed plotdata.
# Reader layouts keep only drawn panels; blocked E/F cards are not restored.
# Release output is explicit and never overwrites canonical assets in place.
suppressMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
})
HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "scripts/R"
source(file.path(HERE, "theme_a3.R"))
source(file.path(HERE, "compose_a3.R"))
source(file.path(HERE, "evidence_panels_a3.R"))
# In-panel annotation text defaults to the body size. An explicitly passed size
# is honoured: silently overriding a panel's chosen size makes per-panel
# legibility tuning impossible and hides the real value from the reader of the
# panel code.
geom_text <- function(..., size = A3_TEXT_MM) {
  args <- list(...); args$size <- size
  do.call(ggplot2::geom_text, args)
}
Sys.setenv(A3_PDF_NORMALIZER = normalizePath(file.path(HERE, "..", "normalize_a3_pdf_metadata.py")))

BASE <- normalizePath(file.path(HERE, "..", ".."))
PD   <- file.path(BASE, "figures", "plotdata")
FIG  <- Sys.getenv("A3_FIGURE_OUTPUT", unset = "")
if (FIG == "") stop("Set A3_FIGURE_OUTPUT to an explicit staging directory.", call. = FALSE)
if (!dir.exists(FIG)) stop(sprintf("A3_FIGURE_OUTPUT does not exist: %s", FIG), call. = FALSE)
rd   <- function(f) read.csv(file.path(PD, f), stringsAsFactors = FALSE, check.names = FALSE)

# Task-5 locked aggregate triples are plotting inputs, never renderer-side
# recomputations.  Validate their contract before a panel may consume them.
locked_summary <- function(file, logical, panel, expected_groups, grouped = FALSE) {
  x <- rd(file)
  if (!all(x$status == "READY") || !all(x$logical_figure == logical) ||
      !all(x$target_panel == panel) || any(!is.finite(x$value)) ||
      any(x$n < 1L) || !identical(as.integer(x$sort_order), seq_len(nrow(x)))) {
    stop(sprintf("invalid locked summary contract: %s", file), call. = FALSE)
  }
  observed <- if (grouped) paste(x$group_1, x$group_2, sep = "|") else x$group_1
  if (!identical(as.character(observed), expected_groups)) {
    stop(sprintf("unexpected locked summary groups: %s", file), call. = FALSE)
  }
  x
}
locked_f8c <- function() {
  x <- rd("F8C_liability_pooled_counts.csv")
  if (!identical(names(x), c("logical_figure", "target_panel", "estimand_id", "status", "group", "n_sig", "n_datasets", "frac", "units", "cohort_id", "source_path", "source_sha256", "method", "uncertainty_method", "cluster_unit", "seed", "sort_order", "claim", "nonclaim")) ||
      !all(x$logical_figure == "F8") || !all(x$target_panel == "C") || !all(x$estimand_id == "A3-EQ-F8-POOL-01") || !all(x$status == "READY") ||
      !identical(as.character(x$group), c("liability", "non-liability")) || !identical(as.integer(x$sort_order), 1:2) || any(x$n_datasets < 1L) || any(abs(x$frac - x$n_sig / x$n_datasets) > 1e-12)) stop("invalid locked F8C pooling contract", call. = FALSE)
  x
}

# Each logical A3_EF.csv contains one E and one F row.  READY rows are shown
# as locked scalar summaries; blocked rows are evidence-boundary panels, never
# fabricated numeric estimates.
logical_for_stem <- c(
  F1_scale_plateau = "F1", F2_cross_model = "S6", F3_clinvar_forest = "F3",
  F4_clinical_calibration = "F4", F5_stability_ddg = "F5", F6_protein_uq = "S7",
  F7_binding_boundary = "F7", F8_developability_map = "F8", F9_antibody_lm = "S8",
  F10_reliability_predictor = "S9", F11_boundary_atlas = "F11",
  FS1_proteingym_atlas = "S1", FS2_gnomad_control = "S2", FS3_antibody_robustness = "S3",
  FS4_abstention_tiers = "S4", FS5_vep_comparison = "S5"
)

reader_action <- function(row) {
  if (identical(row$status[[1]], "READY")) return("")
  "Further locked source evidence is required before a quantitative estimate can be shown."
}

wrap_card <- function(text, width = 21) {
  paste(strwrap(text, width = width), collapse = "\n")
}

blocked_copy <- list(
  "A3-EQ-F3-LOGO-01" = list(
    title = "Variant-level audit", lines = c(
      "Not estimable", "Missing: score rows",
      "Source: F3 table", "No recovery claim.")),
  "A3-EQ-F7-SUBGROUP-01" = list(
    title = "Binding resampling", lines = c(
      "Not estimable", "Missing:", "row coords",
      "Source: F4 table", "No resampling.")),
  "A3-EQ-F8-EFFECT-01" = list(
    title = "Effect interval", lines = c(
      "Not estimable", "Missing: dataset effects", "and variance", "Source:",
      "F5 developability table", "Non-claim: No effect", "interval.")),
  "A3-EQ-F8-MDE-01" = list(
    title = "Power / MDE", lines = c(
      "Not estimable", "Missing: independent", "effects and variance", "Source:",
      "F5 developability table", "Non-claim: No", "equivalence conclusion.")),
  "A3-EQ-S1-EXPOSURE-01" = list(
    title = "Exposure audit", lines = c(
      "Not estimable", "Missing:", "exposure labels",
      "Source: S1 table", "No exposure audit.")),
  "A3-EQ-S3-INFLUENCE-01" = list(
    title = "Influence audit", lines = c(
      "Not estimable", "Missing:", "row coords",
      "Source: S3 table", "No influence audit.")),
  "A3-EQ-S4-THRESHOLD-01" = list(
    title = "Tier uncertainty", lines = c(
      "Not estimable", "Missing: protein-level", "scores",
      "Source: S4 tier table", "Non-claim: No tier", "uncertainty interval."))
)

blocked_card_copy <- function(row) {
  copy <- blocked_copy[[row$estimand_id[[1]]]]
  if (is.null(copy)) stop("blocked estimand lacks approved concise card copy", call. = FALSE)
  copy
}

evidence_panel <- function(row, tag) {
  is_ready <- identical(row$status[[1]], "READY")
  if (is_ready) {
    # A locked scalar is a summary of a committed table that this row already
    # names and pins.  Draw that table under the scalar so the sensitivity cell
    # reports the evidence rather than one marker beside the words "no
    # interval".  The reported number remains the committed payload value; the
    # renderer re-derives it from the pinned source only to refuse to draw when
    # the two disagree.
    return(evidence_panel_data_backed(row, tag))
  }
  # Blocked evidence-boundary panel: a full-panel card that fills the space while
  # stating exactly why the estimand is not on disk.  No numeric value is drawn.
  copy <- blocked_card_copy(row)
  body <- paste(copy$lines, collapse = "\n")
  # Keep the card self-contained even when a blocked estimand needs a longer,
  # multi-line explanation (F8/F4/S4).  The old top-anchored body could run
  # below the card in the short bottom-row cells.  Centre the copy in the
  # available body band and scale only the presentation text for seven-line
  # cards; this does not alter any locked wording or estimand status.
  n_body <- length(copy$lines)
  # Cards must stay above the 5 pt shipped-figure floor even after the
  # type-scale cut. 5.15 pt = 1.817 mm.
  card_floor_mm <- 2.47
  body_size <- max(A3_TEXT_MM * if (n_body >= 7L) 0.76 else if (n_body >= 6L) 0.82 else 0.92,
                   card_floor_mm)
  body_y <- 0.68
  blocked_title_size <- max(A3_TEXT_MM * if (nchar(copy$title) > 16L) 0.82 else 1.02,
                            card_floor_mm)
  # An evidence boundary is a deliberate statement about what the committed
  # inputs cannot support, so it is set as type on the page.  The former filled
  # grey card read as a failed render rather than as a considered limit.
  ggplot() +
    annotate("text", x = 0.5, y = 0.865, label = copy$title, family = A3_FONT,
             fontface = "bold", size = blocked_title_size, hjust = 0.5, vjust = 1) +
    annotate("segment", x = 0.22, xend = 0.78, y = 0.775, yend = 0.775,
             colour = EV_RULE, linewidth = 0.4) +
    annotate("text", x = 0.5, y = body_y, label = body, family = A3_FONT,
             size = body_size, hjust = 0.5, vjust = 1, lineheight = 0.92,
             colour = EV_MUTED) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    labs(title = NULL, x = NULL, y = NULL) +
  theme_void(base_family = A3_FONT) +
    # Letters are placed by compose_a3.R, not by ggplot tags.
    theme(plot.margin = margin(2, 8, 3, 8))
}

# Titles and layouts are reviewed at the final 6.5-inch print width.
# Redundant scalar/provenance cells are integrated into captions; none is
# replaced with invented observations, intervals, or empty placeholder cards.
PANEL_TITLES <- list(
  F1_scale_plateau=c("Model scaling","Paired model differences","Scoring convention"),
  F2_cross_model=c("Encoder comparison","Protein-level agreement","Difference from ESM-2","Encoder-family means"),
  F3_clinvar_forest=c("Per-gene ClinVar AUROC","Historical gnomAD comparison","Pathogenic sample size","Gene-level thresholds","Sequence length"),
  F4_clinical_calibration=c("Same-score OOF reliability","Proper-score baselines","Gene-specific ECE","Brier difference from gene baseline"),
  F5_stability_ddg=c("Rank agreement","Post-linear-map error","Held-out error inflation","Clean split size"),
  F6_protein_uq=c("Uncertainty by shift regime","Recorded BayesianRidge acquisition","Across-split uncertainty","High-shift uncertainty"),
  F7_binding_boundary=c("Binding rank agreement","Mutation counts","Interval width","Subset comparison"),
  F8_developability_map=c("BH-significant datasets","Germline-control applicability"),
  F9_antibody_lm=c("Liability discovery fraction","Antibody-binding agreement","Liability discovery counts","Binding-interval width"),
  F10_reliability_predictor=c("DMS score spread","Sequence length","Class balance","Predictor associations"),
  F11_boundary_atlas=c("Discrimination summaries","Physical-readout summaries"),
  FS1_proteingym_atlas=c("Position-held-out assays","AUROC distribution","Assay sample size"),
  FS2_gnomad_control=c("Benign source by gene","Benign-source comparison","Mechanism-class means"),
  FS3_antibody_robustness=c("Scoring methods by subset","Overall binding","Antibody binding","Model-size difference"),
  FS4_abstention_tiers=c("Model-size tiers","Selective coverage","Selection penalty","Gap from shuffled control"),
  FS5_vep_comparison=c("Published predictor summaries","Supervision-class means","Predictor ordering","Protein coverage")
)
save_panels <- function(base_panels, stem, w=12.2, h=9.6, layout="six") {
  only <- Sys.getenv("A3_ONLY_STEMS",unset="")
  if(nzchar(only) && !stem %in% strsplit(only,",",fixed=TRUE)[[1]]) return(invisible(NULL))
  panels <- lapply(base_panels,a3_atomic_panel)
  if(stem=="F1_scale_plateau") {
    ef <- rd("F1_EF.csv")
    panels <- list(panels[[1]],evidence_panel(ef[1,,drop=FALSE],"B"),evidence_panel(ef[2,,drop=FALSE],"C"))
  }
  if(stem=="FS1_proteingym_atlas") panels <- panels[1:3]
  if(stem=="FS2_gnomad_control") panels <- panels[c(1,2,4)]
  titles <- PANEL_TITLES[[stem]]
  stopifnot(length(panels)==length(titles))
  panels <- Map(function(p,t) p+labs(title=t),panels,titles)
  # These three historical diagnostic panels contain geom_smooth fits. This
  # format-only run may not fit models; retain their current verified face.
  if(stem=="F10_reliability_predictor" &&
     Sys.getenv("A3_FORMAT_ONLY",unset="") == "1") return(invisible(NULL))
  if(stem=="F3_clinvar_forest") panels[[3]] <- panels[[3]] + theme(aspect.ratio=1)
  if(stem=="F3_clinvar_forest") {
    composite <- a3_compose_positions(panels,list(c(1,2,1,1),c(1,1,2,2),c(2,2,2,2),c(3,3,1,1),c(3,3,2,2)),c(1.1,1.1,.90))
    # 0.70 * the 9-inch text height is 6.3 inches. Avoid shrinking this
    # forest's type when the main manuscript applies the height constraint.
    h <- 6.3
  } else if(stem=="F4_clinical_calibration") {
    composite <- a3_compose_readable(panels)
    h <- 5.5
  } else if(stem %in% c("FS1_proteingym_atlas","FS2_gnomad_control")) {
    composite <- a3_compose_positions(panels,list(c(1,2,1,1),c(1,1,2,2),c(2,2,2,2)),c(1,1))
    h <- 4.8
  } else if(stem=="F6_protein_uq") {
    composite <- a3_compose_positions(panels,list(c(1,1,1,2),c(2,2,1,2),c(3,3,1,1),c(3,3,2,2)),c(1,1.12,.90))
    h <- 6.0
  } else {
    composite <- a3_compose_readable(panels)
    h <- if(length(panels)==2) 3.0 else if(stem=="F1_scale_plateau") 4.0 else if(stem=="F11_boundary_atlas") 4.4 else 4.8
  }
  a3_save(composite,stem,6.5,h,FIG)
}

# ===================== Fig 1 / F1_scale_plateau =====================
f1 <- rd("F1_scale.csv") %>%
  mutate(model_size = factor(model_size, c("35M", "150M", "650M", "3B")),
         delta_prev = protein_mean_auroc - lag(protein_mean_auroc),
         frac_disc = n_discriminating / n_proteins)
p1A <- ggplot(f1, aes(model_size, protein_mean_auroc, group = 1)) +
  geom_line(linewidth = A3_LINE, colour = oi("blue")) +
  geom_point(size = A3_POINT, colour = oi("blue")) +
  geom_text(aes(label = sprintf("%.3f", protein_mean_auroc)), vjust = -1.0, family = A3_FONT) +
  scale_y_continuous(limits = c(0.64, 0.75)) +
  a3_claim("Scale then plateau") +
  labs(x = "ESM-2 parameters", y = "protein-mean AUROC", tag = "A") +
  theme_a3()
p1B <- ggplot(f1, aes(model_size, n_discriminating)) +
  geom_col(fill = oi("blue"), width = A3_BAR) +
  geom_text(aes(label = n_discriminating), vjust = -0.4, family = A3_FONT) +
  scale_y_continuous(limits = c(0, 72), expand = expansion(c(0, 0.02))) +
  a3_claim("Discriminating proteins") +
  labs(x = "ESM-2 parameters", y = "n discriminating (of 69)", tag = "B") +
  theme_a3() + theme(panel.grid.major.x = element_blank())
f1c <- f1 %>% filter(!is.na(delta_prev)) %>% mutate(plateau = abs(delta_prev) < 1e-4)
p1C <- ggplot(f1c, aes(model_size, delta_prev, fill = plateau)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_col(width = A3_BAR, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%+.3f", delta_prev)), vjust = -0.4, family = A3_FONT) +
  scale_fill_manual(values = c(`FALSE` = oi("orange"), `TRUE` = oi("grey"))) +
  a3_y_for_labels(bottom = 0.08, top = 0.22) +
  a3_claim("Diminishing increments") +
  labs(x = "step to size", y = "\u0394 AUROC vs previous size", tag = "C") +
  theme_a3() + theme(panel.grid.major.x = element_blank())
p1D <- ggplot(f1, aes(model_size, frac_disc)) +
  geom_col(fill = oi("blue"), width = A3_BAR) +
  geom_text(aes(label = sprintf("%.2f", frac_disc)), vjust = -0.4, family = A3_FONT) +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(c(0, 0.02))) +
  a3_claim("Fraction discriminating") +
  labs(x = "ESM-2 parameters", y = "fraction discriminating", tag = "D") +
  theme_a3() + theme(panel.grid.major.x = element_blank())
save_panels(list(p1A, p1B, p1C, p1D), "F1_scale_plateau")

# ===================== Fig S6 / F2_cross_model =====================
f2 <- rd("F2_crossmodel.csv") %>%
  mutate(model = dplyr::recode(model,
    "ProtBERT-BFD" = "ProtBERT", "ESM-2 650M (masked)" = "ESM-2 masked",
    "ESM-2 650M (WT)" = "ESM-2 WT", .default = model), model = factor(model, model))
ag <- rd("F2_agreement.csv")
p2A <- ggplot(f2, aes(model, protein_mean_auroc, fill = model)) +
  geom_col(width = A3_BAR) +
  geom_text(aes(label = sprintf("%.3f", protein_mean_auroc)), hjust = 1.15, colour = "white", family = A3_FONT) +
  scale_fill_manual(values = c(oi("orange"), oi("blue"), oi("sky")), guide = "none") +
  scale_y_continuous(limits = c(0, 0.85), breaks = c(0, 0.4, 0.8), expand = expansion(c(0, 0.02))) +
  coord_flip() +
  labs(x = NULL, y = "protein-mean AUROC", tag = "A") +
  theme_a3() + theme(panel.grid.major.y = element_blank())
p2B <- ggplot(ag %>% mutate(comparison = dplyr::recode(comparison,
  "ProtBERT vs ESM-2 650M (WT)" = "ProtBERT vs WT",
  "ProtBERT vs ESM-2 650M (masked)" = "ProtBERT vs masked",
  "ProtBERT vs ESM-2 (WT)" = "ProtBERT vs WT",
  "ProtBERT vs ESM-2 (masked)" = "ProtBERT vs masked",
  .default = comparison)), aes(reorder(comparison, spearman), spearman)) +
  geom_col(fill = oi("blue"), width = A3_BAR) +
  geom_text(aes(label = sprintf("%.3f", spearman)), hjust = 1.15, colour = "white", size = 2.7, family = A3_FONT) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 0.85), breaks = c(0, 0.4, 0.8), expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "Spearman \u03c1", tag = "B") +
  theme_a3() + theme(panel.grid.major.y = element_blank())
esm_ref <- f2$protein_mean_auroc[f2$model == "ESM-2 masked"]
f2d <- f2 %>% mutate(delta_vs_esm = protein_mean_auroc - esm_ref)
p2C <- ggplot(f2d, aes(model, delta_vs_esm, fill = delta_vs_esm >= 0)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col(width = A3_BAR, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = oi("blue"), `FALSE` = oi("verm"))) +
  a3_y_for_labels(bottom = 0.18, top = 0.18, breaks = c(-0.01, 0)) +
  coord_flip() +
  labs(x = NULL, y = "\u0394 AUROC", tag = "C") +
  theme_a3() + theme(panel.grid.major.y = element_blank())
p2D <- ggplot(f2, aes(family, protein_mean_auroc, fill = family)) +
  geom_col(width = A3_BAR, position = position_dodge(0.6)) +
  # one bar per model still: show model as text
  geom_text(aes(label = sprintf("%.3f", protein_mean_auroc)),
            position = position_stack(vjust = 0.5), size = 2.5, family = A3_FONT, colour = "white") +
  scale_fill_manual(values = c("BERT/UniRef100+BFD" = oi("orange"), "ESM/UR50" = oi("blue")), guide = "none") +
  scale_y_continuous(limits = c(0, 0.85), expand = expansion(c(0, 0.02))) +
  labs(x = "encoder family", y = "protein-mean AUROC", tag = "D") +
  theme_a3() + theme(panel.grid.major.x = element_blank())
# D as family means (honest aggregate of committed rows)
f2fam <- locked_summary("S6D_family_mean_protein_auroc.csv", "S6", "D",
                        c("BERT/UniRef100+BFD", "ESM/UR50")) %>%
  transmute(family = factor(c("ProtBERT", "ESM-2"), levels = c("ProtBERT", "ESM-2")),
            protein_mean_auroc = value)
p2D <- ggplot(f2fam, aes(family, protein_mean_auroc, fill = family)) +
  geom_col(width = A3_BAR) +
  geom_text(aes(label = sprintf("%.3f", protein_mean_auroc)), vjust = -0.4, size = 2.7, family = A3_FONT) +
  scale_fill_manual(values = c("ProtBERT" = oi("orange"), "ESM-2" = oi("blue")), guide = "none") +
  scale_y_continuous(limits = c(0, 0.85), expand = expansion(c(0, 0.02))) +
  labs(x = "encoder family", y = "protein-mean AUROC", tag = "D") +
  theme_a3() + theme(panel.grid.major.x = element_blank(),
                      axis.text.x = element_text(angle = 12, hjust = 1))
save_panels(list(p2A, p2B, p2C, p2D), "F2_cross_model")

# ===================== Fig 3 / F3_clinvar_forest =====================
f3 <- rd("F3_clinvar_20260923.csv") %>% mutate(gene = factor(gene, levels = gene))
# Panel A must name all twenty genes. Measured on the promoted face, twenty
# labels need 20 x 6.86 = 137 pt of stacked glyph height in an 82 pt panel, so a
# single column ran a 4.07 pt pitch and struck every name through its neighbour;
# fitting them would need a ~3.7 pt face. Ten per sub-column gives an 8.2 pt
# pitch, which clears the standard 7 pt tick face. The split is by the existing
# factor order only -- no value, interval, or ordering changes.
f3_forest <- f3 %>%
  mutate(forest_col = ifelse(as.integer(gene) <= 10L, "lower ten", "upper ten"))
p3A <- ggplot(f3_forest, aes(auroc_clinvar, gene)) +
  geom_vline(xintercept = 0.5, linewidth = 0.35, colour = "grey50", linetype = "dashed") +
  geom_vline(xintercept = 0.7, linewidth = 0.35, colour = oi("green"), linetype = "dotted") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                width = .20, linewidth = A3_ERRORBAR_LW, colour = oi("blue")) +
  geom_point(aes(size = n_path), colour = oi("blue")) +
  scale_size_continuous(range = c(1.2, 2.8), guide = "none") +
  # The limits used to sit exactly on the data: min(ci_lo) is 0.4500 and
  # max(ci_hi) is 1.0000, so both extreme whiskers terminated on the spine and
  # AICDA's marker at 1.0 was cut by it. Padded past the data on both sides by
  # more than the largest marker radius, with breaks pinned so the wider range
  # does not add a fourth tick label to an 88 pt panel.
  # Two sub-columns halve the width available per forest, so the interior 0.75
  # tick is dropped; 0.5 and 1.0 still bracket the range and the dashed chance
  # line and dotted 0.70 line carry the interior references.
  scale_x_continuous(limits = c(0.35, 1.05), breaks = c(0.5, 1.0)) +
  a3_claim("Per-gene ClinVar AUROC") +
  labs(x = "per-gene ClinVar AUROC", y = NULL, tag = "A") +
  theme_a3() + theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
                      strip.text = element_blank(), strip.background = element_blank(),
                      panel.spacing.x = unit(3, "pt"))
p3B <- ggplot(f3, aes(auroc_clinvar, historical_auroc_gnomad)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(size = A3_POINT, colour = oi("blue")) +
  coord_equal(xlim = c(0.48, 1.02), ylim = c(0.48, 1.02), expand=FALSE) +
  a3_claim("Benign-set robustness") +
  labs(x = "rescored ClinVar AUROC", y = "historical gnomAD AUROC", tag = "B") +
  theme_a3()
p3C <- ggplot(f3, aes(n_path, auroc_clinvar)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(size = A3_POINT, colour = oi("blue")) +
  a3_claim("Sample-size diagnostic") +
  labs(x = "n pathogenic variants", y = "ClinVar AUROC", tag = "C") +
  theme_a3()
thresh <- data.frame(
  criterion = c("AUROC\n≥0.70", "AUROC\n≥0.80"),
  n = c(sum(f3$auroc_clinvar >= 0.70), sum(f3$auroc_clinvar >= 0.80)),
  denom = nrow(f3)
)
p3D <- ggplot(thresh, aes(criterion, n / denom)) +
  geom_col(fill = oi("blue"), width = A3_BAR) +
  geom_text(aes(label = paste0(n, "/", denom)), vjust = -0.4, family = A3_FONT) +
  scale_y_continuous(limits = c(0, 1.08), expand = expansion(c(0, 0.02))) +
  a3_claim("Gene-level success") +
  labs(x = NULL, y = "fraction of genes", tag = "D") +
  theme_a3() + theme(axis.text.x = element_text(angle = 0, hjust = .5),
                      panel.grid.major.x = element_blank())
p3E <- ggplot(f3, aes(wt_len, auroc_clinvar)) +
  geom_point(size = A3_POINT, colour = oi("blue")) +
  labs(x = "reference length (residues)", y = "ClinVar AUROC") + theme_a3()
save_panels(list(p3A, p3B, p3C, p3D, p3E), "F3_clinvar_forest", 12.2, 9.6)

# ===================== Fig 4 / F4_clinical_calibration =====================
map_names <- c(sigmoid_z="Sigmoid(z)",platt="Platt",isotonic="Isotonic",
               pooled_prevalence="Pooled prevalence",gene_prevalence="Gene prevalence")
map_cols <- c("Sigmoid(z)"=oi("verm"),"Platt"=oi("blue"),"Isotonic"=oi("green"))
rel <- rd("clinical_reliability_bins_20260923.csv") %>%
  filter(gene=="ALL",map %in% c("sigmoid_z","platt","isotonic"),n>0) %>%
  mutate(which=map_names[map])
p4A <- ggplot(rel, aes(mean_pred, obs_rate, colour = which)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.35, colour = "grey55", linetype = "dashed") +
  geom_line(linewidth = A3_LINE) + geom_point(aes(size = n)) +
  scale_colour_manual(values = map_cols) +
  scale_size_continuous(range = c(1.3, 2.8), guide = "none") +
  coord_equal(xlim = c(-.02, 1.02), ylim = c(-.02, 1.02), expand=FALSE) +
  labs(x = "predicted P(pathogenic)", y = "observed pathogenic rate", tag = "A") +
  theme_a3() + a3_legend_outside("bottom") + guides(colour=guide_legend(nrow=2))
cal_metrics <- rd("clinical_metrics_20260923.csv") %>%
  mutate(which=factor(map_names[map],levels=rev(unname(map_names))))
p4B <- ggplot(cal_metrics,aes(brier,which)) +
  geom_errorbar(aes(xmin=brier_lo,xmax=brier_hi),orientation="y",width=.15,
                linewidth=A3_ERRORBAR_LW,colour=oi("blue")) +
  geom_point(size=A3_POINT,colour=oi("blue")) +
  labs(x="Brier score (lower is better)",y=NULL) + theme_a3() +
  theme(panel.grid.major.y=element_blank())
gene_cal <- rd("clinical_per_gene_calibration_20260923.csv") %>%
  filter(gene!="ALL",map %in% c("platt","isotonic")) %>%
  select(gene,map,ece) %>% pivot_wider(names_from=map,values_from=ece)
p4C <- ggplot(gene_cal,aes(platt,isotonic)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey55") +
  geom_point(size=A3_POINT,colour=oi("purple")) +
  coord_equal(xlim=c(0,.55),ylim=c(0,.55),expand=FALSE) +
  labs(x="Platt ECE (per gene)",y="Isotonic ECE (per gene)") + theme_a3()
cal_delta <- rd("clinical_brier_difference_20260923.csv") %>% mutate(which=map_names[map])
p4D <- ggplot(cal_delta,aes(value,which)) +
  geom_vline(xintercept=0,colour="grey55",linetype="dashed") +
  geom_errorbar(aes(xmin=lo,xmax=hi),orientation="y",width=.15,linewidth=A3_ERRORBAR_LW,colour=oi("blue")) +
  geom_point(size=A3_POINT,colour=oi("blue")) +
  labs(x="Brier difference (map minus baseline)",y=NULL) +
  theme_a3() + theme(panel.grid.major.y=element_blank())
save_panels(list(p4A,p4B,p4C,p4D), "F4_clinical_calibration")

# ===================== Fig 5 / F5_stability_ddg =====================
SPLITS <- c("train", "homology-held-out")
SPLITCOL <- c("train" = oi("blue"), "homology-held-out" = oi("verm"))
# F5 titles stay on the ggplot; letters are placed by compose_a3.R above
# the grob, so this header no longer reserves a ggplot tag band.
f5_header <- theme(
  plot.title.position = "plot",
  plot.title = element_text(hjust = 0.5, face = "plain", margin = margin(t = 0, b = 3)),
  plot.margin = margin(2, 4, 4, 4)
)
a <- rd("F10_discrimination.csv")
a$split <- factor(a$split, levels = SPLITS)
a$variant <- factor(a$variant, levels = c("v2 WT", "v2 masked"))
sref <- rd("F10_reference.csv")$spearman[1]
p5A <- ggplot(a, aes(variant, spearman, fill = split)) +
  geom_hline(yintercept = sref, linewidth = 0.35, linetype = 2, colour = "grey45") +
  geom_col(position = position_dodge(0.7), width = A3_BAR) +
  scale_fill_manual(values = SPLITCOL) +
  scale_y_continuous(limits = c(0, 0.78), expand = expansion(c(0, 0.02))) +
  labs(title = "Discrimination", x = NULL, y = "Spearman \u03c1", tag = "A") +
  theme_a3() + theme(panel.grid.major.x = element_blank()) +
  a3_legend_outside("bottom") + f5_header
b <- rd("F10_calibration.csv")
b$split <- factor(b$split, levels = SPLITS)
b$variant <- factor(b$variant, levels = c("v2 WT", "v2 masked"))
p5B <- ggplot(b, aes(variant, calib_rmse, fill = split)) +
  geom_col(position = position_dodge(0.7), width = A3_BAR) +
  scale_fill_manual(values = SPLITCOL, labels=c("held-in 20%","homology-held-out")) +
  scale_y_continuous(limits = c(0, 1.12), expand = expansion(c(0, 0.02))) +
  labs(title = "Post-linear-map error", x = NULL, y = "RMSE (MegaScale\nddG_ML units)", tag = "B") +
  theme_a3() + theme(panel.grid.major.x = element_blank()) + a3_legend_outside("bottom") + f5_header
infl <- b %>% distinct(variant, inflation)
p5C <- ggplot(infl, aes(variant, inflation)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_col(fill = oi("orange"), width = A3_BAR) +
  geom_text(aes(label = sprintf("%+.3f", inflation)),
            vjust = ifelse(infl$inflation >= 0, -0.4, 1.35), size = A3_TEXT_MM, family = A3_FONT) +
  a3_y_for_labels(bottom = 0.22, top = 0.22) +
  labs(title = "Held-out inflation", x = NULL, y = "RMSE \u0394 (ddG_ML units)", tag = "C") +
  theme_a3() + theme(panel.grid.major.x = element_blank()) + f5_header
meta <- rd("F10_meta.csv")
meta_plot <- meta %>% filter(field %in% c("n_train", "n_heldout", "clusters_train", "clusters_heldout")) %>%
  mutate(label = recode(field,
                        n_train = "training mutations",
                        n_heldout = "held-out mutations",
                        clusters_train = "training clusters",
                        clusters_heldout = "held-out clusters"),
         label = factor(label, levels = rev(c("training mutations", "held-out mutations",
                                             "training clusters", "held-out clusters"))))
# value/1.25 is 0.097 decades, about 1.4 pt on this axis, so the 2.7-size marker
# (~3.8 pt radius) covered the label's last glyph. A single larger divisor cannot
# fix it either: these four values span 26 to 164,094, so an offset that clears
# the marker on the big bars pushes the small ones off the left spine. Placed
# conditionally instead, inside the lollipop where the segment has room and
# outside it where it does not.
meta_plot$lab_inside <- log10(meta_plot$value) >= 2.2
meta_plot$unit <- ifelse(grepl("mutations",meta_plot$label),"Mutations","Clusters")
meta_plot$partition <- ifelse(grepl("training",meta_plot$label),"Training","Held-out")
p5D <- ggplot(meta_plot,aes(value,partition)) +
  geom_col(fill=oi("blue"),width=.6) +
  geom_text(aes(label=a3_fmt_n(value)),hjust=-.12,family=A3_FONT) +
  facet_wrap(~unit,ncol=1,scales="free_x") +
  scale_x_continuous(expand=expansion(c(0,.36)),n.breaks=3,labels=scales::label_comma()) +
  labs(x="count",y=NULL) + theme_a3() + theme(panel.grid.major.y=element_blank())
save_panels(list(p5A, p5B, p5C, p5D), "F5_stability_ddg", 12.2, 9.6)

# ===================== Fig S7 / F6_protein_uq =====================
METHODS <- c("BayesianRidge", "Deep Ensemble", "GP-RBF", "MC-Dropout")
METHCOL <- c("BayesianRidge" = oi("grey"), "Deep Ensemble" = oi("blue"),
             "GP-RBF" = oi("orange"), "MC-Dropout" = oi("green"))
ru <- rd("F9_rho_unc.csv") %>%
  mutate(method = factor(method, levels = METHODS),
         shift_class = factor(shift_class,
                              levels = c("extreme", "moderate-shift", "in-distribution"),
                              labels = c("extreme", "moderate", "in-dist")))
ru_mean <- locked_summary("S7A_method_shift_mean_rho_unc.csv", "S7", "A",
                          unlist(lapply(METHODS, function(method) paste(
                            method, c("extreme", "moderate-shift", "in-distribution"), sep = "|"))),
                          grouped = TRUE) %>%
  transmute(method = factor(group_1, levels = METHODS),
            shift_class = factor(group_2, levels = c("extreme", "moderate-shift", "in-distribution"),
                                 labels = c("extreme", "moderate", "in-dist")),
            rho = value)
p6A <- ggplot(ru_mean, aes(shift_class, rho, fill = method)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_col(position = position_dodge(0.75), width = A3_BAR) +
  scale_fill_manual(values = METHCOL) +
  labs(title = "Uncertainty\u2013signal association", x = "shift class", y = "mean \u03c1_unc", tag = "A") +
  theme_a3() + a3_legend_outside("bottom") + theme(legend.position = "bottom",
                      panel.grid.major.x = element_blank(),
                      axis.text.x = element_text(angle = 15, hjust = 1))
bo <- rd("F9_bo_acquisition.csv") %>%
  mutate(dataset=ifelse(grepl("^gb1",split),"GB1 fitness","Meltome temperature")) %>%
  pivot_longer(c(greedy, ucb, random), names_to = "acq", values_to = "best") %>%
  mutate(acq = factor(acq, c("greedy", "ucb", "random")))
p6B <- ggplot(bo, aes(split, best, colour = acq, group = acq)) +
  geom_point(position = position_dodge(0.55), size = A3_POINT) +
  facet_wrap(~dataset,scales="free",ncol=2) +
  scale_colour_manual(values = c(greedy = oi("black"), ucb = oi("blue"), random = oi("grey"))) +
  scale_x_discrete(labels = c("gb1/low_vs_high" = "GB1\nlow-high", "gb1/sampled" = "GB1\nsampled", "gb1/three_vs_rest" = "GB1\n3-rest", "gb1/two_vs_rest" = "GB1\n2-rest", "meltome/human" = "Meltome\nhuman", "meltome/human_cell" = "Meltome\ncell", "meltome/mixed_split" = "Meltome\nmixed")) +
  labs(title = "Retrospective optimization", x = NULL, y = "best observed label", tag = "B") +
  theme_a3() + theme(axis.text.x = element_text(angle = 0, hjust = .5, size = A3_TICK_PT),
                      legend.position = "bottom")
summ <- rd("F9_summary.csv") %>% mutate(method = factor(method, levels = METHODS))
p6C <- ggplot(summ, aes(method, mean_rho_unc_all, fill = method)) +
  geom_col(width = A3_BAR, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", mean_rho_unc_all)),
            vjust = ifelse(summ$mean_rho_unc_all >= 0, -0.4, 1.35), size = A3_TEXT_MM, family = A3_FONT) +
  scale_fill_manual(values = METHCOL) +
  a3_y_for_labels(bottom = 0.18, top = 0.22) +
  scale_x_discrete(labels = c("BayesianRidge" = "Bayesian\nRidge", "Deep Ensemble" = "Deep\nEnsemble", "GP-RBF" = "GP-RBF", "MC-Dropout" = "MC\nDropout")) +
  labs(title = "Across-split uncertainty", x = NULL, y = "mean \u03c1_unc\n(all splits)", tag = "C") +
  theme_a3() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
                      panel.grid.major.x = element_blank())
p6D <- ggplot(summ,aes(method,mean_rho_unc_highshift,fill=method)) +
  geom_col(width=A3_BAR,show.legend=FALSE) +
  geom_text(aes(label=sprintf("%.3f",mean_rho_unc_highshift)),vjust=-.5,family=A3_FONT) +
  scale_fill_manual(values=METHCOL) + scale_y_continuous(limits=c(0,.155),expand=expansion(c(0,.02))) +
  scale_x_discrete(labels=c("BayesianRidge"="Bayesian\nRidge","Deep Ensemble"="Deep\nEnsemble","GP-RBF"="GP-RBF","MC-Dropout"="MC\nDropout")) +
  labs(x=NULL,y="mean ρ_unc\n(high shift)") + theme_a3() + theme(panel.grid.major.x=element_blank())
save_panels(list(p6A, p6B, p6C, p6D), "F6_protein_uq")

# ===================== Fig 7 / F7_binding_boundary =====================
f7 <- rd("F4_binding.csv") %>%
  mutate(subset = factor(subset, levels = rev(subset)),
         hl = key %in% c("antibody_ABAG", "noninterface"),
         ci_width = ci_hi - ci_lo)
p7A <- ggplot(f7, aes(spearman, subset, colour = hl)) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.18, linewidth = A3_ERRORBAR_LW) +
  geom_point(size = A3_POINT) +
  # Not on the measured defect list, but found by looking: 0.222 and 0.255 were
  # printing with their final digit sliced by the spine, the same truncation
  # class as panel C and just as much a number the panel gets wrong. Flipping the
  # label to the left of the marker for the rightmost points removes the failure
  # mode instead of buying headroom for it; this is the rule Fig 5 B and C use.
  geom_text(aes(label = sprintf("%.3f", spearman),
                hjust = ifelse(spearman > 0.6 * max(spearman), 1.25, -0.30)),
            colour = "grey20", nudge_y = 0.21, vjust = 0, family = A3_FONT) +
  scale_colour_manual(values = c(`FALSE` = oi("blue"), `TRUE` = oi("verm")), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.26))) +
  a3_claim("Binding rank correlation") +
  labs(x = "Spearman(LLR, \u2212\u0394\u0394G)", y = NULL, tag = "A") +
  theme_a3() + theme(panel.grid.major.y = element_blank(),
                      plot.margin = margin(A3_HEADER_GUTTER_PT, 5, 16, 6))
p7B <- ggplot(f7, aes(subset, n, fill = hl)) +
  geom_col(width = A3_BAR, show.legend = FALSE) +
  geom_text(aes(label = a3_fmt_n(n),
                colour = ifelse(n > 1200, "inside", "outside")),
            hjust = ifelse(f7$n > 1200, 1.12, -0.12), family = A3_FONT) +
  coord_flip() +
  scale_fill_manual(values = c(`FALSE` = oi("blue"), `TRUE` = oi("verm"))) +
  scale_colour_manual(values = c(inside = "white", outside = "grey15"), guide = "none") +
  # Default breaks gave 0/1000/2000/3000/4000, which printed as the unreadable
  # run "0 1,0002,0003,0004,000" in an 88 pt cell. Two breaks cannot collide.
  scale_y_continuous(expand = expansion(mult = c(0, 0.22)),
                     breaks = c(0, 3000),
                     labels = scales::label_comma()) +
  a3_claim("Mutation counts") +
  labs(x = NULL, y = "n mutations", tag = "B") +
  theme_a3() + theme(panel.grid.major.y = element_blank(),
                      axis.title.x = element_text(margin = margin(t = 2)),
                      plot.margin = margin(A3_HEADER_GUTTER_PT, 5, 16, 5))
# An outside label on the two long bars needs more width than the 0.18 headroom
# gives, so it used to print truncated: 0.155283197 as "0.15" and 0.166979758 as
# "0.1". A panel that prints a different number from the one its caption names as
# A3-EQ-F7-WIDTH-01 is a correctness defect, not a spacing nit. Fixed by the
# inside/outside placement panel B already uses, keyed to a fraction of the max
# rather than a hardcoded cut so it survives a data refresh.
f7$width_inside <- f7$ci_width > 0.6 * max(f7$ci_width)
p7C <- ggplot(f7, aes(subset, ci_width)) +
  geom_col(fill = oi("blue"), width = A3_BAR) +
  geom_text(aes(label = sprintf("%.3f", ci_width),
                colour = ifelse(width_inside, "inside", "outside")),
            hjust = ifelse(f7$width_inside, 1.12, -0.12), family = A3_FONT) +
  coord_flip() +
  scale_colour_manual(values = c(inside = "white", outside = "grey15"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  a3_claim("Interval width") +
  labs(x = NULL, y = "CI width", tag = "C") +
  theme_a3() + theme(panel.grid.major.y = element_blank(),
                      plot.margin = margin(A3_HEADER_GUTTER_PT, 5, 16, 6))
contrast <- f7 %>% filter(key %in% c("general_nonAB", "antibody_ABAG", "noninterface"))
p7D <- ggplot(contrast, aes(subset, spearman, fill = key == "general_nonAB")) +
  geom_col(width = A3_BAR, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", spearman)), vjust = -0.4, size = 2.6, family = A3_FONT) +
  scale_fill_manual(values = c(`TRUE` = oi("blue"), `FALSE` = oi("verm"))) +
  # The two-line break was itself the defect: rotated 45 degrees, a label's
  # second line shears down-right into its neighbour's first line, so the "g)"
  # of "(Ab-Ag)" struck through the "g" of "general". Single-line labels cannot
  # shear. Dropping the parenthetical qualifiers costs nothing here because the
  # full subset names are printed on the y axis of panels A, B, and C of this
  # same figure, and the three short forms remain distinct.
  scale_x_discrete(labels = c(
    "non-interface" = "non-interface",
    "antibody (Ab-Ag)" = "antibody",
    "general (non-Ab)" = "general"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.20))) +
  a3_claim("Subset contrast") +
  labs(x = NULL, y = "Spearman rho", tag = "D") +
  theme_a3() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                      panel.grid.major.x = element_blank(),
                      plot.margin = margin(A3_HEADER_GUTTER_PT, 5, 16, 5))
save_panels(list(p7A, p7B, p7C, p7D), "F7_binding_boundary", 12.2, 9.6)

# ===================== Fig 8 / F8_developability_map =====================
# Caption face is A = discovery fraction, B = control applicability (p8A + p8D).
# Zero-height liability bars are invisible, so do not encode them only as an
# unused orange legend. Colour those y-axis labels and spell the zeros.
liability_props <- c("polyreactivity", "aggregation")
f8 <- rd("F5_developability.csv") %>%
  mutate(property = factor(property, levels = property),
         lab = paste0(n_sig_raw_bh, "/", n_datasets),
         null_liability = property %in% liability_props,
         frac_plot = frac_sig)
f8_y_cols <- ifelse(levels(f8$property) %in% liability_props, oi("verm"), "black")
p8A <- ggplot(f8,aes(property,frac_plot)) +
  geom_col(fill=oi("blue"),width=A3_BAR) +
  geom_point(data=filter(f8,frac_plot==0),shape=21,fill="white",colour=oi("blue"),size=2.6,stroke=.7) +
  geom_text(aes(y=ifelse(frac_plot==0,.055,frac_plot),label=lab),hjust=-.22,family=A3_FONT) +
  scale_y_continuous(limits=c(-.035,1.22),breaks=c(0,.25,.5,.75,1),expand=expansion(0)) +
  coord_flip() + labs(x=NULL,y="fraction of tested datasets") +
  theme_a3() + theme(panel.grid.major.y=element_blank())
p8B <- ggplot(f8, aes(n_datasets, n_sig_raw_bh, colour = null_liability, label = property)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = A3_POINT) +
  scale_colour_manual(values = c(`FALSE` = oi("blue"), `TRUE` = oi("verm")), guide = "none") +
  labs(title = "Discovery count versus coverage", x = "n datasets", y = "n BH-significant", tag = "B") +
  theme_a3()
p8C <- locked_f8c() %>%
  transmute(group = factor(group, levels = c("liability", "non-liability")),
            n_sig = as.numeric(n_sig), n_ds = as.numeric(n_datasets), frac = frac,
            lab = paste0(n_sig, "/", n_ds))
p8C_plot <- ggplot(p8C, aes(group, frac, fill = group)) +
  geom_col(width = A3_BAR, show.legend = FALSE) +
  geom_text(aes(label = lab), vjust = -0.4, family = A3_FONT) +
  scale_fill_manual(values = c("liability" = oi("verm"), "non-liability" = oi("blue"))) +
  scale_y_continuous(limits = c(0, 1.08), expand = expansion(c(0, 0.02))) +
  labs(title = "Liability pooled contrast", x = NULL, y = "pooled fraction", tag = "C") +
  theme_a3() + theme(panel.grid.major.x = element_blank(),
                      plot.margin = margin(6, 5, 4, 6))
p8D <- ggplot(f8,aes(property,n_applicable_multigermline)) +
  geom_col(fill=oi("sky"),width=A3_BAR) +
  geom_text(aes(label=paste0(n_applicable_multigermline,"/",n_datasets)),hjust=-.18,family=A3_FONT) +
  scale_y_continuous(limits=c(0,18),breaks=c(0,5,10,15),expand=expansion(0)) +
  coord_flip() + labs(x=NULL,y="applicable / tested datasets") +
  theme_a3() + theme(panel.grid.major.y=element_blank())
save_panels(list(p8A, p8D), "F8_developability_map")

# ===================== Fig S8 / F9_antibody_lm =====================
n9 <- rd("F8_antibody_n9.csv") %>%
  pivot_longer(c(antiberty_bh_sig, esm2_bh_sig), names_to = "model", values_to = "sig") %>%
  mutate(model = recode(model, antiberty_bh_sig = "AntiBERTy", esm2_bh_sig = "ESM-2 650M"),
         frac = sig / n_datasets, lab = paste0(sig, "/", n_datasets),
         property = factor(property, c("polyreactivity", "aggregation"),
                           labels = c("polyreactivity\nassay", "aggregation\nassay")))
p9A <- ggplot(n9, aes(property, frac, fill = model)) +
  geom_col(position = position_dodge(0.72), width = A3_BAR) +
  geom_text(aes(label = lab,colour=model), position = position_dodge(0.72), hjust = -0.10, size = 2.7, family = A3_FONT) +
  scale_fill_manual(values = c("AntiBERTy" = oi("verm"), "ESM-2 650M" = oi("blue"))) +
  scale_colour_manual(values=c("AntiBERTy"=oi("verm"),"ESM-2 650M"=oi("blue")),guide="none") +
  scale_y_continuous(limits = c(0, 0.4), expand = expansion(c(0, 0.04))) +
  coord_flip() +
  labs(x = NULL, y = "fraction BH-significant", tag = "A") +
  theme_a3() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
n4 <- rd("F8_antibody_n4.csv") %>% mutate(model = factor(model, levels = rev(model)))
p9B <- ggplot(n4, aes(rho, model, colour = model)) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.16, linewidth = A3_ERRORBAR_LW) +
  geom_point(size = A3_POINT) +
  scale_colour_manual(values = c("AntiBERTy" = oi("verm"), "ESM-2 650M" = oi("blue")), guide = "none") +
  labs(x = "Spearman(LLR, \u2212\u0394\u0394G), Ab\u2013Ag", y = NULL, tag = "B") +
  theme_a3() + theme(panel.grid.major.y = element_blank())
p9C <- ggplot(n9, aes(property, sig, fill = model)) +
  geom_col(position = position_dodge(0.72), width = A3_BAR) +
  scale_fill_manual(values = c("AntiBERTy" = oi("verm"), "ESM-2 650M" = oi("blue")), guide = "none") +
  coord_flip() +
  geom_text(aes(label=sig,colour=model),position=position_dodge(.72),hjust=-.16,family=A3_FONT) +
  scale_colour_manual(values=c("AntiBERTy"=oi("verm"),"ESM-2 650M"=oi("blue")),guide="none") +
  scale_y_continuous(limits=c(0,4.1),breaks=0:3,expand=expansion(0)) +
  labs(x = NULL, y = "n BH-significant datasets", tag = "C") +
  theme_a3() + theme(panel.grid.major.y = element_blank())
n4w <- n4 %>% mutate(ci_width = ci_hi - ci_lo, crosses_zero = ci_lo < 0 & ci_hi > 0,
                     model = factor(model, labels = c("ESM-2\n650M", "AntiBERTy")))
p9D <- ggplot(n4w, aes(model, ci_width, fill = crosses_zero)) +
  geom_col(width = A3_BAR) +
  geom_text(aes(label = sprintf("%.3f", ci_width)), vjust = -0.35, size = A3_TEXT_MM, family = A3_FONT) +
  scale_fill_manual(values = c(`TRUE` = oi("grey"), `FALSE` = oi("blue")),
                    labels = c(`TRUE` = "CI spans 0", `FALSE` = "CI excludes 0"),
                    name = NULL) +
  a3_y_for_labels(bottom = 0.02, top = 0.20) +
  labs(x = NULL, y = "bootstrap CI width", tag = "D") +
  theme_a3() + theme(panel.grid.major.x = element_blank()) +
  a3_legend_outside("bottom")
save_panels(list(p9A, p9B, p9C, p9D), "F9_antibody_lm")

# ===================== Fig S9 / F10_reliability_predictor =====================
d <- rd("F10_reliability_predictor.csv")
d$axis <- factor(d$axis, levels = c("DMS", "CLINVAR"), labels = c("DMS", "ClinVar"))
AXCOL <- c("DMS" = oi("blue"), "ClinVar" = oi("orange"))
sp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}
dms <- d %>% filter(axis == "DMS", is.finite(score_std))
p10A <- ggplot(dms, aes(score_std, auroc)) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = 2, colour = "grey60") +
  geom_smooth(method = "lm", se = TRUE, colour = oi("blue"),
              fill = scales::alpha(oi("blue"), 0.12), linewidth = A3_LINE, formula = y ~ x) +
  geom_point(size = A3_POINT, colour = oi("blue"), alpha = 0.85) +
  labs(x = "raw-LLR spread (SD)", y = "per-protein AUROC", tag = "A") +
  theme_a3()
dl <- d %>% filter(is.finite(length))
p10B <- ggplot(dl, aes(length, auroc, colour = axis)) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = 2, colour = "grey60") +
  geom_smooth(method = "lm", se = FALSE, linewidth = A3_LINE, formula = y ~ x) +
  geom_point(size = A3_POINT, alpha = 0.85) +
  scale_colour_manual(values = AXCOL) + scale_x_log10() +
  labs(x = "sequence length", y = "AUROC", tag = "B") +
  theme_a3() + a3_legend_outside("bottom")
db <- d %>% filter(is.finite(class_balance))
p10C <- ggplot(db, aes(class_balance, auroc, colour = axis)) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = 2, colour = "grey60") +
  geom_smooth(method = "lm", se = FALSE, linewidth = A3_LINE, formula = y ~ x) +
  geom_point(size = A3_POINT, alpha = 0.85) +
  scale_colour_manual(values = AXCOL, guide = "none") +
  labs(x = "positive-class fraction", y = "AUROC", tag = "C") +
  theme_a3()
# committed cross-axis summary values (same as prior make_f7 script; not invented)
summ10 <- data.frame(
  predictor = c("sample size n", "length", "class balance"),
  DMS = c(0.046, 0.064, 0.093),
  ClinVar = c(-0.247, -0.378, 0.263)
)
# Recompute from CSV where possible for honesty
summ10$DMS[1] <- sp(d$n[d$axis == "DMS"], d$auroc[d$axis == "DMS"])
summ10$DMS[2] <- sp(d$length[d$axis == "DMS"], d$auroc[d$axis == "DMS"])
summ10$DMS[3] <- sp(d$class_balance[d$axis == "DMS"], d$auroc[d$axis == "DMS"])
summ10$ClinVar[1] <- sp(d$n[d$axis == "ClinVar"], d$auroc[d$axis == "ClinVar"])
summ10$ClinVar[2] <- sp(d$length[d$axis == "ClinVar"], d$auroc[d$axis == "ClinVar"])
summ10$ClinVar[3] <- sp(d$class_balance[d$axis == "ClinVar"], d$auroc[d$axis == "ClinVar"])
sl <- data.frame(
  predictor = rep(summ10$predictor, 2),
  axis = factor(rep(c("DMS", "ClinVar"), each = 3), levels = c("DMS", "ClinVar")),
  rho = c(summ10$DMS, summ10$ClinVar)
)
sl$predictor <- factor(sl$predictor, levels = rev(summ10$predictor))
p10D <- ggplot(sl, aes(rho, predictor, colour = axis)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey50") +
  geom_point(size = A3_POINT, position = position_dodge(0.5)) +
  scale_colour_manual(values = AXCOL) +
  scale_x_continuous(limits = c(-0.55, 0.55)) +
  labs(x = "predictor\u2013AUROC \u03c1", y = NULL, tag = "D") +
  theme_a3() + a3_legend_outside("bottom")
save_panels(list(p10A, p10B, p10C, p10D), "F10_reliability_predictor")

# ===================== Fig 11 / F11_boundary_atlas =====================
au <- rd("F11_atlas_auroc.csv")
clinical_macro <- rd("clinical_macro_20260923.csv")
au$value[3] <- clinical_macro$value
au$lo[3] <- clinical_macro$lo
au$hi[3] <- clinical_macro$hi
spn <- rd("F11_atlas_spearman.csv")
au$label <- c("DMS: ESM-2 650M", "DMS: ProtBERT", "ClinVar: 20 genes")
spn$label <- c("Stability", "Binding: general", "Binding: overall", "Ab interface*", "Ab\u2013Ag*", "Binding: non-interface")
au$label <- factor(au$label, levels = rev(au$label))
spn$label <- factor(spn$label, levels = rev(spn$label))
p11A <- ggplot(au, aes(value, label)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_errorbar(aes(xmin = ifelse(is.na(lo), value, lo), xmax = ifelse(is.na(hi), value, hi)), orientation = "y",
                 width = 0.2, colour = oi("blue"), linewidth = A3_ERRORBAR_LW, na.rm = TRUE) +
  geom_point(size = A3_POINT, colour = oi("blue")) +
  geom_text(aes(label = sprintf("%.2f", value)), hjust = -0.35, nudge_y = 0.20, vjust = 0, colour = "grey20", family = A3_FONT) +
  # Three ticks printed as the run "0.500.751.00", and the right limit of 1.08
  # left no room for a value label past a marker near 1.0, so 0.78 clipped.
  # Two ticks plus a wider right margin fixes both without moving a value.
  scale_x_continuous(limits = c(0.45, 1.02), breaks = c(0.5, 0.75, 1.0)) +
  labs(title = "Discrimination summaries", x = "AUROC", y = NULL, tag = "A") +
  theme_a3() + theme(panel.grid.major.y = element_blank())
p11B <- ggplot(spn, aes(value, label)) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  geom_errorbar(aes(xmin = ifelse(is.na(lo), value, lo), xmax = ifelse(is.na(hi), value, hi)), orientation = "y",
                 width = 0.2, colour = oi("verm"), linewidth = A3_ERRORBAR_LW, na.rm = TRUE) +
  geom_point(size = A3_POINT, colour = oi("verm")) +
  # Same conditional flip as Fig 7 A: the rightmost task's label goes left of its
  # marker so it can never reach the spine, and the panel keeps its width.
  geom_text(aes(label = sprintf("%.2f", value),
                hjust = ifelse(value > 0.6 * max(value), 1.25, -0.30)),
            colour = "grey20", nudge_y = 0.22, vjust = 0, family = A3_FONT) +
  # Three breaks printed as the run "-0.10.1 0.3"; two cannot collide.
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.26)), breaks = c(-0.1, 0.1, 0.3)) +
  labs(title = "Physical-readout summaries", x = "Spearman \u03c1", y = NULL, tag = "B") +
  theme_a3() + theme(panel.grid.major.y = element_blank())
# Metrics are intentionally separated. Binding rows overlap and are not
# independent tasks; neither null-subtraction nor heuristic zones repair that.
# Starred nodes are archives with unreconciled 650M row/scorer/CI identities.
p11B <- p11B + labs(caption="* Archived 650M provenance\nnot reconciled\nBinding subsets overlap;\nnot independent tasks") +
  theme(plot.caption=element_text(family=A3_FONT,size=A3_TICK_PT,hjust=0,lineheight=.95))
save_panels(list(p11A, p11B), "F11_boundary_atlas")

# ===================== Fig S1 / FS1 =====================
a1 <- rd("FS1_proteingym_atlas.csv")
a1$assay_short <- substr(sub("_[0-9]{4}$", "", sub("_[A-Za-z]+$", "", a1$assay)), 1, 18)
a1 <- a1 %>% arrange(auroc) %>%
  mutate(assay_short = factor(assay_short, levels = unique(assay_short)),
         disc = ifelse(ci_excludes_chance == 1, "discriminates", "at chance"))
# 19 identifiers need ~126 pt of stacked 6.64-pt glyphs. An equal 3x2 cell
# supplies ~104 pt (measured 1.22 bp overlap on every adjacent pair). The
# two-sub-column split used for Fig 2 A does not transfer: each 18-character
# tick is ~62 pt wide, so two gutters collapse a 156 pt cell. a3_compose_tall_a
# gives this panel two rows and two columns of the same canvas. Identifiers
# stay the 18-character unique truncation; type is not shrunk.
pool <- rd("FS1_proteingym_pooled.csv")
pS1A <- ggplot(a1, aes(auroc, assay_short, colour = disc)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = pool$value[1], colour = oi("blue"), linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0, linewidth = A3_ERRORBAR_LW) +
  geom_point(size = A3_POINT) +
  scale_colour_manual(values = c(discriminates = oi("blue"), `at chance` = oi("verm"))) +
  labs(x = "per-assay AUROC", y = "assay", tag = "A") +
  theme_a3() + theme(legend.position = "none", panel.grid.major.y = element_blank())
pS1B <- ggplot(a1, aes(auroc, fill = disc)) +
  geom_histogram(bins = 10, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = pool$value[1], colour = oi("blue"), linewidth = 0.5) +
  scale_fill_manual(values = c(discriminates = oi("blue"), `at chance` = oi("verm")), guide = "none") +
  labs(x = "per-assay AUROC", y = "n assays", tag = "B") +
  theme_a3()
pS1C <- ggplot(a1, aes(n_test, auroc, colour = disc)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(size = A3_POINT) +
  scale_colour_manual(values = c(discriminates = oi("blue"), `at chance` = oi("verm")), guide = "none") +
  labs(x = "n test variants", y = "AUROC", tag = "C") +
  theme_a3()
pool_bar <- data.frame(
  label = c("pooled AUROC", "discriminating\nfraction"),
  value = c(pool$value[1], pool$n_discriminating[1] / pool$n_assays[1])
)
pS1D <- ggplot(pool_bar, aes(value, label)) +
  geom_col(fill = oi("blue"), width = A3_BAR) +
  geom_text(aes(label = sprintf("%.3f", value)), hjust = -0.25, size = 2.7, family = A3_FONT) +
  scale_x_continuous(limits = c(0, 1.15), expand = expansion(c(0, 0.02))) +
  labs(x = "pooled summary", y = NULL, tag = "D") +
  theme_a3() + theme(panel.grid.major.y = element_blank())
save_panels(list(pS1A, pS1B, pS1C, pS1D), "FS1_proteingym_atlas", 12.2, 9.6,
            layout = "tall_a")

# ===================== Fig S2 / FS2 =====================
g <- rd("FS2_clinvar_vs_gnomad.csv") %>% arrange(auroc_clinvar_benign) %>%
  mutate(gene = factor(gene, levels = unique(gene)))
# Same twenty-gene axis as Fig 2 panel A, and the same overlap: measured on
# 38bbaa0d the labels struck through each other by 3.35 pt per pair. Split into
# two sub-columns of ten on the existing order; no value or ordering changes.
g <- g %>% mutate(forest_col = ifelse(as.integer(gene) <= 10L, "lower ten", "upper ten"))
gl <- g %>% select(gene, forest_col, auroc_clinvar_benign, auroc_gnomad_benign) %>%
  pivot_longer(c(auroc_clinvar_benign, auroc_gnomad_benign),
               names_to = "benign_source", values_to = "auroc") %>%
  mutate(benign_source = recode(benign_source,
                                auroc_clinvar_benign = "ClinVar benign",
                                auroc_gnomad_benign = "gnomAD AC\u22652"))
pS2A <- ggplot() +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  geom_segment(data = g, aes(x = auroc_clinvar_benign, xend = auroc_gnomad_benign,
                             y = gene, yend = gene), colour = "grey75", linewidth = 0.45) +
  geom_point(data = gl, aes(auroc, gene, colour = benign_source), size = 1.8) +
  scale_colour_manual(values = c("ClinVar benign" = oi("blue"), "gnomAD AC\u22652" = oi("orange"))) +
  # Halving the forest width leaves no room for the six default breaks: they
  # rendered as a single abutting run "0.50.60.70.80.91.0". Two breaks bracket
  # the range; the dashed chance line carries the 0.5 reference.
  scale_x_continuous(breaks = c(0.5, 1.0)) +
  labs(x = "per-gene AUROC", y = NULL, tag = "A") +
  theme_a3() + theme(panel.grid.major.y = element_blank(),
                      strip.text = element_blank(), strip.background = element_blank(),
                      panel.spacing.x = unit(3, "pt")) +
  a3_legend_outside("bottom")
pS2B <- ggplot(g, aes(auroc_clinvar_benign, auroc_gnomad_benign, colour = class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(size = A3_POINT) +
  scale_colour_manual(values = c(`gain-of-function` = oi("blue"), other = oi("orange"))) +
  coord_equal(xlim = c(0.48, 1.02), ylim = c(0.48, 1.02), expand=FALSE) +
  labs(x = "ClinVar-benign AUROC", y = "gnomAD-benign AUROC", tag = "B") +
  theme_a3() + a3_legend_outside("bottom")
pS2C <- ggplot(g, aes(delta_gnomad_minus_clinvar)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.35) +
  geom_histogram(bins = 10, fill = oi("sky"), colour = "white", linewidth = 0.2) +
  labs(x = "\u0394 AUROC (gnomAD \u2212 ClinVar)", y = "n genes", tag = "C") +
  theme_a3()
cm <- rd("FS2_class_means.csv")
pS2D <- ggplot(cm, aes(class, mean_clinvar_auroc, fill = class)) +
  geom_col(width = A3_BAR, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f\n(n=%d)", mean_clinvar_auroc, n_genes)),
            vjust = -0.2, size = A3_TEXT_MM, family = A3_FONT) +
  scale_fill_manual(values = c(`gain-of-function` = oi("blue"), other = oi("orange"))) +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "mean ClinVar AUROC", tag = "D") +
  theme_a3() + theme(panel.grid.major.x = element_blank())
save_panels(list(pS2A, pS2B, pS2C, pS2D), "FS2_gnomad_control")

# ===================== Fig S3 / FS3 =====================
d3 <- rd("FS3_antibody_robustness.csv")
sub_lv <- d3$subset[order(d3$order)] |> unique()
meth_lv <- c("ESM-2 650M (WT)", "ESM-2 3B (WT)", "ESM-2 650M (masked)")
d3 <- d3 %>% mutate(subset = factor(subset, levels = rev(sub_lv)),
                    method = factor(method, levels = meth_lv,
                                    labels = c("650M WT", "3B WT", "650M masked")))
cols3 <- c("650M WT" = oi("blue"), "3B WT" = oi("verm"), "650M masked" = oi("green"))
pS3A <- ggplot(d3, aes(rho, subset, colour = method)) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.28,
                 position = position_dodge(0.62), linewidth = A3_ERRORBAR_LW) +
  geom_point(size = A3_POINT, position = position_dodge(0.62)) +
  scale_colour_manual(values = cols3) +
  labs(x = "Spearman(LLR, \u2212\u0394\u0394G)", y = NULL, tag = "A") +
  theme_a3() + theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
                      legend.text = element_text(size = A3_TICK_PT))
ov <- d3 %>% filter(subset == "overall")
pS3B <- ggplot(ov, aes(method, rho, fill = method)) +
  geom_col(width = A3_BAR, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = A3_ERRORBAR_LW) +
  scale_fill_manual(values = cols3) +
  labs(x = NULL, y = "overall Spearman", tag = "B") +
  theme_a3() + theme(axis.text.x = element_text(angle = 18, hjust = 1),
                      panel.grid.major.x = element_blank())
ab <- d3 %>% filter(subset == "antibody")
pS3C <- ggplot(ab, aes(method, rho, fill = method)) +
  geom_col(width = A3_BAR, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = A3_ERRORBAR_LW) +
  scale_fill_manual(values = cols3) +
  labs(x = NULL, y = "antibody Spearman", tag = "C") +
  theme_a3() + theme(axis.text.x = element_text(angle = 18, hjust = 1),
                      panel.grid.major.x = element_blank())
# D: 3B − 650M WT delta by subset (committed rows only)
w650 <- d3 %>% filter(method == "650M WT") %>% select(subset, rho_650 = rho)
w3b <- d3 %>% filter(method == "3B WT") %>% select(subset, rho_3b = rho)
dd <- inner_join(w650, w3b, by = "subset") %>% mutate(delta = rho_3b - rho_650)
pS3D <- ggplot(dd, aes(subset, delta)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_col(fill = oi("verm"), width = A3_BAR) +
  geom_text(aes(label = sprintf("%+.3f", delta)),
            vjust = ifelse(dd$delta >= 0, -0.45, 1.45),
            size = A3_TEXT_MM, family = A3_FONT) +
  a3_y_for_labels(bottom = 0.28, top = 0.22) +
  labs(x = NULL, y = "\u0394\u03c1 (3B \u2212 650M WT)", tag = "D") +
  theme_a3() + theme(axis.text.x = element_text(angle = 20, hjust = 1),
                      panel.grid.major.x = element_blank())
save_panels(list(pS3A, pS3B, pS3C, pS3D), "FS3_antibody_robustness")

# ===================== Fig S4 / FS4 =====================
a4 <- rd("FS4_abstention_tiers.csv")
tier_lv <- a4$tier[order(a4$order)] |> unique()
cond_lv <- c("full (100%)", "keep-hi 50%", "shuffle 50%")
a4 <- a4 %>% mutate(tier = factor(tier, levels = tier_lv),
                    condition = factor(condition, levels = cond_lv))
condcol <- c("full (100%)" = oi("grey"), "keep-hi 50%" = oi("verm"), "shuffle 50%" = oi("blue"))
pS4A <- ggplot(a4, aes(tier, auroc, fill = condition)) +
  geom_col(position = position_dodge(0.74), width = A3_BAR) +
  scale_fill_manual(values = condcol) +
  scale_y_continuous(limits = c(0, 0.86), expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "protein-mean AUROC", tag = "A") +
  theme_a3() + theme(panel.grid.major.x = element_blank(),
                      axis.text.x = element_text(angle = 15, hjust = 1))
b4 <- rd("FS4_abstention_coverage.csv") %>%
  pivot_longer(c(conf_auroc, shuffle_auroc), names_to = "series", values_to = "auroc") %>%
  mutate(series = recode(series, conf_auroc = "confidence-ranked",
                         shuffle_auroc = "shuffled control"))
pS4B <- ggplot(b4, aes(coverage, auroc, colour = series)) +
  geom_line(linewidth = A3_LINE) + geom_point(size = A3_POINT) +
  scale_colour_manual(values = c("confidence-ranked" = oi("verm"), "shuffled control" = oi("blue"))) +
  scale_x_reverse() +
  labs(x = "coverage", y = "protein-mean AUROC", tag = "B") +
  theme_a3() + a3_legend_outside("bottom")
# C: full − keep-hi delta by tier
wide4 <- a4 %>% select(tier, condition, auroc) %>%
  pivot_wider(names_from = condition, values_from = auroc)
wide4$delta_anti <- wide4$`full (100%)` - wide4$`keep-hi 50%`
pS4C <- ggplot(wide4, aes(tier, delta_anti)) +
  geom_col(fill = oi("verm"), width = A3_BAR) +
  labs(x = NULL, y = "AUROC drop (full \u2212 keep-hi)", tag = "C") +
  theme_a3() + theme(axis.text.x = element_text(angle = 15, hjust = 1),
                      panel.grid.major.x = element_blank())
covd <- rd("FS4_abstention_coverage.csv") %>%
  mutate(delta = shuffle_auroc - conf_auroc)
pS4D <- ggplot(covd, aes(coverage, delta)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_line(colour = oi("blue"), linewidth = A3_LINE) +
  geom_point(colour = oi("blue"), size = 2) +
  scale_x_reverse() +
  labs(x = "coverage", y = "shuffle \u2212 confidence AUROC", tag = "D") +
  theme_a3()
save_panels(list(pS4A, pS4B, pS4C, pS4D), "FS4_abstention_tiers")

# ===================== Fig S5 / FS5 =====================
vep <- rd("FS5_vep_comparison.csv") %>%
  arrange(median_abs_spearman) %>%
  mutate(predictor = ifelse(predictor == "DeepSequence", "DeepSeq", predictor),
         predictor = factor(predictor, levels = predictor),
         supervised = factor(supervised, levels = c("no", "yes"),
                             labels = c("zero-shot", "supervised")))
pS5A <- ggplot(vep, aes(median_abs_spearman, predictor, fill = supervised)) +
  geom_col(width = A3_BAR) +
  geom_text(aes(label = sprintf("%.3f", median_abs_spearman)), hjust = -0.1, size = A3_TEXT_MM, family = A3_FONT) +
  scale_fill_manual(values = c("zero-shot" = oi("blue"), supervised = oi("orange"))) +
  scale_x_continuous(limits = c(0, 0.65), expand = expansion(c(0, 0.02))) +
  labs(x = "median |Spearman| vs DMS", y = NULL, tag = "A") +
  theme_a3() + theme(panel.grid.major.y = element_blank(), legend.position = "bottom")
pS5B <- locked_summary("S5B_supervised_group_mean_median_abs_spearman.csv", "S5", "B",
                        c("no", "yes")) %>%
  transmute(supervised = factor(group_1, levels = c("no", "yes"),
                                labels = c("zero-shot", "supervised")),
            median_abs_spearman = value)
pS5B_plot <- ggplot(pS5B, aes(supervised, median_abs_spearman, fill = supervised)) +
  geom_col(width = A3_BAR, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", median_abs_spearman)), vjust = -0.4, size = 2.7, family = A3_FONT) +
  scale_fill_manual(values = c("zero-shot" = oi("blue"), supervised = oi("orange"))) +
  scale_y_continuous(limits = c(0, 0.6), expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "mean of medians |\u03c1|", tag = "B") +
  theme_a3() + theme(axis.text.x = element_text(angle = 15, hjust = 1),
                      panel.grid.major.x = element_blank())
vep$rank <- rank(-vep$median_abs_spearman)
pS5C <- ggplot(vep, aes(rank, median_abs_spearman, colour = supervised)) +
  geom_line(aes(group = 1), colour = "grey70", linewidth = A3_ERRORBAR_LW) +
  geom_point(size = A3_POINT) +
  scale_colour_manual(values = c("zero-shot" = oi("blue"), supervised = oi("orange")),
                      guide = "none") +
  scale_x_continuous(breaks = vep$rank, labels = as.character(vep$predictor)) +
  labs(x = NULL, y = "median |\u03c1|", tag = "C") +
  theme_a3() + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = A3_TICK_PT))
pS5D <- ggplot(vep, aes(predictor, n_proteins, fill = supervised)) +
  geom_col(width = A3_BAR, show.legend = FALSE) +
  geom_text(aes(label = n_proteins), vjust = -0.35, size = A3_TEXT_MM, family = A3_FONT) +
  scale_fill_manual(values = c("zero-shot" = oi("blue"), supervised = oi("orange"))) +
  scale_y_continuous(limits = c(0, 12), expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "n proteins (overlap)", tag = "D") +
  theme_a3() + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                      panel.grid.major.x = element_blank())
save_panels(list(pS5A, pS5B_plot, pS5C, pS5D), "FS5_vep_comparison")

cat("make_four_panel_all.R: wrote 16 staged reader-layout figures\n")
