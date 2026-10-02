# Data-backed evidence panels for the A3 six-panel compositor.
#
# The E/F cells report locked scalar sensitivity estimands.  Rendering each one
# as a single marker spends a sixth of every figure on one number and prints
# "no interval" in 23 of 24 READY cells, which reads as missing work rather than
# as the robustness evidence it actually is.
#
# Every locked scalar is a summary of a committed table that the E/F row already
# names in `source_path` and pins in `source_sha256`.  These renderers draw that
# table underneath the scalar.  The reported number is still the committed
# payload value -- it is never recomputed here -- but the reader now sees the
# constituent observations it summarises.
#
# Fail-closed contract: a READY estimand must have a registry entry, its source
# file must hash to the pinned digest, and the drawn observations must be
# consistent with the committed value.  Any breach stops the build.
suppressMessages({
  library(ggplot2)
  library(digest)
})

EV_INK <- "grey20"      # primary text
EV_MUTED <- "grey38"    # explanatory prose on blocked evidence-boundary cards
EV_RULE <- "grey72"     # reference rules
EV_HAIR <- "grey88"     # backdrop marks

# A locked interval, a contrast level and an interval endpoint are data, not
# footnotes: Fig 1's caption commits to printing this figure's single named
# interval "exactly as panel E prints it".  These annotations therefore take the
# primary ink at near-body size and stay subordinate to the headline scalar by
# weight and size ratio alone, not by being faded out.
EV_ANNOT_MM <- A3_TEXT_MM   # numeric secondary annotations
EV_LEVEL_MM <- A3_TEXT_MM   # the categorical level a number belongs to

# ---------------------------------------------------------------- data access

# Load the exact table the E/F row pins.  A digest mismatch means the drawn
# observations would no longer be the ones behind the committed scalar.
.ev_load <- function(row) {
  src <- basename(row$source_path[[1]])
  path <- file.path(PD, src)
  if (!file.exists(path)) {
    stop(sprintf("evidence panel source missing: %s", src), call. = FALSE)
  }
  pinned <- row$source_sha256[[1]]
  actual <- digest::digest(file = path, algo = "sha256")
  if (!identical(actual, pinned)) {
    stop(sprintf("evidence panel source digest drift for %s (%s)", src, row$estimand_id[[1]]),
         call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# The drawn observations must reproduce the committed scalar.  Tolerance is
# loose enough for the rounding actually present in the committed payloads
# (several are stored to 3--4 significant figures) and tight enough to catch a
# wrong column or a wrong subset.
.ev_agrees <- function(derived, committed, tol = 5e-3) {
  is.finite(derived) && is.finite(committed) && abs(derived - committed) <= tol
}

.ev_check <- function(derived, committed, id, tol = 5e-3) {
  if (!.ev_agrees(derived, committed, tol)) {
    stop(sprintf("evidence panel disagrees with committed value for %s: drawn %.6g vs committed %.6g",
                 id, derived, committed), call. = FALSE)
  }
  invisible(TRUE)
}

# ------------------------------------------------------------------- helpers

# Deterministic dot stacking.  Random jitter would break byte-reproducible
# rendering, so observations are binned and stacked in a fixed order.
.ev_stack <- function(x, nbins = 20L) {
  x <- as.numeric(x)
  rng <- range(x, finite = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) <= 0) return(seq_along(x))
  edges <- seq(rng[1], rng[2], length.out = nbins + 1L)
  b <- as.integer(cut(x, breaks = edges, include.lowest = TRUE, labels = FALSE))
  y <- integer(length(x))
  seen <- integer(nbins + 1L)
  for (i in order(b, x)) {
    seen[b[i]] <- seen[b[i]] + 1L
    y[i] <- seen[b[i]]
  }
  y
}

# A single extreme observation can claim most of the axis and squeeze the
# informative mass into a spike a few points wide.  Dropping, clipping or
# winsorising it would be dishonest, and a smooth compressing transform
# (asinh, signed-sqrt) would silently move it nearer the mass without telling
# the reader.  So the axis is broken instead: bulk and tail get separate,
# separately-ticked linear segments with a physical gap between them, every
# observation stays drawn at its own value, and the jump in the tick sequence
# declares the discontinuity.  Distance is never comparable across the gap --
# which is exactly the property a broken axis asserts.
#
# The rule is deterministic and deliberately reluctant: break only when the
# largest gap in the sorted observations is a genuine chasm and only a handful
# of observations sit beyond it.  Otherwise return NA and draw one linear axis.
.ev_tail_cut <- function(x, min_gap_frac = 0.30, max_tail_frac = 0.06,
                         min_n = 12L) {
  x <- sort(as.numeric(x[is.finite(x)]))
  n <- length(x)
  span <- if (n) diff(range(x)) else 0
  if (n < min_n || !is.finite(span) || span <= 0) return(NA_real_)
  gaps <- diff(x)
  i <- which.max(gaps)
  if (gaps[i] / span < min_gap_frac || (n - i) / n > max_tail_frac) return(NA_real_)
  mean(c(x[i], x[i + 1L]))
}

# Choose an inventory width that leaves the last row as full as possible; a
# single orphan square reads as a rendering fault rather than as a count.
.ev_grid_ncol <- function(n) {
  cand <- 4:12
  cand <- cand[cand <= max(4L, n)]
  rem <- n %% cand
  fullness <- ifelse(rem == 0, 1, rem / cand)
  target <- sqrt(n * 1.6)
  best <- order(-fullness, abs(cand - target))[1]
  as.integer(cand[best])
}

.ev_fmt <- function(v, digits = 3) {
  if (!is.finite(v)) return("--")
  if (abs(v) >= 1000 || (abs(v) > 0 && abs(v) < 1e-3)) {
    formatC(v, format = "g", digits = digits)
  } else {
    formatC(v, format = "fg", digits = digits, flag = "#")
  }
}

.ev_signed <- function(v, digits = 3) {
  s <- .ev_fmt(abs(v), digits)
  paste0(if (v < 0) "\u2212" else "+", s)
}

# Shared chrome for every evidence panel: the committed scalar sits in a fixed
# top-left slot and the sample descriptor in a fixed top-right slot.  Consistent
# placement lets a reader scan the sensitivity cells across all 16 figures.
#
# The E/F cells are roughly two inches wide, so a value and an interval cannot
# share one line.  Any secondary string goes on a second, lighter line below.
#
# `facet` names the single facet the header belongs to when the panel's x axis
# is broken; without it a bare annotate() would repeat the header inside every
# facet, each anchored to its own left edge.
.ev_header <- function(p, value_label, xlim, ytop, sub_label = NULL,
                       line_drop = 0.23, facet = NULL) {
  slot <- function(y, label) {
    df <- data.frame(x = xlim[1], y = y, label = label)
    if (!is.null(facet)) df$seg <- factor(facet, levels = c("bulk", "tail"))
    df
  }
  p <- p +
    ggplot2::geom_text(data = slot(ytop, value_label),
                       aes(x = x, y = y, label = label), inherit.aes = FALSE,
                       family = A3_FONT, fontface = "plain",
                       size = A3_TEXT_MM, hjust = 0, vjust = 1,
                       colour = EV_INK)
  if (!is.null(sub_label)) {
    p <- p + ggplot2::geom_text(data = slot(ytop * (1 - line_drop), sub_label),
                                aes(x = x, y = y, label = label),
                                inherit.aes = FALSE, family = A3_FONT,
                                size = EV_ANNOT_MM, hjust = 0, vjust = 1,
                                lineheight = 0.92, colour = EV_INK)
  }
  p
}

# The sample descriptor rides on the axis title rather than inside the cell.
# In a two-inch panel any in-panel right-hand slot eventually collides with the
# value, and the descriptor is a statement about the axis anyway.
.ev_xlab <- function(units, n_label) {
  paste0(wrap_card(units, 24), "\n(", n_label, ")")
}

.ev_theme <- function(tag, xlab) {
  list(
    labs(title = NULL, x = xlab, y = NULL),
    theme_a3(),
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linewidth = 0.22, colour = EV_HAIR),
          # Top margin matches the shared header gutter: data-backed evidence
          # panels use the same margin-band tags as A--D, so they need the
          # same tag room above the panel.
          plot.margin = margin(2, 6, 3, 6))
  )
}

# --------------------------------------------------------------- shape: strip

# One-dimensional display of the constituent observations, with either the
# committed mean (and its locked interval, when one exists) or the committed
# min--max span drawn over them.
.ev_strip <- function(row, tag, values, summary = c("mean", "span"), n_label,
                      null_at = 0) {
  summary <- match.arg(summary)
  value <- as.numeric(row$value[[1]])
  lo <- suppressWarnings(as.numeric(row$interval_lo[[1]]))
  hi <- suppressWarnings(as.numeric(row$interval_hi[[1]]))
  has_ci <- is.finite(lo) && is.finite(hi)

  derived <- if (summary == "mean") mean(values) else diff(range(values))
  .ev_check(derived, value, row$estimand_id[[1]])

  cut <- .ev_tail_cut(values)
  broken <- is.finite(cut)
  seg <- if (broken) {
    factor(ifelse(values > cut, "tail", "bulk"), levels = c("bulk", "tail"))
  } else factor(rep("bulk", length(values)), levels = c("bulk", "tail"))
  # Stack inside each axis segment.  Binning over the whole range would put the
  # entire bulk into two or three bins of a range the tail defines.
  stack <- integer(length(values))
  for (level in levels(seg)) {
    keep <- seg == level
    if (any(keep)) stack[keep] <- .ev_stack(values[keep])
  }
  df <- data.frame(x = values, y = stack, seg = seg)
  ymax <- max(stack)
  # Reserve the upper part of the cell for the header slots; a two-line interval
  # sub-label needs a taller band than a bare scalar.
  ytop <- if (has_ci) ymax * 2.15 + 1.4 else ymax * 1.72 + 1.1
  # The committed summary rides in its own reserved lane just above the swarm.
  # Drawn behind the points it was invisible: a 0.031-wide CI under a 69-point
  # stack is covered by the very observations it summarises.
  yb <- ymax * 1.22 + 0.6
  ytick <- ymax * 0.10 + 0.25

  bulk <- values[seg == "bulk"]
  pts <- c(bulk, if (has_ci) c(lo, hi), if (is.finite(null_at)) null_at,
           if (summary == "mean") value)
  rng <- range(pts, finite = TRUE)
  pad <- if (diff(rng) > 0) diff(rng) * 0.10 else max(abs(value) * 0.5, 0.05)
  xlim <- c(rng[1] - pad, rng[2] + pad)
  in_bulk <- function(v) is.finite(v) & v >= xlim[1] & v <= xlim[2]

  # `bulk` rows carry the summary marks; when the axis is broken the facet
  # column restricts each mark to the segment it belongs in.
  mark <- function(...) {
    df <- data.frame(...)
    df$seg <- factor("bulk", levels = c("bulk", "tail"))
    df
  }
  p <- ggplot(df, aes(x, y))
  if (in_bulk(null_at)) {
    p <- p + geom_vline(data = mark(xintercept = null_at), aes(xintercept = xintercept),
                        inherit.aes = FALSE, linetype = "dashed",
                        colour = EV_RULE, linewidth = 0.4)
  }
  p <- p + geom_point(size = A3_POINT, colour = oi("blue"), alpha = 0.82)
  if (summary == "span") {
    # Bracket the committed span across the observations it is taken over.
    p <- p +
      geom_segment(data = mark(x = min(values), xend = max(values), y = yb, yend = yb),
                   aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                   colour = oi("blue"), linewidth = 0.55) +
      geom_segment(data = mark(x = c(min(values), max(values)),
                               xend = c(min(values), max(values)),
                               y = yb - ytick, yend = yb + ytick),
                   aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                   colour = oi("blue"), linewidth = 0.55)
  } else {
    if (has_ci) {
      p <- p + geom_rect(data = mark(xmin = lo, xmax = hi,
                                     ymin = yb - ytick, ymax = yb + ytick),
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                         inherit.aes = FALSE, fill = oi("blue"), alpha = 0.22)
    }
    p <- p + geom_segment(data = mark(x = value, xend = value,
                                      y = yb - ytick * 1.7, yend = yb + ytick * 1.7),
                          aes(x = x, xend = xend, y = y, yend = yend),
                          inherit.aes = FALSE, colour = oi("blue"), linewidth = 0.8)
  }

  p <- p + scale_y_continuous(limits = c(0, ytop * 1.02), expand = expansion(0))
  # Two inches of cell cannot hold the label and the interval on one line at a
  # size a reader will actually read, and Fig 1's caption undertakes to print
  # this interval exactly.  So it wraps rather than shrinks or truncates.
  sub <- if (has_ci) {
    sprintf("95%% CI [%s, %s]", .ev_fmt(lo, 2), .ev_fmt(hi, 2))
  } else NULL
  if (!broken) {
    p <- p + scale_x_continuous(limits = xlim, n.breaks = 4)
    return(.ev_header(p, .ev_fmt(value, 3), xlim, ytop, sub_label = sub) +
             .ev_theme(tag, .ev_xlab(sub(" range$", "", row$units[[1]]), n_label)))
  }

  # Declared per-segment ranges: expansion is zero so `space = "free_x"` divides
  # the cell in exactly the ratio of the two declared ranges, and a lone tail
  # observation still gets a segment wide enough to carry its own tick.
  tail <- values[seg == "tail"]
  tail_half <- max(diff(xlim) / 8, diff(range(tail)) * 0.6, abs(cut) * 1e-3)
  tail_lim <- mean(range(tail)) + c(-1, 1) * tail_half
  edges <- data.frame(
    seg = factor(rep(c("bulk", "tail"), each = 2), levels = c("bulk", "tail")),
    x = c(xlim, tail_lim), y = 0
  )
  tail_break <- signif(mean(range(tail)), 3)
  # Classic break marks at the two inner segment ends, so the gap cannot be
  # misread as two sub-panels sharing one continuous scale.
  nick <- diff(xlim) * 0.022
  nicks <- rbind(
    data.frame(seg = "bulk", x = xlim[2] - nick * c(1.6, 0.6),
               xend = xlim[2] - nick * c(0.6, -0.4)),
    data.frame(seg = "tail", x = tail_lim[1] + nick * c(-0.4, 0.6),
               xend = tail_lim[1] + nick * c(0.6, 1.6))
  )
  nicks$seg <- factor(nicks$seg, levels = c("bulk", "tail"))
  nicks$y <- 0
  nicks$yend <- ytop * 0.055
  p <- p +
    geom_blank(data = edges, aes(x = x, y = y), inherit.aes = FALSE) +
    geom_segment(data = nicks, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, colour = "black", linewidth = 0.4) +
    scale_x_continuous(
      expand = expansion(mult = 0.02),
      breaks = function(limits) {
        if (limits[[1]] >= cut) tail_break else scales::extended_breaks(n = 4)(limits)
      }
    ) +
    facet_grid(cols = vars(seg), scales = "free_x", space = "free_x")
  .ev_header(p, .ev_fmt(value, 3), xlim, ytop, sub_label = sub, facet = "bulk") +
    .ev_theme(tag, .ev_xlab(sub(" range$", "", row$units[[1]]), n_label)) +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          panel.spacing.x = unit(7, "pt"))
}

# ------------------------------------------------------------- shape: scatter

# Two-dimensional display for rank-association estimands.  No fit is drawn: the
# committed statistic is a rank correlation, so an overlaid least-squares line
# would assert a model the estimand does not claim.
.ev_scatter <- function(row, tag, x, y, xlab, ylab, n_label) {
  value <- as.numeric(row$value[[1]])
  derived <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  .ev_check(derived, value, row$estimand_id[[1]])
  df <- data.frame(x = x, y = y)
  xr <- range(x, finite = TRUE); yr <- range(y, finite = TRUE)
  xpad <- diff(xr) * 0.08; ypad <- diff(yr) * 0.10
  xlim <- c(xr[1] - xpad, xr[2] + xpad)
  ylim <- c(yr[1] - ypad, yr[2] + ypad * 3.2)
  ggplot(df, aes(x, y)) +
    geom_point(size = A3_POINT, colour = oi("blue"), alpha = 0.80) +
    annotate("text", x = xlim[1], y = ylim[2], family = A3_FONT, fontface = "bold",
             size = A3_TEXT_MM, hjust = 0, vjust = 1, colour = EV_INK,
             label = sprintf("\u03c1 = %s", .ev_fmt(value, 3))) +
    scale_x_continuous(limits = xlim, n.breaks = 4) +
    scale_y_continuous(limits = ylim, n.breaks = 4) +
    labs(title = NULL, x = .ev_xlab(xlab, n_label), y = ylab, tag = tag) +
    theme_a3() +
    theme(panel.grid.major = element_line(linewidth = 0.22, colour = EV_HAIR),
          plot.margin = margin(2, 6, 3, 6))
}

# ------------------------------------------------------------ shape: contrast

# Two committed group summaries with the committed difference bracketed between
# them, so the reader sees the levels rather than only the gap.
.ev_contrast <- function(row,tag,labels,values,n_label,digits=3,value_fmt=.ev_signed) {
  value <- as.numeric(row$value[[1]])
  .ev_check(values[1]-values[2],value,row$estimand_id[[1]])
  df <- data.frame(label=factor(labels,levels=rev(labels)),value=values)
  ggplot(df,aes(value,label)) +
    geom_col(fill=oi("blue"),width=.55) +
    geom_text(aes(label=sprintf("%.3f",value)),hjust=-.15,family=A3_FONT) +
    scale_x_continuous(limits=c(0,1.02),breaks=c(0,.25,.5,.75,1),expand=expansion(0)) +
    labs(x="protein-mean AUROC",y=NULL) + theme_a3() +
    theme(panel.grid.major.y=element_blank())
}

# ------------------------------------------------------------ shape: fraction

# Unit-square inventory for proportion estimands.  Each square is one committed
# record, so the denominator is visible instead of implied.
.ev_fraction <- function(row, tag, flags, n_label, ncol = NULL) {
  value <- as.numeric(row$value[[1]])
  n <- length(flags)
  k <- sum(flags)
  .ev_check(k / n, value, row$estimand_id[[1]])
  if (is.null(ncol)) ncol <- .ev_grid_ncol(n)
  # Marked records fill from the top left.  A proportion inventory carries no
  # record order, and grouping the marked squares makes the count legible at
  # thumbnail size.
  flags <- sort(flags, decreasing = TRUE)
  idx <- seq_len(n) - 1L
  df <- data.frame(col = idx %% ncol, row = -(idx %/% ncol), hit = flags)
  nrow_grid <- max(-df$row) + 1
  # coord_fixed keeps the inventory squares square; without it the tiles stretch
  # with the cell and read as bars rather than as counted records.
  ggplot(df, aes(col, row)) +
    geom_tile(aes(fill = hit), width = 0.84, height = 0.84,
              colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = c(`TRUE` = oi("blue"), `FALSE` = EV_HAIR), guide = "none") +
    annotate("text", x = -0.55, y = 1.15, family = A3_FONT, fontface = "bold",
             size = A3_TEXT_MM, hjust = 0, vjust = 0.5, colour = EV_INK,
             label = sprintf("%d / %d", k, n)) +
    coord_fixed(ratio = 1, xlim = c(-0.7, ncol - 0.3),
                ylim = c(-nrow_grid + 0.4, 1.7), clip = "off") +
    labs(title = NULL, x = .ev_xlab(sub(" range$", "", row$units[[1]]), n_label), y = NULL, tag = tag) +
    theme_a3() +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.grid.major = element_blank(),
          plot.margin = margin(2, 6, 3, 6))
}

# ------------------------------------------------------------ shape: interval

# The locked interval itself, with its committed width bracketed.
.ev_interval <- function(row, tag, lo, hi, point, n_label) {
  value <- as.numeric(row$value[[1]])
  .ev_check(hi - lo, value, row$estimand_id[[1]])
  rng <- range(c(lo, hi, 0), finite = TRUE)
  pad <- diff(rng) * 0.16
  xlim <- c(rng[1] - pad, rng[2] + pad)
  ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = EV_RULE, linewidth = 0.4) +
    annotate("rect", xmin = lo, xmax = hi, ymin = 0.82, ymax = 1.18,
             fill = oi("blue"), alpha = 0.18) +
    annotate("segment", x = lo, xend = hi, y = 1, yend = 1,
             colour = oi("blue"), linewidth = A3_LINE) +
    annotate("point", x = point, y = 1, size = A3_POINT, colour = oi("blue")) +
    annotate("segment", x = lo, xend = lo, y = 0.80, yend = 1.20,
             colour = oi("blue"), linewidth = A3_LINE) +
    annotate("segment", x = hi, xend = hi, y = 0.80, yend = 1.20,
             colour = oi("blue"), linewidth = A3_LINE) +
    annotate("text", x = lo, y = 0.62, label = .ev_fmt(lo, 2), family = A3_FONT,
             size = EV_ANNOT_MM, hjust = 0.5, colour = EV_INK) +
    annotate("text", x = hi, y = 0.62, label = .ev_fmt(hi, 2), family = A3_FONT,
             size = EV_ANNOT_MM, hjust = 0.5, colour = EV_INK) +
    annotate("text", x = xlim[1], y = 1.62, label = paste0("CI width ", .ev_fmt(value, 3)), family = A3_FONT,
             fontface = "bold", size = A3_TEXT_MM, hjust = 0, vjust = 1,
             colour = EV_INK) +
    scale_x_continuous(limits = xlim, n.breaks = 4) +
    scale_y_continuous(limits = c(0.42, 1.66), expand = expansion(0)) +
    .ev_theme(tag, .ev_xlab("Spearman rho", n_label))
}

# ------------------------------------------------------------------- registry

# One entry per READY estimand.  Each closure receives the committed E/F row and
# the pinned source table, and returns the panel.
EV_REGISTRY <- list(

  "A3-EQ-F1-PAIRED-01" = function(row, tag, d)
    .ev_strip(row, tag, d$delta_3b_minus_650m, "mean", sprintf("%d proteins", nrow(d))),

  "A3-EQ-F1-SCORE-01" = function(row, tag, d) {
    row$units <- "protein-mean AUROC"
    wt <- d$protein_mean_auroc[d$scoring == "WT-marginal"]
    mk <- d$protein_mean_auroc[d$scoring == "masked-marginal"]
    .ev_contrast(row, tag, c("WT-marg.", "masked-marg."), c(wt, mk), "217-assay scoring")
  },

  "A3-EQ-F11-SEPARATION-01" = function(row, tag, d)
    .ev_strip(row, tag, d$value, "span", sprintf("%d atlas axes", nrow(d)), null_at = NA),

  "A3-EQ-F11-PROVENANCE-01" = function(row, tag, d)
    .ev_fraction(row, tag, rep(TRUE, nrow(d)), "AUROC atlas rows", ncol = 3L),

  "A3-EQ-F3-SIZE-01" = function(row, tag, d)
    .ev_scatter(row, tag, d$wt_len, d$auroc_clinvar, "wild-type length (aa)",
                "ClinVar AUROC", sprintf("%d genes", nrow(d))),

  "A3-EQ-F4-P53-01" = function(row, tag, d) {
    r <- d[grepl("^P53", d$gene), , drop = FALSE]
    .ev_contrast(row, tag, c("ESM-2", "ESM-1v"), c(r$auroc_esm2[1], r$auroc_esm1v[1]),
                 sprintf("P53; n=%d overlap", r$n_overlap[1]))
  },

  "A3-EQ-F4-COVERAGE-01" = function(row, tag, d)
    .ev_strip(row, tag, d$accuracy, "span", sprintf("%d coverages", nrow(d)), null_at = NA),

  "A3-EQ-F5-DISC-01" = function(row, tag, d) {
    tr <- mean(d$spearman[d$split == "train"])
    ho <- mean(d$spearman[d$split == "homology-held-out"])
    .ev_contrast(row, tag, c("held-out", "train"), c(ho, tr),
                 sprintf("%d locked variants", length(unique(d$variant))))
  },

  "A3-EQ-F5-CAL-01" = function(row, tag, d) {
    tr <- mean(d$calib_rmse[d$split == "train"])
    ho <- mean(d$calib_rmse[d$split == "homology-held-out"])
    .ev_contrast(row, tag, c("held-out", "train"), c(ho, tr),
                 sprintf("%d locked variants", length(unique(d$variant))))
  },

  "A3-EQ-F7-WIDTH-01" = function(row, tag, d) {
    r <- d[d$key == "antibody_interface_ag" | grepl("ab.*ag|antigen", d$key, ignore.case = TRUE), , drop = FALSE]
    if (nrow(r) != 1L) {
      w <- d$ci_hi - d$ci_lo
      r <- d[which.min(abs(w - as.numeric(row$value[[1]]))), , drop = FALSE]
    }
    .ev_interval(row, tag, r$ci_lo[1], r$ci_hi[1], r$spearman[1],
                 sprintf("Ab\u2013Ag; n=%d", r$n[1]))
  },

  "A3-EQ-S1-CI-01" = function(row, tag, d)
    .ev_fraction(row, tag, d$ci_excludes_chance == 1, "ProteinGym assays", ncol = 7L),

  "A3-EQ-S2-GENE-01" = function(row, tag, d)
    .ev_strip(row, tag, d$delta_gnomad_minus_clinvar, "mean", sprintf("%d genes", nrow(d))),

  "A3-EQ-S2-CLASS-01" = function(row, tag, d) {
    gof <- d$mean_clinvar_auroc[d$class == "gain-of-function"]
    oth <- d$mean_clinvar_auroc[d$class != "gain-of-function"]
    .ev_contrast(row, tag, c("other", "gain-of-fn"), c(oth[1], gof[1]),
                 sprintf("%d genes", sum(d$n_genes)))
  },

  "A3-EQ-S3-PAIRED-01" = function(row, tag, d) {
    wide <- reshape(d[, c("subset", "method", "rho")], idvar = "subset",
                    timevar = "method", direction = "wide")
    cn <- colnames(wide)
    c3b <- cn[grepl("3B", cn)][1]
    c650 <- cn[grepl("650M", cn)][1]
    .ev_strip(row, tag, wide[[c3b]] - wide[[c650]], "mean",
              sprintf("%d subsets", nrow(wide)))
  },

  "A3-EQ-S4-RISK-01" = function(row, tag, d)
    .ev_strip(row, tag, d$conf_auroc - d$shuffle_auroc, "span",
              sprintf("%d coverages", nrow(d))),

  "A3-EQ-S5-RANGE-01" = function(row, tag, d)
    .ev_strip(row, tag, d$median_abs_spearman, "span",
              sprintf("%d predictors", nrow(d)), null_at = NA),

  "A3-EQ-S5-MISSING-01" = function(row, tag, d)
    .ev_fraction_count(row, tag, d$n_proteins < 10, "VEP predictors", ncol = 4L),

  "A3-EQ-S6-PAIR-01" = function(row, tag, d)
    .ev_strip(row, tag, d$protein_mean_auroc, "span", sprintf("%d models", nrow(d)),
              null_at = NA),

  "A3-EQ-S6-DIRECTION-01" = function(row, tag, d) {
    wt <- d$spearman[grepl("WT", d$comparison)]
    mk <- d$spearman[grepl("masked", d$comparison)]
    .ev_contrast(row, tag, c("vs WT", "vs masked"), c(wt[1], mk[1]), "ProtBERT agreement")
  },

  "A3-EQ-S7-SHIFT-01" = function(row, tag, d)
    .ev_strip(row, tag, d$mean_rho_unc_highshift, "span",
              sprintf("%d UQ methods", nrow(d))),

  "A3-EQ-S7-NEGATIVE-01" = function(row, tag, d)
    .ev_strip(row, tag, d$ucb - d$greedy, "mean", sprintf("%d splits", nrow(d))),

  "A3-EQ-S8-EFFECT-01" = function(row, tag, d) {
    e <- d$rho[grepl("ESM", d$model)]
    a <- d$rho[grepl("AntiBERTy", d$model)]
    .ev_contrast(row, tag, c("ESM-2", "AntiBERTy"), c(e[1], a[1]), "n=4 antibody sets")
  },

  "A3-EQ-S8-FDR-01" = function(row, tag, d)
    .ev_contrast(row, tag, c("AntiBERTy", "ESM-2"),
                 c(sum(d$antiberty_bh_sig), sum(d$esm2_bh_sig)),
                 sprintf("%d datasets", sum(d$n_datasets)), digits = 2,
                 value_fmt = function(v, digits) sprintf("%+d", as.integer(round(v)))),

  "A3-EQ-S9-PARTIAL-01" = function(row, tag, d) {
    keep <- is.finite(d$length) & is.finite(d$auroc)
    .ev_scatter(row, tag, d$length[keep], d$auroc[keep], "record length",
                "AUROC", sprintf("%d complete cases", sum(keep)))
  },

  "A3-EQ-S9-MISSING-01" = function(row, tag, d)
    .ev_fraction(row, tag, tolower(as.character(d$length_is_proxy)) == "true",
                 "reliability records")
)

# A count estimand reports the numerator, not the ratio; reuse the inventory
# display but validate against the count.
.ev_fraction_count <- function(row, tag, flags, n_label, ncol = NULL) {
  value <- as.numeric(row$value[[1]])
  .ev_check(sum(flags), value, row$estimand_id[[1]], tol = 1e-9)
  stub <- row
  stub$value[[1]] <- sum(flags) / length(flags)
  .ev_fraction(stub, tag, flags, n_label, ncol)
}

# ---------------------------------------------------------------- entry point

evidence_panel_data_backed <- function(row, tag) {
  id <- row$estimand_id[[1]]
  fn <- EV_REGISTRY[[id]]
  if (is.null(fn)) {
    stop(sprintf("READY estimand %s has no data-backed evidence panel", id), call. = FALSE)
  }
  fn(row, tag, .ev_load(row))
}
