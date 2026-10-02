# Shared A3 ggplot2 theme. Publication style, Arial only.
# Font: Arial is the required release font. A publication render may never
# silently fall back: metric changes from a substitute font invalidate both
# spacing checks and the visual review.
#
# This file owns type, colour, geom defaults, and the save helper.
# Panel letters and six-panel assembly live in compose_a3.R so ggplot
# automation cannot put tags on spines or align y-labels across siblings.
suppressMessages({
  library(ggplot2)
  library(systemfonts)
  library(ragg)
  library(scales)
  library(grid)
})

.a3_required_font <- function() {
  faces <- systemfonts::system_fonts()
  required <- faces$family == "Arial" & !faces$italic
  if (!all(c("normal", "bold") %in% as.character(faces$weight[required]))) {
    stop("A3 release rendering requires regular and bold Arial; refusing font fallback.", call. = FALSE)
  }
  "Arial"
}
A3_FONT <- .a3_required_font()

# Type scale at the printed 6.5-inch width (each of six cells ~2.1 in).
# 8.5/8.0/10 and 7.5/7.0/9 both left lettering dominant and data marks small.
# Reader-scale 8 / 7.5 / 10 pt typography; layouts allocate room before type is reduced.
A3_AXIS_PT <- 9.0
A3_TICK_PT <- 8.0
A3_TEXT_PT <- 8.0
A3_TEXT_MM <- 2.82
A3_PANEL_LABEL_PT <- 11
# Dedicated compositor band above each cell (compose_a3.R). Not a ggplot margin:
# ggplot plot.tag previously sat in this space and crowded the title/ylab.
A3_LABEL_BAND_PT <- 20
A3_LABEL_INSET_PT <- 1
A3_CELL_GUTTER_PT <- 14
# Horizontal gutters protect long y labels; row gaps need less empty space.
A3_ROW_GUTTER_PT <- 8
A3_PAGE_PAD_PT <- 2
# Retained name: some panel writers still pass it as a top inset. The
# compositor owns the letter band, so the theme itself uses a 2 pt top.
A3_HEADER_GUTTER_PT <- 2

# Data-ink defaults. Writers may still override, but the default mark is
# larger than the 2.2--2.8 pt points / 0.6-wide bars that starved the graphics.
A3_POINT <- 2.2
A3_LINE <- 0.8
A3_BAR <- 0.64
A3_ERRORBAR_LW <- 0.65
A3_ERRORBAR_W <- 0.22

# Okabe-Ito colorblind-safe palette (identical hex to the matplotlib scripts).
# oi() returns UNNAMED hex so it can be used as an element inside a named values
# vector (c(label = oi("blue"))) without R contaminating the element name.
.OI <- c(blue = "#0072B2", orange = "#E69F00", verm = "#D55E00", green = "#009E73",
         purple = "#CC79A7", sky = "#56B4E9", yellow = "#F0E442", black = "#000000",
         grey = "#999999")
oi <- function(k) unname(.OI[[k]])

a3_apply_data_ink <- function() {
  ggplot2::update_geom_defaults("point", list(size = A3_POINT))
  ggplot2::update_geom_defaults("line", list(linewidth = A3_LINE))
  ggplot2::update_geom_defaults("smooth", list(linewidth = A3_LINE))
  ggplot2::update_geom_defaults("errorbar", list(linewidth = A3_ERRORBAR_LW,
                                                width = A3_ERRORBAR_W))
  invisible(TRUE)
}
a3_apply_data_ink()

# Legends must sit OUTSIDE the coordinate spine / plot box — never over data
# marks and never clipped by the axis frame. Prefer the margin between panels
# ("axis between the completions"). Numeric inside-panel legend.position = c(x, y)
# is forbidden on official A3 figures.
a3_legend_outside <- function(where = "top") {
  horiz <- where %in% c("top", "bottom")
  theme(
    # ggplot2 4.x: "top"/"bottom"/"left"/"right" are the only outside slots.
    # legend.location = "plot" keeps the key in the plot margin, not the panel.
    legend.position = where,
    legend.location = "plot",
    legend.direction = if (horiz) "horizontal" else "vertical",
    legend.justification = "center",
    legend.justification.top = "center",
    legend.justification.bottom = "center",
    legend.box.just = "center",
    legend.box = if (horiz) "horizontal" else "vertical",
    legend.box.spacing = unit(4, "pt"),
    legend.margin = margin(1, 3, 1, 3),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA)
  )
}

# Extra y expansion so geom_text value labels sit clear of the spine / panel edge.
a3_y_for_labels <- function(bottom = 0.06, top = 0.20, ...) {
  scale_y_continuous(expand = expansion(mult = c(bottom, top)), ...)
}

# Thousands separators for count axes / bar labels. Presentation only.
a3_fmt_n <- function(x) {
  format(as.integer(round(as.numeric(x))), big.mark = ",", scientific = FALSE, trim = TRUE)
}

# Claim-first panel title. Face stays plain: only A--F tags may be bold.
# Renderer titles are always explicit: axis labels are never promoted
# automatically, because on a 2x3 page that duplicates scientific text and
# can cross a neighbouring cell boundary.
a3_claim <- function(title) {
  labs(title = title)
}

theme_a3 <- function() {
  theme_classic(base_size = A3_AXIS_PT, base_family = A3_FONT) +
    theme(
      axis.title = element_text(size = A3_AXIS_PT, face = "plain"),
      axis.text = element_text(size = A3_TICK_PT, face = "plain", colour = "black"),
      axis.line = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = A3_TICK_PT, face = "plain"),
      legend.key.size = unit(0.8, "lines"),
      legend.position = "top",
      legend.location = "plot",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.box.just = "center",
      legend.box.spacing = unit(4, "pt"),
      legend.margin = margin(1, 3, 1, 3),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      plot.title = element_text(size = A3_AXIS_PT, hjust = 0.5, face = "plain",
                                family = A3_FONT,
                                margin = margin(t = 0, b = 2)),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      # Panel letters are NOT a ggplot tag. compose_a3.R places them in a
      # reserved band outside this grob. Leaving plot.tag live re-introduces
      # the topleft crowding the compositor exists to prevent.
      plot.tag = element_blank(),
      strip.text = element_text(size = A3_TICK_PT, face = "plain", family = A3_FONT),
      strip.background = element_blank(),
      panel.grid.major = element_line(linewidth = 0.22, colour = "grey88"),
      # Tight furniture: the compositor band holds the letter, so the ggplot
      # itself should spend its cell on the coordinate region.
      plot.margin = margin(2, 6, 3, 4),
      axis.title.y = element_text(margin = margin(r = 4)),
      axis.title.x = element_text(margin = margin(t = 2))
    )
}

a3_draw <- function(plot) {
  if (inherits(plot, "gtable") || inherits(plot, "grob")) {
    grid::grid.newpage()
    grid::grid.draw(plot)
  } else {
    print(plot)
  }
}

# Save a ggplot/patchwork/gtable as both a 300-dpi PNG (ragg) and a vector PDF
# (cairo_pdf, fonts embedded) into the figures/ directory.
a3_save <- function(plot, stem, w, h, figdir) {
  if (!dir.exists(figdir)) {
    stop(sprintf("Explicit A3 output directory does not exist: %s", figdir), call. = FALSE)
  }
  # Render at the printed text width (6.5in = \linewidth) so text prints at
  # authored point size. Scale the authored height by 6.5/w to keep aspect.
  h <- h * 6.5 / w
  w <- 6.5
  png <- file.path(figdir, paste0(stem, ".png"))
  pdf <- file.path(figdir, paste0(stem, ".pdf"))
  ragg::agg_png(png, width = w, height = h, units = "in", res = 300, background = "white")
  a3_draw(plot); invisible(grDevices::dev.off())
  grDevices::cairo_pdf(pdf, width = w, height = h, family = A3_FONT)
  a3_draw(plot); invisible(grDevices::dev.off())
  # Cairo inserts a wall-clock creation date.  Normalize only the fixed-width
  # PDF metadata value so identical inputs yield byte-identical staged PDFs.
  normalizer <- Sys.getenv("A3_PDF_NORMALIZER", unset = "")
  if (normalizer == "" || !file.exists(normalizer) ||
      system2("python3", c(shQuote(normalizer), shQuote(pdf))) != 0L) {
    stop("A3 staged PDF metadata normalization failed", call. = FALSE)
  }
  cat("wrote", basename(png), "and", basename(pdf), "\n")
}
