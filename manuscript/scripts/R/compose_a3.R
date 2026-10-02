# A3 explicit panel compositor.
# Panel bodies keep their own axis-label gutters. Uppercase letters and centred
# titles occupy a separate header band, drawn after each body. The current
# position-based layouts give the 20-gene and 19-assay forests two row heights.
# Legacy helpers are retained below for archived source compatibility.

suppressMessages({
  library(ggplot2)
  library(grid)
  library(gtable)
})

if (!exists("A3_FONT")) {
  stop("compose_a3.R requires theme_a3.R to be sourced first.", call. = FALSE)
}

a3_atomic_panel <- function(panel) {
  if (!inherits(panel, "ggplot") || inherits(panel, "patchwork")) {
    stop("canonical compositor requires atomic ggplot panels, never nested patchwork",
         call. = FALSE)
  }
  # Decouple from ggplot tag automation. A leftover labs(tag=) must not draw.
  panel + labs(tag = NULL) + theme(plot.tag = element_blank())
}

# ggplotGrob measures text on the current device. The default PDF device
# does not know Arial, so ylab/tick widths would be fallback-font metrics
# and the independent-cell layout would be wrong. Measure on ragg instead.
a3_ggplot_grob <- function(panel) {
  tmp <- tempfile(fileext = ".png")
  ragg::agg_png(tmp, width = 2.1, height = 2.1, units = "in", res = 150,
                background = "white")
  on.exit({
    try(grDevices::dev.off(), silent = TRUE)
    unlink(tmp)
  }, add = TRUE)
  ggplotGrob(a3_atomic_panel(panel))
}

# One cell: letter in a clearance band, then the panel grob underneath.
# The letter is a grid textGrob with a fixed inset, not plot.tag.
a3_labelled_cell <- function(panel, letter) {
  title <- panel$labels$title
  if (is.null(title) || !nzchar(title)) stop(paste("Missing panel title", letter))
  inner <- a3_ggplot_grob(panel + labs(title = NULL))
  header <- grid::grobTree(
    grid::textGrob(letter, x = unit(1,"pt"), y = unit(.55,"npc"),
      just = c("left","center"), gp = grid::gpar(fontfamily=A3_FONT,
      fontsize=A3_PANEL_LABEL_PT, fontface="bold")),
    grid::textGrob(title, x = unit(.5,"npc"), y = unit(.55,"npc"),
      just = "center", gp = grid::gpar(fontfamily=A3_FONT, fontsize=A3_AXIS_PT))
  )
  gt <- gtable::gtable(widths=unit(1,"null"),
      heights=unit.c(unit(A3_LABEL_BAND_PT,"pt"),unit(1,"null")))
  gt <- gtable::gtable_add_grob(gt,inner,t=2,l=1,clip="off")
  gtable::gtable_add_grob(gt,header,t=1,l=1,clip="off")
}

# Explicit cell geometry. The forest consumes two row heights; diagnostic
# plots keep a readable coordinate region without shrinking their typography.
a3_compose_positions <- function(panels, positions, weights) {
  pad <- unit(A3_PAGE_PAD_PT,"pt")
  heights <- unit.c(pad)
  for(i in seq_along(weights)) {
    heights <- unit.c(heights,unit(weights[i],"null"))
    if(i<length(weights)) heights <- unit.c(heights,unit(A3_ROW_GUTTER_PT,"pt"))
  }
  gt <- gtable::gtable(widths=unit.c(pad,unit(1,"null"),unit(A3_CELL_GUTTER_PT,"pt"),unit(1,"null"),pad), heights=unit.c(heights,pad))
  for(i in seq_along(panels)) {
    q <- positions[[i]]
    gt <- gtable::gtable_add_grob(gt,a3_labelled_cell(panels[[i]],LETTERS[i]),
      t=q[1]*2,b=q[2]*2,l=q[3]*2,r=q[4]*2,clip="off")
  }
  gt
}

# 2x3 of independent grobs. Equal cell weights; no shared ylab column.
a3_compose_six <- function(panels) {
  if (!is.list(panels) || length(panels) != 6L) {
    stop("canonical compositor accepts exactly 6 panels", call. = FALSE)
  }
  cells <- Map(a3_labelled_cell, panels, as.list(LETTERS[1:6]))
  # [A][g][B][g][C]
  # [g]
  # [D][g][E][g][F]
  gutter <- unit(A3_CELL_GUTTER_PT, "pt")
  pad <- unit(A3_PAGE_PAD_PT, "pt")
  # Outer pad keeps the letter band off the page edge; gutters keep cells
  # from sharing a ylab column (the patchwork over-alignment failure mode).
  gt <- gtable::gtable(
    widths = unit.c(pad, unit(1, "null"), gutter, unit(1, "null"), gutter, unit(1, "null"), pad),
    heights = unit.c(pad, unit(1, "null"), unit(A3_ROW_GUTTER_PT, "pt"), unit(1, "null"), pad),
    name = "a3-six"
  )
  cols <- c(2L, 4L, 6L, 2L, 4L, 6L)
  rows <- c(2L, 2L, 2L, 4L, 4L, 4L)
  for (i in seq_along(cells)) {
    gt <- gtable::gtable_add_grob(
      gt, cells[[i]], t = rows[[i]], l = cols[[i]],
      clip = "off", name = paste0("a3-panel-", LETTERS[[i]])
    )
  }
  gt
}

# Same canvas and type scale as a3_compose_six. Panel A takes the top-left
# 2x2 so a long identifier forest has both stacked-glyph height and tick
# width. B and C stack on the remaining top-right third; D--F keep the
# bottom row. Used only where an equal 3x2 cell is geometrically too small.
a3_compose_tall_a <- function(panels) {
  if (!is.list(panels) || length(panels) != 6L) {
    stop("tall-A compositor accepts exactly 6 panels", call. = FALSE)
  }
  cells <- Map(a3_labelled_cell, panels, as.list(LETTERS[1:6]))
  gutter <- unit(A3_CELL_GUTTER_PT, "pt")
  pad <- unit(A3_PAGE_PAD_PT, "pt")
  gt <- gtable::gtable(
    widths = unit.c(pad, unit(1, "null"), gutter, unit(1, "null"), gutter, unit(1, "null"), pad),
    heights = unit.c(pad, unit(1, "null"), unit(A3_ROW_GUTTER_PT, "pt"),
                     unit(1, "null"), unit(A3_ROW_GUTTER_PT, "pt"), unit(1, "null"), pad),
    name = "a3-tall-a"
  )
  # [ A spanning 2x2 ] [B]
  # [ A spanning 2x2 ] [C]
  # [D] [E] [F]
  places <- list(
    list(t = 2L, l = 2L, b = 4L, r = 4L),
    list(t = 2L, l = 6L, b = 2L, r = 6L),
    list(t = 4L, l = 6L, b = 4L, r = 6L),
    list(t = 6L, l = 2L, b = 6L, r = 2L),
    list(t = 6L, l = 4L, b = 6L, r = 4L),
    list(t = 6L, l = 6L, b = 6L, r = 6L)
  )
  for (i in seq_along(cells)) {
    spec <- places[[i]]
    gt <- gtable::gtable_add_grob(
      gt, cells[[i]], t = spec$t, l = spec$l, b = spec$b, r = spec$r,
      clip = "off", name = paste0("a3-panel-", LETTERS[[i]])
    )
  }
  gt
}

# Reader-oriented layouts: two columns at final width; no requirement to fill
# a panel with unavailable statistics. A three/five-panel figure gives A the
# full top row; every data panel retains the same type size.
a3_compose_readable <- function(panels, tall_first = FALSE) {
  n <- length(panels)
  stopifnot(n >= 2L, n <= 6L)
  cells <- Map(a3_labelled_cell, panels, as.list(LETTERS[seq_len(n)]))
  spanning <- n %% 2L == 1L
  nrows <- ceiling(n / 2)
  weights <- rep(1, nrows)
  if (tall_first) weights[1] <- 1.8
  pad <- unit(A3_PAGE_PAD_PT, "pt")
  heights <- unit.c(pad)
  for (r in seq_len(nrows)) {
    heights <- unit.c(heights, unit(weights[r], "null"))
    if (r < nrows) heights <- unit.c(heights, unit(A3_ROW_GUTTER_PT, "pt"))
  }
  heights <- unit.c(heights,pad)
  gt <- gtable::gtable(widths=unit.c(pad,unit(1,"null"),unit(A3_CELL_GUTTER_PT,"pt"),unit(1,"null"),pad), heights=heights)
  for (i in seq_len(n)) {
    if (spanning && i==1L) {r<-1L;c<-2L;right<-4L} else {
      k<-if(spanning) i-2L else i-1L
      r<-k%/%2L+if(spanning) 2L else 1L
      c<-if(k%%2L==0L) 2L else 4L;right<-c
    }
    gt<-gtable::gtable_add_grob(gt,cells[[i]],t=2L*r,l=c,r=right,clip="off",name=paste0("reader-panel-",LETTERS[i]))
  }
  gt
}
