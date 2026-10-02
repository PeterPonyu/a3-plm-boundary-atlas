#!/usr/bin/env Rscript
# GA-only entry point retained for the existing manuscript build. The native
# vector renderer reads current saved evidence; it never infers model scores,
# re-fits a map, changes another figure, or rebinds the manuscript contract.
args <- grep('^--file=', commandArgs(FALSE), value = TRUE)
if (length(args) != 1L) stop('Expected one graphical-abstract entry point')
HERE <- dirname(normalizePath(sub('^--file=', '', args[[1]])))
FIG <- Sys.getenv('A3_FIGURE_OUTPUT', unset = '')
if (!nzchar(FIG) || !dir.exists(FIG)) {
  stop('Explicit A3 graphical-abstract output directory is required')
}
helper <- normalizePath(file.path(HERE, '..', 'make_a3_graphical_abstract_vector.py'))
status <- system2('python3', c('-B', shQuote(helper), '--output', shQuote(normalizePath(FIG))))
if (status != 0L) stop('A3 graphical-abstract vector export failed')
