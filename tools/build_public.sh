#!/usr/bin/env bash
# Explicit public rebuild from saved inputs; no model inference or scientific fits.
# Run only when PDF regeneration is intended, in a working copy.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 -B tools/verify_public_release.py
python3 -B analysis/revision_20260923/build_clinical_plotdata.py
cd manuscript
export SOURCE_DATE_EPOCH=0 FORCE_SOURCE_DATE=1 TZ=UTC
export A3_FIGURE_OUTPUT="$PWD/figures"
export A3_PDF_NORMALIZER="$PWD/scripts/normalize_a3_pdf_metadata.py"
Rscript scripts/R/make_four_panel_all.R
Rscript scripts/R/make_a3_graphical_abstract.R
python3 -B scripts/make_binding_objects.py
python3 -B scripts/verify_figure_release.py --output "$PWD/figures"
latexmk -g -pdf -interaction=nonstopmode -halt-on-error manuscript.tex
latexmk -g -pdf -interaction=nonstopmode -halt-on-error SI.tex
cd ..
# Refresh refuses scientific/result/build-code drift. Compiler logs remain temporary.
python3 -B tools/update_public_manifest.py
python3 -B tools/verify_public_release.py
