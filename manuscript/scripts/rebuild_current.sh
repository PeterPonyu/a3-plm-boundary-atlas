#!/usr/bin/env bash
# Current manuscript build. Results are read from locked inputs; no model inference.
set -euo pipefail
cd "$(dirname "$0")/.."
export SOURCE_DATE_EPOCH=0 FORCE_SOURCE_DATE=1 TZ=UTC
python ../reproducibility/verify_method_contract.py
python ../analysis/revision_20260923/validate_clinical.py
python ../analysis/revision_20260923/build_clinical_plotdata.py
# A deterministic export must leave the locked input set unchanged.
python ../reproducibility/verify_method_contract.py
export A3_FIGURE_OUTPUT="$PWD/figures"
export A3_PDF_NORMALIZER="$PWD/scripts/normalize_a3_pdf_metadata.py"
Rscript scripts/R/make_four_panel_all.R
Rscript scripts/R/make_a3_graphical_abstract.R
python scripts/make_binding_objects.py
python scripts/verify_figure_release.py --output "$PWD/figures"
latexmk -g -pdf -interaction=nonstopmode -halt-on-error manuscript.tex
latexmk -g -pdf -interaction=nonstopmode -halt-on-error SI.tex
python ../analysis/revision_20260923/bind_author_review.py
python ../analysis/revision_20260923/verify_integration.py
