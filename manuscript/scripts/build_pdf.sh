#!/usr/bin/env bash
# Compatibility entry point for the complete, fail-closed A3 submission build.
#
# SOURCE OF TRUTH = manuscript.tex (hand-maintained LaTeX).  The release command
# regenerates figures and compiles a complete manuscript only in external staging;
# it never writes canonical figures, PDFs, logs, or submission assets in this tree.
set -euo pipefail
exec bash "$(dirname "$0")/build_submission_pdf.sh" "$@"
