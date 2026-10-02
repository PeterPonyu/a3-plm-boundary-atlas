#!/usr/bin/env bash
# Fail-closed A3 release-asset build.  Canonical figures are never overwritten.
set -euo pipefail

cd "$(dirname "$0")/.."
SOURCE_ROOT="$(pwd -P)"
CAPSULE_ROOT="$(cd .. && pwd -P)"
# Keep the direct figure-only entry point self-contained: the R theme requires
# this normalizer to make Cairo PDF metadata deterministic, just like the
# complete submission builder does for its staged worktree.
export A3_PDF_NORMALIZER="$SOURCE_ROOT/scripts/normalize_a3_pdf_metadata.py"
if [[ $# -ne 1 ]]; then
  echo "usage: $0 /absolute/empty-output-directory" >&2
  exit 64
fi
OUTPUT="$1"
case "$OUTPUT" in
  /*) ;;
  *) echo "output directory must be absolute" >&2; exit 64 ;;
esac

# Resolve before mkdir: a lexical /tmp path can still traverse or symlink into
# the capsule, where a failed render must never be allowed to remove assets.
OUTPUT_LEXICAL="$(realpath -m -- "$OUTPUT")"
OUTPUT_PARENT="$(dirname -- "$OUTPUT_LEXICAL")"
if [[ ! -d "$OUTPUT_PARENT" ]]; then
  echo "output parent does not exist: $OUTPUT_PARENT" >&2
  exit 64
fi
OUTPUT_PARENT_REAL="$(realpath -e -- "$OUTPUT_PARENT")"
OUTPUT_REAL="$OUTPUT_PARENT_REAL/$(basename -- "$OUTPUT_LEXICAL")"
case "$OUTPUT_REAL" in
  "$CAPSULE_ROOT"|"$CAPSULE_ROOT"/*)
    echo "refusing output inside source/canonical/package tree: $OUTPUT_REAL" >&2
    exit 64
    ;;
esac
if [[ -e "$OUTPUT_REAL" || -L "$OUTPUT_REAL" ]]; then
  echo "refusing to overwrite output path: $OUTPUT_REAL" >&2
  exit 64
fi
mkdir -- "$OUTPUT_REAL"
CANONICAL_HASHES=""
CURRENT_HASHES=""
cleanup_failed_release() {
  [[ -z "$CANONICAL_HASHES" ]] || rm -f -- "$CANONICAL_HASHES"
  [[ -z "$CURRENT_HASHES" ]] || rm -f -- "$CURRENT_HASHES"
  rm -rf -- "$OUTPUT_REAL"
}
on_signal() {
  cleanup_failed_release
  exit "$1"
}
trap cleanup_failed_release ERR
trap 'on_signal 130' INT
trap 'on_signal 143' TERM

CANONICAL_HASHES="$(mktemp)"
find figures -maxdepth 1 -type f \( -name '*.pdf' -o -name '*.png' \) -print0 \
  | sort -z | xargs -0 sha256sum > "$CANONICAL_HASHES"

A3_FIGURE_OUTPUT="$OUTPUT_REAL" Rscript scripts/R/make_four_panel_all.R
A3_FIGURE_OUTPUT="$OUTPUT_REAL" Rscript scripts/R/make_a3_graphical_abstract.R
python3 scripts/verify_figure_release.py --output "$OUTPUT_REAL"

CURRENT_HASHES="$(mktemp)"
find figures -maxdepth 1 -type f \( -name '*.pdf' -o -name '*.png' \) -print0 \
  | sort -z | xargs -0 sha256sum > "$CURRENT_HASHES"
if ! cmp -s "$CANONICAL_HASHES" "$CURRENT_HASHES"; then
  rm -f -- "$CURRENT_HASHES"
  echo "canonical figure hashes changed during staging build" >&2
  exit 1
fi
rm -f -- "$CURRENT_HASHES" "$CANONICAL_HASHES"
trap - ERR INT TERM
echo "A3 release assets staged at $OUTPUT_REAL (canonical figures unchanged)."
