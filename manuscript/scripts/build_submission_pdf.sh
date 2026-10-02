#!/usr/bin/env bash
# Build a complete A3 submission package in an external, atomically published directory.
# Usage: build_submission_pdf.sh /absolute/nonexistent-output-directory
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
SOURCE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd -P)"
CAPSULE_ROOT="$(cd "$SOURCE_ROOT/.." && pwd -P)"

if [[ $# -ne 1 ]]; then
  echo "usage: $0 /absolute/nonexistent-output-directory" >&2
  exit 64
fi
OUTPUT="$1"
case "$OUTPUT" in
  /*) ;;
  *) echo "output directory must be absolute" >&2; exit 64 ;;
esac

# Canonicalize the parent before creating anything.  This rejects lexical traversal
# and symlink routes into the capsule, including sibling submission directories.
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
    echo "refusing output inside A3 capsule tree: $OUTPUT_REAL" >&2
    exit 64
    ;;
esac
if [[ -e "$OUTPUT_REAL" || -L "$OUTPUT_REAL" ]]; then
  echo "refusing to overwrite output path: $OUTPUT_REAL" >&2
  exit 64
fi

STAGING="$(mktemp -d "$OUTPUT_PARENT_REAL/.a3-submission.XXXXXXXX")"
WORK="$STAGING/work"
PACKAGE="$STAGING/package"
cleanup() {
  rm -rf -- "$STAGING"
}
on_signal() {
  cleanup
  exit "$1"
}
# EXIT covers explicit validation failures as well as command failures.  Signal
# handlers call cleanup eagerly, then EXIT makes that cleanup idempotent.
trap cleanup EXIT
trap 'on_signal 130' INT
trap 'on_signal 143' TERM

# Stage only the committed, hash-bound transitive closure.  The builder never
# reads canonical figures or arbitrary siblings from the mutable source tree.
python3 "$SCRIPT_DIR/stage_release_inputs.py" \
  --source-root "$SOURCE_ROOT" --spec "$SOURCE_ROOT/release-spec.json" --stage "$WORK"
mkdir -- "$PACKAGE"
mkdir -- "$PACKAGE/attestation"
cd "$WORK"
export SOURCE_DATE_EPOCH=0
export FORCE_SOURCE_DATE=1
export TZ=UTC
export A3_PDF_NORMALIZER="$WORK/scripts/normalize_a3_pdf_metadata.py"
mkdir -p -- figures
if [[ -n "${A3_FROZEN_FIGURE_STAGE:-}" ]]; then
  FROZEN_STAGE="$(realpath -e -- "$A3_FROZEN_FIGURE_STAGE")"
  case "$FROZEN_STAGE" in
    "$CAPSULE_ROOT"|"$CAPSULE_ROOT"/*)
      echo "refusing frozen figure stage inside A3 capsule tree: $FROZEN_STAGE" >&2
      exit 64
      ;;
  esac
  FROZEN_MANIFEST="$WORK/registries/frozen_task6_stage.sha256"
  [[ -f "$FROZEN_MANIFEST" ]] || { echo "missing frozen figure manifest" >&2; exit 1; }
  while read -r expected filename; do
    [[ -z "$expected" || "$expected" == \#* ]] && continue
    [[ "$filename" != */* && "$filename" != *F7_20gene* ]] || { echo "invalid frozen figure manifest entry: $filename" >&2; exit 1; }
    [[ -f "$FROZEN_STAGE/$filename" ]] || { echo "missing frozen figure asset: $filename" >&2; exit 1; }
    actual="$(sha256sum "$FROZEN_STAGE/$filename" | awk '{print $1}')"
    [[ "$actual" == "$expected" ]] || { echo "frozen figure hash mismatch: $filename" >&2; exit 1; }
    cp "$FROZEN_STAGE/$filename" "$WORK/figures/$filename"
  done < "$FROZEN_MANIFEST"
  echo "using hash-checked external frozen figure overlay: $FROZEN_STAGE"
else
  A3_FIGURE_OUTPUT="$WORK/figures" Rscript scripts/R/make_four_panel_all.R
fi
A3_FIGURE_OUTPUT="$WORK/figures" Rscript scripts/R/make_a3_graphical_abstract.R
python3 scripts/verify_figure_release.py --output "$WORK/figures"
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error SI.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error SI.tex >/dev/null

ACTUAL_PAGES="$(pdfinfo manuscript.pdf | awk -F: '/^Pages:/ {gsub(/^[[:space:]]+/, "", $2); print $2}')"
ACTUAL_SI_PAGES="$(pdfinfo SI.pdf | awk -F: '/^Pages:/ {gsub(/^[[:space:]]+/, "", $2); print $2}')"
EXPECTED_PAGES="$(python3 -c 'import json; print(json.load(open("release-spec.json", encoding="utf-8"))["documents"]["main"]["pages"])')"
EXPECTED_SI_PAGES="$(python3 -c 'import json; print(json.load(open("release-spec.json", encoding="utf-8"))["documents"]["si"]["pages"])')"
if [[ "$ACTUAL_PAGES" != "$EXPECTED_PAGES" ]]; then
  echo "release manuscript has $ACTUAL_PAGES pages; release specification requires $EXPECTED_PAGES" >&2
  exit 1
fi
if [[ "$ACTUAL_SI_PAGES" != "$EXPECTED_SI_PAGES" ]]; then
  echo "release SI has $ACTUAL_SI_PAGES pages; release specification requires $EXPECTED_SI_PAGES" >&2
  exit 1
fi
if rg -q '^(LaTeX Font Warning|! |Emergency stop|Output loop)' manuscript.log SI.log; then
  echo "release LaTeX log contains a configured fatal diagnostic" >&2
  exit 1
fi
if pdffonts manuscript.pdf | awk 'NR > 2 && ($0 ~ /Type 3/ || $0 ~ /[[:space:]]no[[:space:]]+no[[:space:]]/) {exit 1}'; then
  :
else
  echo "release manuscript contains unembedded or Type 3 fonts" >&2
  exit 1
fi
if pdffonts SI.pdf | awk 'NR > 2 && ($0 ~ /Type 3/ || $0 ~ /[[:space:]]no[[:space:]]+no[[:space:]]/) {exit 1}'; then
  :
else
  echo "release SI contains unembedded or Type 3 fonts" >&2
  exit 1
fi

cp manuscript.pdf SI.pdf "$PACKAGE/"
cp manuscript.log "$PACKAGE/attestation/manuscript.log"
cp SI.log "$PACKAGE/attestation/SI.log"
cp "$SOURCE_ROOT/attestation/log-policy.json" "$PACKAGE/attestation/log-policy.json"
cp figures/A3_graphical_abstract.pdf figures/A3_graphical_abstract.png "$PACKAGE/"
cp submission/flat_Proteins/cover_letter_Proteins.md "$PACKAGE/"
cp submission/flat_Proteins/submission_checklist.md "$PACKAGE/"
printf '%s\n' \
  'A3 externally staged candidate submission package' \
  'Main article PDF carries seven figures. Supporting information is a separate SI.pdf (Figures S1--S9).' \
  'The post-compile visual-registry.json is the only page-location authority.' \
  'Candidate-only F7_20gene assets are excluded.' > "$PACKAGE/MANIFEST.txt"
python3 scripts/verify_submission_package.py \
  --package "$PACKAGE" \
  --source-spec release-spec.json \
  --source-registry visual-review-source.json \
  --builder "$WORK/scripts/build_submission_pdf.sh"
(
  cd "$PACKAGE"
  printf '%s\n' manuscript.pdf SI.pdf A3_graphical_abstract.pdf A3_graphical_abstract.png \
    cover_letter_Proteins.md submission_checklist.md MANIFEST.txt visual-registry.json \
    | LC_ALL=C sort | xargs sha256sum > SHA256SUMS.txt
)
(
  cd "$PACKAGE"
  sha256sum -c SHA256SUMS.txt >/dev/null
)

# Narrow test seam for the release-contract regression: it exercises a failure
# after all staged validation, before publication, without touching canonical
# sources or an output target.
if [[ "${A3_RELEASE_TEST_FAIL_AFTER_STAGING:-0}" == "1" ]]; then
  echo "intentional post-staging validation failure" >&2
  exit 1
fi

# PACKAGE and OUTPUT_REAL share an already-canonicalized parent filesystem, so mv is
# the publish operation; a failure before this point leaves no output target behind.
mv -- "$PACKAGE" "$OUTPUT_REAL"
trap - ERR INT TERM
rm -rf -- "$STAGING"
echo "A3 complete ${ACTUAL_PAGES}-page manuscript + ${ACTUAL_SI_PAGES}-page SI submission package published at $OUTPUT_REAL"
