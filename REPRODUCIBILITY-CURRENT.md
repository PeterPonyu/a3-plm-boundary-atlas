# A3 current author-review revision — 23 September 2026

Status: `AUTHOR_REVIEW_SAME_SCORE_INTEGRATED_NOT_SUBMISSION_READY`.
The clinical score-version gap is closed for this 4,677-variant panel.
Independent scientific review and visual acceptance remain open.

The current working paper is `manuscript/manuscript.tex` and its matching
`manuscript/manuscript.pdf`; `manuscript/SI.tex` and `manuscript/SI.pdf` are the
matching supplement. `manuscript/release-current.json` binds their actual
bytes. The older `submission/flat_Proteins` payload is historical, not the
current upload authority. No external deposit or submission has been changed.

The original near-ready paper has been revised in place. Its DMS, architecture,
stability, binding and developability content is retained. Fig. 2 now uses
recovered row scores and newly generated class-stratified AUROC intervals.
Fig. 3 uses five-fold OOF predictions from the same table, with fit-fold-only
standardization, pooled Platt, tie-aggregated isotonic, pooled- and
gene-prevalence baselines, and 1,000 identity-grouped refitting draws.
Fig. 7 and the graphical abstract read the same macro-AUROC and gene interval.
The old score table (mean AUROC 0.6560), six-gene reliability illustration,
degradation/Cox/Holm claims and old intervals are not current evidence.
The gnomAD comparator and S9 predictor exploration remain explicitly
historical summary-level analyses; gnomAD variants were not newly rescored.

The predictor section now distinguishes missing ClinVar score-spread data
from a failed transfer test, and names the actual DMS OLS covariates.

## Rebuild and check

From the capsule root, `bash manuscript/scripts/rebuild_current.sh` checks
locked inputs, independently reconciles clinical outputs, deterministically
exports current plotdata/TeX values, renders figures and compiles both
documents without network access or model inference. It then binds the
actual source/PDF versions and runs the integration checks.

After an intentional source edit, review the scientific change first, then
run `python analysis/revision_20260923/bind_author_review.py --lock-inputs`.
This accepts input hashes only and does **not** attest existing PDFs.
Run the complete rebuild after that. Ordinary rebuilds refuse input drift;
they do not quietly update the method contract. The old submission bundle
is never written by the current build.

Binding hashes only attest version identity, never scientific clearance.
Behavior tests for scoring identities, isotonic ties and grouped bootstrap
are in the revision directories. No dependencies were installed.

## Completed scoring and analysis

`analysis/revision_20260922/fetch_official_checkpoint.py` resumes separately
validated HTTP ranges from the official fair-esm source. It never overwrites
the truncated torch cache. The expected checkpoint length is 2,604,537,549
bytes. All 156 parts are complete with no remaining partial file; the
assembled checkpoint passed all ZIP-member CRCs, hash binding and loading.
SHA-256: `ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c`.
This is a locally computed integrity hash, not an independently published
upstream signature. Official HTTP headers were checked again on 23 September
2026: status 200, matching length and range support.

Despite its historical filename, `rescore_clinvar_candidate.py` has now
actually produced `analysis/revision_20260923/row-scores/scores.csv`.
Its 4,677 variant rows use the validated 20 reference sequences and source
labels. All 24 N6c full-sequence/window gates match to four decimal places;
repeat full-sequence forwards were bitwise identical for all 20 genes.
`row-scores/provenance.json` and `run-binding.json` retain the original run
record (the candidate/HOLD status there describes the stage at scoring time,
not the later manuscript integration). The immutable download record's
`not yet inferred` field likewise describes download completion only.

`analysis/revision_20260923/clinical-unified/` contains the finished analysis:
macro-AUROC 0.7840632, gene-resampling interval [0.7428235, 0.8237764];
Platt/isotonic Brier 0.1411603/0.1407112, versus gene-prevalence 0.1516821.
All figures and new clinical TeX values derive from these outputs. The
1,000 refits condition on gene and class composition and are not external,
gene-disjoint or position-disjoint validation. Original OOF points and refit
means/percentile ranges are separately reported, including ranges that do
not include the original point.

Rerunning inference uses the existing `cofold5090` environment and the
complete local checkpoint. Rerunning clinical analysis requires a fresh
output directory because `analyze_clinical.py` refuses overwrites. No
dependency was installed.

Numerical verification is recorded in `clinical-validation.json` (75 checks).
`integration-validation.json` separately records document binding and
compile checks. Rendered page PNGs and bbox QA are technical checks only:
`view_image` returned `view_image is not allowed because you do not support
image inputs`; therefore **visual review was not performed** and no visual
PASS is claimed.

The pre-edit original manuscript and contracts are additionally backed up in
`analysis/revision_20260923/pre-integration.08EJ4V/`. The immediate
pre-unified-integration manuscript, scripts and contracts are also preserved
in `analysis/revision_20260923/pre-unified-integration.VNp8RS/`.
