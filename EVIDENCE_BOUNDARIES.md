# Current evidence, historical snapshots and scientific HOLD

This file is the public export's interpretation boundary. Copied old documentation,
filenames and manifests retain their historical wording as provenance; that wording
must not be promoted over the current source manuscript or this explicit boundary.

## Current clinical authority

`analysis/revision_20260923/row-scores/scores.csv` has 4,677 unique identities in
20 genes and SHA-256 `0fac1a086b5b4797b336484f8f0a6ee6eb00334b83607c353b4ea4cd5788eccd`.
`clinical-unified/` is the matching analysis. Current figure inputs are
`F3_clinvar_20260923.csv` and `clinical_*_20260923.csv`; current TeX clinical
values/table are generated from that analysis. All 24 recorded N6c rounded
AUROC gates are retained. Recorded repeat-forward diagnostics concern the
completed run, not inference performed during public export.

The immutable row-score record's candidate/HOLD and checkpoint download record's
`not yet inferred` describe their creation stage. They are not silently rewritten.
Later integration does not clear scientific HOLD.

OOF uses five within-panel folds, fit-only standardization, unpenalized logistic
MLE and tie-aggregated isotonic PAV. Saved refit records condition on gene/class
composition, keep identity multiplicities and are not gene/position-disjoint or
external clinical validation. An OOF point outside a refit percentile interval
remains visible; it is not adjusted to fit the interval.

## Historical clinical and exploratory records

- `analysis/xfit_calibration_20260918/inputs/manifest_*.csv`: source identities,
  labels and **old scores**, not the current row-score values. Retained because
  the actual analysis/fold tests depend on these identities and old methods.
- `analysis/xfit_calibration_20260918/outputs/`: old calibration, rounds 2/3 and
  repaired identity-grouped candidate analyses. Not current clinical estimates.
  The original withdrawn nested-interval notice is retained explicitly.
- `results-snapshot/N6_clinvar/clinical_calibration_full20.json`, `F7_*` clinical
  plotdata, `F3_clinvar.csv`, older EF/registry records and old claim ledger:
  historical only. No retired 7/20 degradation, Cox/Holm or six-gene headline is
  restored by their inclusion.
- `FS2_*` and current figure's gnomAD comparator: historical summary-only;
  gnomAD variants were not newly rescored.
- `F10_reliability_predictor.csv`, S9 and the corresponding old artifact-index
  claims: historical exploratory data. Clinical score spread was missing and
  transfer was not evaluated. New scores do not retrospectively perform that test.

## Retained nonclinical result snapshots

`results-snapshot/N4_antibody_ddg/results.json`, `results-snapshot/a3_meta/*.csv`,
the saved SKEMPI 3B table/robustness outputs, and nonclinical manuscript plotdata
are legacy research snapshots consumed by the current paper. No new inference,
scientific fit or confidence-interval study was done during export. Plotting
summaries alone do not reconstruct all upstream analyses or missing coefficients.

Antibody interval identity is unresolved between the N4 archive and the separate
antibody-language-model archive. The N4 archive's positive lower endpoint and
the other archive's slightly negative lower endpoint remain distinct. No common
scorer/resampling identity, pooling, significance or subgroup-difference claim
is inferred. General/non-Ab and Ab-antigen are complementary parts of one SKEMPI
set, not independent experiments; interface partitions are separate definitions.

## Manifests and build authority

All 256 source method-locked files and 237 source build inputs are preserved.
`provenance/export-allowlist.json` records source hashes and roles. The original
method/release manifests bind the source capsule at export time and remain
source-version provenance. Old registries and legacy writer docs are historical;
they are not a public-release inventory or evidence that excluded submission
assets are present. The original private submission build path is not the public
entry point. Original compiled temporary logs are deliberately excluded.

`PUBLIC_RELEASE_MANIFEST.json` and `SHA256SUMS.txt` are the public package's fresh
inventory. After coordinator DOI/DAS and compilation changes, refresh them and
rerun verification. Do not relock original capsule manifests merely to silence
drift. A changed scientific analysis source requires independent review.

**HOLD remains:** independent scientific review, separate visual acceptance,
external clinical validation, unresolved antibody provenance and unperformed
gnomAD row-level / clinical-descriptor transfer analyses. A release may expose
research evidence without pretending these scientific questions are resolved.
