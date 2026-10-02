# A3 original-manuscript integration handoff — 23 September 2026

Status: `AUTHOR_REVIEW_SAME_SCORE_INTEGRATED_NOT_SUBMISSION_READY`.

The score-recovery dependency is resolved. Actual new inference, same-score
clinical analysis and original manuscript/SI integration are complete for the
4,677-variant, 20-gene panel. This is not merely a candidate plan or another
review. No submission-readiness or independent-review approval is claimed.

## Checkpoint and inference

- Official checkpoint: `analysis/revision_20260922/candidate-v1-checkpoint/esm2_t33_650M_UR50D.pt`.
- 156/156 complete HTTP-range parts, exactly 2,604,537,549 bytes, zero bytes
  remaining and no `.partial` files. Download/inference processes finished
  normally; none remains waiting on the network.
- SHA-256: `ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c`.
  This locally computed integrity hash is not represented as an upstream
  published signature. Download hash binding, all ZIP member CRCs and actual
  fp32 model loading passed.
- Actual 651,043,254-parameter inference used the existing `cofold5090`
  environment, CUDA, fp32, TF32 disabled. No dependency was installed.
- `row-scores/scores.csv`: 4,677 unique variants, 20 genes, source-matched
  labels and reference-sequence hashes. Whole and mutation-centred-window
  LLRs are retained with exact window bounds. Raw scores were not inferred
  from AUROC or mixed with the old calibration table.
- All 24 full-sequence/window AUROCs agree with N6c at recorded four-decimal
  precision. Repeat full-sequence forwards are bitwise identical for all
  20 genes. Per-window repetition is not claimed.
- The immutable scoring provenance still says candidate/HOLD because that
  was the stage at inference time; the download record's `not yet inferred`
  field likewise describes download completion. The current release and
  this handoff record describe the later integration.

Official source checking used the primary fair-esm model endpoint. The web
tool returned no usable text from GitHub; a direct read-only HEAD request to
the official checkpoint endpoint returned HTTP 200, content-length
2604537549, `accept-ranges: bytes`, last-modified 18 August 2022, and ETag
`12a18098227c0ff911354647d25d494d-311`. This network limitation did not block
completed downloading or inference.

## Same-score results actually computed

`clinical-unified/summary.json` and its CSVs are the current clinical evidence.

- Mean per-gene AUROC: 0.7840631995; 10,000-draw gene-resampling percentile
  interval [0.7428235431, 0.8237763521]. Per-gene intervals use 2,000
  class-stratified variant draws. AICDA's empirical interval is [1,1],
  reflecting perfect observed separation, not zero population uncertainty.
- Five-fold OOF analysis uses seed 20260923, fitting-fold-only standardization,
  pooled unpenalized logistic MLE and distinct-score weighted PAV isotonic
  regression with step predictions. Ties do not depend on row order.
- Platt: ECE 0.0211221, Brier 0.1411603, log loss 0.4440988.
- Isotonic: ECE 0.0083165, Brier 0.1407112, log loss 0.4445450.
- Gene-prevalence baseline: Brier 0.1516821, log loss 0.4657190.
  Pooled prevalence has ECE approximately 7.4e-8 but Brier 0.1807806;
  the paper explicitly rejects ECE-only claims of usefulness.
- Paired Brier differences from the gene baseline: Platt -0.0105218,
  refit range [-0.0160730, -0.0072239]; isotonic -0.0109709,
  [-0.0159787, -0.0070170].
- 1,000 refitting draws resample within gene and class, retain bootstrap
  multiplicity as weights and keep every original variant in one fold.
  These remain within-panel, not gene/position/homology/external validation.
- OOF points and refit means/percentiles are separately reported. Isotonic
  ECE point 0.0083165 falls below the refit interval [0.0090162, 0.0253387];
  it was not shifted or hidden. Split-half endpoints are exported.

## Original manuscript files changed

All paths below are capsule-relative; no shared docs or other paper was written.

- `manuscript/manuscript.tex`, `manuscript/text/clinical_revision.tex`,
  `manuscript/text/abstract.tex`, `manuscript/text/methods.tex`:
  actual recovered discrimination/calibration, precise window protocol,
  bootstrap scope and proper-score baselines integrated into the original
  near-ready paper. Nonclinical evidence retained. The missing clinical
  score-spread comparison remains “not evaluated”, not “failed transfer”.
- `manuscript/SI.tex`: same-score Tables S1/S2 added; historical S2 gnomAD
  and S9 predictor boundaries explicit.
- `manuscript/captions/fig2.tex`, `fig3.tex`, `fig8.tex`, `figs2.tex`:
  same-score figure provenance and interval meaning corrected.
- `manuscript/text/clinical_values_20260923.tex`,
  `clinical_table_20260923.tex`: generated from the analysis, not independently
  typed headline numbers.
- `manuscript/figures/plotdata/F3_clinvar_20260923.csv` and
  `clinical_*_20260923.csv`: new version-named exports; old plotdata retained.
- `manuscript/scripts/R/make_four_panel_all.R`,
  `make_a3_graphical_abstract.R`: Fig. 2, Fig. 3, Fig. 7 and the graphical
  abstract consume the new tables. Old six-gene calibration/coverage and
  unmatched model-comparison inputs are absent from the current Fig. 3.
- `manuscript/scripts/rebuild_current.sh`, `verify_figure_release.py`,
  `manuscript/release-current.json`, `release-spec.json`,
  `reproducibility/method-contract.json`, `REPRODUCIBILITY-CURRENT.md`,
  `reproducibility/METHOD-AUDIT.md`, `manuscript/FIGURE-PROVENANCE.md`:
  current author-review authority and fail-closed rebuild chain.
- New analysis/export/validation/test helpers are in this revision directory.
  Scoring/download improvements remain in `analysis/revision_20260922/`.

The unsupported old per-gene intervals, six-gene reliability display and
old-score degradation/Cox/Holm/clinical-utility conclusions are not current
findings. gnomAD remains summary-only; it was not newly rescored. No claim
that the historical predictor analysis is retroactively completed is made.

## Verification actually executed

- `current-tests.xml` / `current-tests.log`: **42 targeted tests passed**.
  These cover scoring identity/window behavior, incomplete-checkpoint
  rejection, PAV ties and agreement with independent sklearn fitted values,
  fit-only scaling, multiplicity weights, old-copy leakage regression,
  new current-release contracts, shipped figure labels and legibility.
  This is not a claim that the entire legacy submission test suite passed.
- `clinical-validation.json`: **75 numerical checks passed**, including
  independent sklearn AUROC/proper scores, OOF reconstruction, label/fold
  identity, reliability bins and paired percentile calculations.
- `integration-validation.json`: **17 checks passed**, no source/PDF drift,
  no Overfull or undefined-reference/citation warnings.
- Two complete invocations of `manuscript/scripts/rebuild_current.sh`
  succeeded. `rebuild-current-run.log` and `rebuild-current-repeat.log`
  preserve their actual outputs. `rebuild-repeat-validation.json` confirms
  identical article/SI hashes and identical 36 PDF/PNG inventory files
  (17 rebuilt current assets plus retained historical files in that inventory).
- `document-qa.json`: all **27 article + 10 SI pages rendered** at 100 dpi;
  zero extracted text boxes outside page bounds, no blank extracted-text
  pages, all fonts embedded, no Type 3. Figure checker separately passed
  all 17 current assets. These automated checks are not visual inspection.
- `visual_review.status = NOT_PERFORMED`, `pass_attested = false`.
  `functions.view_image` returned exactly: `view_image is not allowed because
  you do not support image inputs`. No visual PASS is claimed.
- **29 existing historical submission files** match the pre-integration
  backup. The old package was not rebuilt or repacked. The earlier
  `HISTORICAL-NOT-CURRENT.md` marker predates the renewed no-package-write
  instruction and was not modified in this continuation.

Current PDFs:

- `manuscript/manuscript.pdf`, 27 pages,
  SHA-256 `02bbd817f10dfc58bce3d9e6605eaac609c49ba6fbef19cf169d668bc9c36115`.
- `manuscript/SI.pdf`, 10 pages,
  SHA-256 `688bee039d041f167fd54627d3a5936d7f435a0451e23da66f19bd194ee213a2`.

## Remaining genuine blockers and limits

1. Visual acceptance is not performed; an image-enabled or human reviewer
   must inspect the rendered documents. Automated bbox/font checks passed
   but cannot detect all layout defects.
2. Independent scientific review of the recovered scoring and the revised
   analysis is still required. In particular, the bootstrap conditions on
   gene/class composition and can be shifted relative to original OOF
   estimates. No statistical interpretation was widened to external use.
3. gnomAD row-level rescoring and clinical descriptor-transfer analysis are
   outside this completed panel chain and remain explicitly unperformed.
   They must not be promoted from historical summaries without new evidence.
4. No new submission package, external deposit or journal upload was made.
   The current draft explicitly says not submission-ready. Historical
   packaging scripts/tests have older release semantics; the tested current
   entrypoint is `scripts/rebuild_current.sh`.

There is **no remaining checkpoint/network/inference blocker** for this panel.
Further inference should not be restarted just because old run records use
candidate-era status labels.

## Final handoff recheck — 23 September 2026, 02:43 PDT

This final check did **not** rerun model inference, the 1,000 clinical refits,
or the manuscript render/build. It verified the completed outputs and their
current contracts, leaving scores, estimates, manuscript source and PDFs
unchanged.

- `reproducibility/verify_method_contract.py`: all **254** locked inputs/code
  files match. `release-current.json` binds **235** manuscript build inputs
  and the two PDFs listed above. `release-spec.json` is identical to the
  current release manifest, not an instruction to repack the old submission.
- `handoff-contract-tests.xml`: **17 tests passed, 0 skipped**, comprising
  current-contract, shipped visual-geometry and legibility tests. These
  overlap the previously recorded 42 tests and must not be added as 59
  independent tests. The prior 42-test record has 0 failures/errors/skips.
- All saved analysis output hashes still match `clinical-validation.json`,
  whose 75 numerical checks passed. The 17 saved integration checks and the
  two-build repeatability record remain bound to the unchanged PDFs.
- The 17-asset figure contract and document bbox/font QA were rerun and
  passed. The 37 retained page PNGs are present and their hashes are recorded.
  `visual_review.status` remains **NOT_PERFORMED**, not a visual PASS.
- The same 29 pre-existing historical submission files still match the
  backup. No submission file was written during this handoff continuation.

Frozen identifiers at handoff:

- Scores CSV SHA-256:
  `0fac1a086b5b4797b336484f8f0a6ee6eb00334b83607c353b4ea4cd5788eccd`.
- `manuscript/release-current.json` SHA-256:
  `5d77bb51c61df53a3060a72b5c62fbbd2771ce1bf724679225014d5e1bab2796`.
- `reproducibility/method-contract.json` SHA-256:
  `33bc4d6dcda97e2b3ad130be483da10e01f2f2353eabbeeee5106c2e013b85af`.

The sole current build entrypoint is:

```bash
bash /home/zeyufu/Desktop/labs/capsules/a3-plm-boundary-atlas/manuscript/scripts/rebuild_current.sh
```

It checks the locked method inputs before writing, regenerates only clinical
plotdata/TeX summaries and manuscript assets from existing evidence, compiles
the original article and SI, then binds and checks the resulting source/PDF
pair. It does not call the PLM scorer or the 1,000-draw analysis and never
writes a submission package. An intentional source change needs explicit
input-lock acceptance before a subsequent rebuild; normal rebuilds do not
silently accept changed scientific inputs.

For the parent: A3 is ready for **integration handoff**, not for submission.
The active acceptance holds are visual inspection and independent scientific
review. Lack of new gnomAD row scores, an unperformed descriptor-transfer
analysis and absence of external clinical validation remain stated scope
limits, not newly asserted positive or negative results. Shared cross-paper
reports and B1 remain with the parent; this agent wrote neither.

Backups: `pre-integration.08EJ4V/` and
`pre-unified-integration.VNp8RS/` preserve prior manuscript/scripts/contracts.
All inherited dirty-tree edits were preserved. No git reset, dependency
installation, other-paper write or shared-doc write was performed.
