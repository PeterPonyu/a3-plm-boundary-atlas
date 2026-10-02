# Current method audit — 23 September 2026

`AUTHOR_REVIEW_SAME_SCORE_INTEGRATED_NOT_SUBMISSION_READY`.

The score-version mismatch is now resolved by actual fp32 inference on all
4,677 variants from the complete official checkpoint. Twenty-four
full-sequence/window AUROC gates reproduce the corrected N6c summaries to
their recorded precision. Labels and reference identities are unchanged;
scores were never reconstructed from summary AUROCs. The run table is bound
to checkpoint, source/sequence and script hashes.

The active original manuscript, Fig. 2/3/7, graphical abstract and SI Tables
S1/S2 now use `A3-CLINICAL-UNIFIED-20260923`. Clinical AUROC intervals are
newly generated, not the old unsupported bounds. Calibration is five-fold,
class-stratified with seed 20260923, fit-only scaling, unpenalized logistic
MLE, distinct-score weighted PAV with step interpolation, and fitting-fold
prevalence baselines. No per-gene significance family is carried forward.

Refitting uses 1,000 draws within gene and class, one entry per original
variant and multiplicity weights. Identity copies cannot cross folds; all
maps and baselines are refitted in every draw. This fixes copy leakage, not
within-position, gene, homology or population dependence. The bootstrap
distribution can be shifted from the original OOF estimate; means, percentile
ranges and split-half endpoints are explicitly retained. Pooled-prevalence
ECE is near zero by design despite worse proper scores, so no ECE-only
usefulness claim is allowed. Brier comparisons to the gene-prevalence baseline
are paired by bootstrap draw.

The old 200 identity-grouped refits used the superseded scores. They remain
historical, as do earlier tie-sensitive isotonic outputs, six-gene
calibration, Cox/Holm, degradation and clinical-utility conclusions. The
historical six-gene seed 1/all-score scaling is disclosed, not reused.
gnomAD and S9 predictor evidence are clearly labelled summary-only; a missing
clinical score-spread variable is not evidence of failed transfer.

Current tests check source-label/reference identity, incomplete-weight
rejection, repeatability, tie handling, perfect-AUROC percentile behavior,
weighted bootstrap multiplicity and original-identity fold separation.
The completed inference and analysis pass an independent 75-check
reconciliation. Behavioral tests compare PAV values to scikit-learn,
validate held-out scaling, bootstrap multiplicity and perfect-AUROC
degeneracy. These do not substitute for independent clinical validation.
Document hashes attest source/PDF identity only. Rendering and bbox/font
checks are separate from visual inspection, which was not possible in the
agent's image-disabled runtime. No scientific-ready attestation is issued.
