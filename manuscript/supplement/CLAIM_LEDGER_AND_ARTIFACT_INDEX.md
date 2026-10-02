# Canonical claim ledger and artifact index

**Scope:** A3 protein-LM mutation-effect boundary-atlas manuscript and supplement.
**Prepared:** 2026-07-11.
**Purpose:** Give reviewers and future maintainers one canonical map from every load-bearing
claim to the committed result artifact(s), the deterministic local command that regenerates it, and
the claim-boundary notes that keep each axis inside its evidence.

This ledger is intentionally conservative. A3's contribution is a **uniform audit protocol**, not a
new model or a SOTA number: it reports where an ESM-2 zero-shot log-likelihood-ratio (LLR) score is
in-panel / abstention-bounded and where it must be abstained from, across six axes, under one raw-score,
bootstrap-everywhere, applicability-gated discipline. The most load-bearing evidence for the
protocol is that it caught and retracted three of its own earlier over-claims; those retractions are
tracked here as first-class ledger rows, not hidden.

## Reproduction environment boundary

Unlike B1 (whose full supplement needs external `~/n11_zenodo` tables), **every A3 figure and headline
regenerates from committed in-repo artifacts with no network and no GPU.** All model inference
(ESM-2/ProtBERT/ESM-1v/AntiBERTy forward passes) was run upstream and its per-variant outputs are
frozen in the committed result JSON/CSV; the figure pipeline is pure tabular re-computation.

| Tier | Network | GPU | External/local input | What it proves |
|---|---:|---:|---|---|
| A3 committed-artifact audit | no | no | none beyond this repo | All main and supplement figures/tables regenerate from committed `research/results/**` JSON/CSV. |
| A3 figure re-render | no | no | none beyond this repo; R + `newtx` LaTeX locally | `scripts/export_plot_data.py` reshapes committed JSON to tidy CSV; `scripts/R/make_figures.R` plots; `pdflatex` builds. |

Recommended interpreter and build:

```bash
PY=~/miniconda3/bin/python
test -x "$PY" || PY="$(command -v python3)"
```

## Canonical local commands

Run from the repository root unless noted. No command uploads data or touches any sibling workspace.

```bash
# Full A3 build: regenerate figure plotdata from committed result JSON, then compile the PDF.
bash paper/A3-plm-boundary-atlas/scripts/build_pdf.sh   # -> manuscript.pdf (Times New Roman)

# Regenerate only the tidy plotdata CSVs consumed by the R figure scripts.
$PY paper/A3-plm-boundary-atlas/scripts/export_plot_data.py

# Synthesis boundary-atlas figure and external VEP-comparison figure (R).
Rscript research/results/a3_meta/make_a3_atlas.R
Rscript research/results/a3_meta/make_vep_comparison.R
```

## Artifact index

| Artifact | Role | Regeneration / provenance |
|---|---|---|
| `research/results/N1_scale/expand80/scale_{35M,150M,650M,3B}_masked.json` | Scale axis (Fig 1): per-size protein-mean AUROC on the 69-protein de-pseudo-replicated subset. | Upstream committed result; consumed by `export_plot_data.py` -> `make_figures.R`. |
| `research/results/N1_scale/family_650M.json` | Full-217-assay ProteinGym corroboration (650M protein-mean AUROC 0.702). | Upstream committed result. |
| `research/results/N5_arch/N5_summary.json` | Architecture axis (**Fig S6**, demoted ex-F2): ProtBERT vs ESM-2 650M mean AUROC and per-protein Spearman; ESM-1v full-benchmark agreement. | Upstream committed result. |
| `research/results/N6_clinvar/n6c_consistent.json` | Clinical axis (Fig 3): per-gene ClinVar AUROC, one fp32 consistent re-score; supersedes the withdrawn long-protein-artifact picture. | Upstream committed result. |
| `research/results/N6_clinvar/clinical_calibration_full20.json` | Clinical calibration/abstention inventory on the full 20-gene panel (printed face uses raw ECE 0.2443 only; calibrated pooled tokens are NOT A PRINTED CLAIM; abstention; leave-P53-out). | Upstream committed result. |
| `research/results/N6_clinvar/clinical_calibration.json` | Clinical calibration on the six-gene reproducibility subset (ECE 0.157->0.047). | Upstream committed result. |
| `research/results/N10_ddg/ddg_calibration_v2.json` | Stability axis (Fig 5): homology-clustered split calibration-RMSE inflation (+0.04) and Spearman 0.358--0.384; clean cluster-disjoint re-run. | Upstream committed result. |
| `research/results/N8_uq/` | Supervised-embedding UQ under shift (**Fig S7**, demoted ex-F6): FLIP splits, deep-ensemble / MC-dropout / GP / BayesianRidge; acquisition 0/7. | Upstream committed result. |
| `research/results/N4_antibody_ddg/results.json` | Binding axis (Fig 7): SKEMPI v2 Spearman(LLR, -ddG) by structural subset, 1000x bootstrap CIs. | Upstream committed result. |
| `research/results/N4_antibody_ddg/results_3b_wt.json` | Binding scale sensitivity (Fig S3): 650M->3B does not rescue interface awareness (0.197->0.148). | Upstream committed result. |
| `research/results/N9_antibody/developability_corrected.json` | Developability axis (Fig 8): per-property BH-significant fractions; single-parent vs multi-germline applicability split. | Upstream committed result (corrected map). |
| `research/results/N9_antibody/antibody_lm_baseline.json` | AntiBERTy control (**Fig S8**, demoted ex-F9): antibody-specific LM on the same data (aggregation 0/7; antibody binding ~0.04; polyreactivity 3/14). | Upstream committed result. |
| `research/results/A3_reliability_predictor.json` | Reliability-predictor axis (**Fig S9**, demoted ex-F10): candidate covariates vs per-task AUROC across DMS and ClinVar; best single Spearman 0.25 (CI touches 0), multivariate R^2 ~0.03. | Upstream committed result. |
| `research/results/a3_meta/a3_atlas_{auroc,spearman}.csv` | Synthesis boundary atlas (Fig 11): AUROC and Spearman panels; a summary of Figs 1,3--5,7--8 and S9, not new data. | `Rscript research/results/a3_meta/make_a3_atlas.R`. |
| `research/results/a3_meta/vep_comparison.csv` | External VEP comparison (Fig S5): ESM-1v median |rho|=0.49 among published VEPs on the Livesey & Marsh compilation. | `Rscript research/results/a3_meta/make_vep_comparison.R`; source Figshare 38248524. |
| `paper/A3-plm-boundary-atlas/figure_ledger.md` | One-claim-per-figure provenance and regeneration commands (authoritative figure map). | Hand-maintained alongside the figures. |

## Claim ledger

| ID | Manuscript claim | Status | Canonical evidence | Local check |
|---|---|---|---|---|
| A3-CL-01 | DMS-fitness discrimination rises 35M->650M then plateaus: raw-score protein-mean AUROC 0.667/0.710/0.725/0.725 (n=69); no gain detected at 3B (paired-bootstrap CI spans 0); full-217 corroborates at 0.702. Bounded as "no gain detected", not "3B provably adds nothing". | Headline (efficiency ceiling) | `N1_scale/expand80/scale_*_masked.json`; `N1_scale/family_650M.json` | Inspect per-size protein-mean AUROC and the 3B-650M paired-bootstrap interval. |
| A3-CL-02 | The signal is shared across independently trained MLM encoders (ProtBERT 0.711; ESM-1v 0.711 / per-protein Spearman 0.73 on the full benchmark), not ESM-specific. Scoped to MLM encoders; autoregressive/structure-conditioned models untested; ESM-1v was purpose-trained for variant effect. | Headline (encoder-robustness, scoped) | `N5_arch/N5_summary.json` | Inspect ProtBERT/ESM-1v mean AUROC and per-protein Spearman vs ESM-2. |
| A3-CL-03 | Zero-shot LLR discriminates ClinVar pathogenic vs benign per gene (mean AUROC 0.78 over 20 genes; 19/20 individually significant; 16/20 >=0.70) under an explicit base-model-leakage caveat; six-gene calibratable under the historical 50/50 protocol (pooled ECE 0.157->0.047; n_held-out=547) and supports selective abstention; historical full-20 printed ECE is uncalibrated (raw, no Platt) 0.2443 (n_held-out=2339); a NEW 5-fold cross-fit (A3-XFIT-CAL-20260918; not a restoration of 0.0207) gives uncalibrated 0.248 [0.237, 0.259] and Platt 0.019 [0.014, 0.032] (n=4677 OOF) with 7/20 genes degrading; the follow-up A3-XFIT-CAL-R2-20260918 attributes that 7/20 to prevalence mismatch (Spearman +0.815) not small n, Holm-significant degraders OTC and PIK3CA only, and per-gene Platt does not significantly beat pooled; robust to dropping P53. | Headline (clinical, leakage-bounded) | `N6_clinvar/n6c_consistent.json`; `clinical_calibration_full20.json`; `clinical_calibration.json`; `analysis/xfit_calibration_20260918/outputs/xfit_calibration.json`; `analysis/xfit_calibration_20260918/outputs/round2/round2_summary.json` | Inspect per-gene AUROC array, six-gene ECE before/after, printed raw full-20 ECE, new xfit bins/CIs, r2 per-gene/Holm/DCA tables, abstention coverage curve, leave-P53-out rows. |
| A3-CL-04 | Stability ddG is a controlled negative (homology-held-out calibration-RMSE inflation +0.04) but weak everywhere (Spearman 0.358--0.384, below supervised ThermoMPNN 0.6--0.7); the absence of a gap is expected, not reassuring, because a foundation model's pretraining is not truly OOD here. | Headline (honest null + hedge) | `N10_ddg/ddg_calibration_v2.json` | Inspect in-distribution vs held-out Spearman and calibration-RMSE inflation. |
| A3-CL-05 | Binding ddG (SKEMPI v2, n=3,424): overall Spearman 0.20 carried by fold destabilization, not interface recognition; antibody-antigen collapses to 0.08 vs general interfaces 0.22; non-interface ~0 (corroborative only, range-restricted); scaling to 3B does not rescue (0.148). | Headline (sharp boundary) | `N4_antibody_ddg/results.json`; `results_3b_wt.json` | Inspect Spearman by subset with CIs; 650M vs 3B overall. |
| A3-CL-06 | Developability: no BH-significant association is detected for polyreactivity (0/14) or aggregation (0/7); these denominators do not exclude small effects. AntiBERTy also has 0/7 aggregation discoveries and binding near 0, while recovering modest polyreactivity (3/14) that ESM-2 misses. On binding, 27/30 datasets are single-parent DMS where germline identity is constant, so the germline-distance control is undefined (cannot confound); 14 of those 27 are BH-significant raw associations (|rho| up to 0.69) reported as germline-unconfounded; only 3 multi-germline binding datasets exist and are collinear (control under-powered). | Bounded liability non-detection + corrected binding sub-claim | `N9_antibody/developability_corrected.json` (`summary.per_property`: binding `n_control_inapplicable_single_parent`=27, `single_parent_sig_raw_bh`=14; polyreactivity/aggregation `sig_raw_bh`=0); `antibody_lm_baseline.json` | Inspect `summary.per_property` counts and per-dataset `germline_constant`/`sig_raw_bh`. |
| A3-CL-07 | External context: ESM-1v zero-shot DMS rank-agreement (median |rho|=0.49) sits at/above a panel of published VEPs (REVEL 0.47, DeepSequence 0.45, VESPA 0.40, CADD 0.37, EVE 0.35) on the Livesey & Marsh compilation (10 proteins) -- situating, not ranking. | Supporting context | `a3_meta/vep_comparison.csv` | Inspect median |rho| per predictor; small external set. |
| A3-CL-08 | No single measurable covariate (length, sample size, score spread, class balance) predicts per-task reliability across DMS/ClinVar (best single Spearman 0.25 with CI touching 0; multivariate R^2 ~0.03); length and sample size flip sign between axes. | Headline (thesis: audit, not heuristic) | `A3_reliability_predictor.json` | Inspect per-predictor cross-axis Spearman and multivariate R^2. |
| A3-CL-09 | The uniform protocol caught and retracted three of its own prior over-claims. Each reversed a reportable conclusion. | Load-bearing method evidence | See A3-CL-09a/b/c | See sub-rows. |
| A3-CL-09a | Retraction 1 -- a clamped-confidence "toolkit reproducibility gap" (~0.05 AUROC) vanished on the raw LLR (raw gaps -0.003/-0.005): raw-score discipline. | Retraction | `N1_scale/expand80/scale_*_masked.json` (raw vs clamped) | Confirm raw sim gap is ~0, not ~0.05. |
| A3-CL-09b | Retraction 2 -- a long-protein (>1024 aa) scoring artifact masqueraded as per-gene clinical unreliability (MET 0.37, BRAF 0.39 "fail"); one fp32 consistent re-score (max repeat |dlog p|=0) raised them to 0.80/0.83: same-input-same-score determinism. | Retraction | `N6_clinvar/n6c_consistent.json` | Confirm MET/BRAF per-gene AUROC ~0.80/0.83 under the consistent re-score. |
| A3-CL-09c | Retraction 3 -- an "antibody binding is germline-confounded (18->4)" claim had auto-coded an inapplicable (undefined, zero-variance) germline control on single-parent DMS as a confound loss: applicability-gated de-circularization. | Retraction | `N9_antibody/developability_corrected.json` | Confirm `control_applicable=false` on single-parent (`germline_constant=true`) datasets; the 14 BH-significant binding rows sit there. |

## Explicit non-claims / overclaim guardrails

- Do **not** upgrade any axis to clinical deployability. The clinical axis is benchmark-only and
  explicitly base-model-leakage-bounded (ESM-2 pretraining subsumes the human proteome); it is not
  an ACMG-grade classifier.
- Do **not** claim architecture-independence beyond MLM encoders. The cross-encoder check covers
  ProtBERT/ESM-1v (sequence-only MLMs) only; autoregressive and structure-conditioned PLMs are
  untested, and ESM-1v was itself purpose-trained for variant effect.
- Do **not** state the moderate discrimination ceiling as *intrinsic*. It is the achievable level
  for a single sequence-only zero-shot pass; ESM-family co-distillation (VESM) and structure-aware
  tokenisation (SaProt) raise it, and are out of this audit's scope.
- Do **not** read the stability "controlled negative" as reassurance. The honest reading is
  "stable yet weak" -- a homology hold-out is not genuinely OOD for a foundation model, so the
  absent gap is expected.
- Do **not** describe non-interface binding ~0 as independent evidence; it is corroborative only,
  because its narrow binding-ddG range restricts the Spearman regardless of model quality.
- Do **not** restore the withdrawn "antibody binding is germline-confounded" framing. On
  single-parent DMS the germline control is undefined, not passing; the 14 BH-significant binding
  associations are reported as germline-unconfounded, and the residual caveat is the 3 collinear
  multi-germline datasets, not a confound.
- Do **not** present Fig 11 as new data or a B1-style independent-dataset meta-analysis; A3's tasks
  are heterogeneous, so it is a task-boundary synthesis of Figs 1,3--5,7--8 and S9 only
  (post-demotion spine; do not claim continuous F1–F11 main).
- Do **not** present the atlas as a physical-world statistical law. Locked physical readouts are
  rank correlations (SKEMPI Spearman with bootstrap CIs) and in-panel ECE tokens. A law paper
  would need paired experimental kcal/mol measurements on a physical hold-out, not another
  AUROC or ECE table. Full-20 per-bin reliability coordinates were never written; do not
  synthesise them by re-fitting Platt.
- **Demotion map (admin):** ex-F2→S6, ex-F6→S7, ex-F9→S8, ex-F10→S9. Main faces remain
  F1, F3, F4, F5, F7, F8, F11.
- Do **not** claim any covariate predicts reliability; the paper's point is the opposite (measure
  per task).

## No-upload / no-cross-workspace confirmation

The commands above are local read/write commands inside `biochem-research` and do not upload data,
reserve DOIs, contact journal systems, or read/mutate any sibling workspace. All inputs are already
committed under `research/results/`; if a path is absent in a given worktree, treat it as an
integration-state issue and stop with a missing-input note rather than re-deriving or downloading
the missing table.
