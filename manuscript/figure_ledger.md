# A3 Figure Ledger — one claim per figure, provenance, regeneration

This ledger pins each A3 figure to one manuscript claim and to the committed result artifact used to regenerate it. Scope is retrospective, public-data, computational-only benchmark evidence on ESM-2 zero-shot variant-effect scores; no clinical-deployment, wet-lab, or prospective claim is made.

**Canonical reader-facing figure spine (independent of legacy logical asset IDs):**

| Reader-facing number | Logical asset ID |
|---:|---|
| 1 | `F1_scale_plateau` |
| 2 | `F3_clinvar_forest` |
| 3 | `F4_clinical_calibration` |
| 4 | `F5_stability_ddg` |
| 5 | `F7_binding_boundary` |
| 6 | `F8_developability_map` |
| 7 | `F11_boundary_atlas` |

Supplement labels remain **Fig S1–S9** and retain their legacy logical stems. SI demotion map:
**ex-F2→S6** (cross-encoder), **ex-F6→S7** (UQ), **ex-F9→S8** (AntiBERTy),
**ex-F10→S9** (reliability predictor). File basenames may still carry legacy `F2`/`F6`/`F9`/`F10`
prefixes; manuscript captions and in-text refs use the reader-facing sequence above.

**Figure authority (2026-09-03).** The release is built from *data plus renderer*, never from the
committed `figures/*.pdf`: those image files are **not** release build inputs. The hash-bound
closure is the 56 `figures/plotdata/*.csv` artifacts plus the R renderers, and
`build_submission_pdf.sh` regenerates every figure into external staging at release time. The
`figures/*.pdf` on disk are therefore a **dev-loop cache** for `build_pdf.sh`, not an authority.
On 2026-09-03 all 16 cached figures had drifted from renderer output, so a quick local build
displayed different panels than the package that would ship; the cache was regenerated from the
renderer and now matches byte-for-byte. If they diverge again, the renderer wins — refresh the
cache with `A3_FIGURE_OUTPUT=<tmp> Rscript scripts/R/make_four_panel_all.R` and copy back, never
the reverse.

**Rendering pipeline.** Statistical figures regenerate via R/ggplot2 in the required **Tinos**
font — `scripts/export_plot_data.py` reshapes committed result JSONs into tidy CSVs under
`figures/plotdata/`, while the canonical renderers are `scripts/R/make_four_panel_all.R`
and `make_a3_graphical_abstract.R`. The retained full-20 candidate renderer is explicitly
candidate-only and is never called by a canonical builder. LaTeX uses `newtx`.
**Source of truth = `manuscript.tex` (hand-maintained LaTeX); no markdown intermediary.**
`bash scripts/build_release_assets.sh /absolute/empty/staging-directory` is the fail-closed
asset entry point; it never overwrites canonical assets. Legacy renderers are quarantined and
not called. Figures are floated `[htbp]` inline at their first in-text reference.

## Main figures

| Figure | One claim carried | Direct committed inputs | Output | Scope guardrail |
|---|---|---|---|---|
| Fig 1 / `F1_scale_plateau` | ESM-2 DMS variant discrimination rises 35M→650M then plateaus; 3B gives no gain (efficient ceiling at 650M). | `research/results/N1_scale/expand80/scale_{35M,150M,650M,3B}_masked.json` (n=69) | `figures/F1_scale_plateau.{png,pdf}` | Masked-marginal, n=69 de-pseudo-replicated subset; full-217 benchmark corroborates 650M≈0.70. |
| Fig 2 / `F3_clinvar_forest` | Zero-shot LLR discriminates ClinVar pathogenic vs benign per gene (mean 0.78; 19/20 nominal CI exclude chance; 16/20 ≥0.70). | `research/results/N6_clinvar/n6c_consistent.json` (one fp32 consistent re-score) | `figures/F3_clinvar_forest.{png,pdf}` | Base-model leakage (ESM-2 pretraining includes human proteome); benchmark-only, NOT a clinical classifier. Supersedes the withdrawn long-protein-artifact N6 per-gene picture. |
| Fig 3 / `F4_clinical_calibration` | Clinical axis is calibratable (pooled ECE drop) and supports selective abstention on the study panel. | `research/results/N6_clinvar/clinical_calibration_full20.json`; `clinical_calibration.json` | `figures/F4_clinical_calibration.{png,pdf}` | In-panel calibratability only; not ACMG-grade / not transport to new genes. |
| Fig 4 / `F5_stability_ddg` | Stability ΔΔG is a controlled negative — stable across homology hold-out, but weak everywhere ("stable yet weak"). | `research/results/N10_ddg/ddg_calibration_v2.json` | `figures/F5_stability_ddg.{png,pdf}` | Homology hold-out is not truly OOD for a foundation model; absent gap is expected, not reassuring. |
| Fig 5 / `F7_binding_boundary` | ESM-2 LLR tracks binding ΔΔG only as a weak fold-destabilization proxy (overall 0.20), collapsing for antibody–antigen (0.08) and non-interface (≈0). | `research/results/N4_antibody_ddg/results.json` (1000× bootstrap CIs) | `figures/F7_binding_boundary.{png,pdf}` | Non-interface ≈0 is corroborative only (range-restricted). Benchmark-only. |
| Fig 6 / `F8_developability_map` | No BH-significant liability associations (polyreactivity 0/14, aggregation 0/7); binding sub-claim is germline-unconfounded on single-parent DMS (14/27). | `research/results/N9_antibody/developability_corrected.json` (BH q<0.05; corrected) | `figures/F8_developability_map.{png,pdf}` | Do **not** restore withdrawn "binding germline-confounded" framing. |
| Fig 7 / `F11_boundary_atlas` | `atlas-not-law`: one-glance boundary/reliability synthesis of main axes + S9; moderate naturalness/pathogenicity discrimination cannot be generalized into a law, and this is not an independent data axis. | `research/results/a3_meta/a3_atlas_{auroc,spearman}.csv` | `figures/F11_boundary_atlas.{png,pdf}` | Task-boundary synthesis, not B1-style independent-dataset meta-analysis or a new experiment. |

## Supplement figures (including demoted former main faces)

| Figure | Legacy file id | One claim | Direct committed inputs | Output |
|---|---|---|---|---|
| Fig S1 | `FS1_proteingym_atlas` | Position-held-out ProteinGym calibration sub-study (19 assays). | B05 / ProteinGym calib artifacts via export pipeline | `figures/FS1_proteingym_atlas.{png,pdf}` |
| Fig S2 | `FS2_gnomad_control` | Clinical transfer robust to ClinVar-benign vs gnomAD negative set. | `N6_clinvar/n6c_consistent.json` (+ gnomAD control rows) | `figures/FS2_gnomad_control.{png,pdf}` |
| Fig S3 | `FS3_antibody_robustness` | Antibody/interface binding boundary robust to scale and scoring. | `N4_antibody_ddg/results.json`; `results_3b_wt.json` | `figures/FS3_antibody_robustness.{png,pdf}` |
| Fig S4 | `FS4_abstention_tiers` | Clamped confidence is anti-helpful for selective abstention. | N1 / abstention tier artifacts via export pipeline | `figures/FS4_abstention_tiers.{png,pdf}` |
| Fig S5 | `FS5_vep_comparison` | ESM-1v zero-shot DMS rank-agreement sits among/above published VEPs (situating, not ranking). | `a3_meta/vep_comparison.csv` (Figshare 38248524) | `figures/FS5_vep_comparison.{png,pdf}` |
| Fig S6 (ex-F2) | `F2_cross_model` | Signal shared across MLM encoders (ProtBERT ≈ ESM-2); AR/structure untested. | `research/results/N5_arch/N5_summary.json` | `figures/F2_cross_model.{png,pdf}` |
| Fig S7 (ex-F6) | `F6_protein_uq` | Learned UQ under shift is weak; UCB acquisition never beats greedy (0/7). | `research/results/N8_uq/` | `figures/F6_protein_uq.{png,pdf}` |
| Fig S8 (ex-F9) | `F9_antibody_lm` | AntiBERTy does not rescue aggregation/binding; modest polyreactivity only. | `research/results/N9_antibody/antibody_lm_baseline.json` | `figures/F9_antibody_lm.{png,pdf}` |
| Fig S9 (ex-F10) | `F10_reliability_predictor` | No single tested covariate predicts per-task reliability in this sample. | `research/results/A3_reliability_predictor.json` | `figures/F10_reliability_predictor.{png,pdf}` |

## Self-refutation record (carried in Discussion)

The uniform protocol caught three of its own over-claims, each reversing a reportable conclusion:
1. **Clamped-confidence "toolkit gap" (N1)** → raw-score discipline; raw gap −0.003.
2. **Long-protein per-gene "unreliability" (N6)** → same-input-same-score determinism (fp32); MET 0.37→0.80, BRAF 0.39→0.83.
3. **"Binding germline-confounded" (N9)** → applicability-gated de-circularization; control undefined on single-parent DMS.

## Verification

```bash
bash paper/A3-plm-boundary-atlas/manuscript/scripts/build_release_assets.sh /absolute/empty/staging-directory
```

All inputs are committed under `research/results/`; no network or GPU is needed.
