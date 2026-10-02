# FIGURE-PROVENANCE — A3 PLM boundary atlas

## Current authority: 23 September 2026

The original working source/PDF pair in `manuscript/` is now the author-review
revision, bound by `release-current.json`. The same-score clinical chain has
been integrated but the manuscript is not submission-ready. The older
authority claims, registry tables and payload hashes below are historical.

Current Fig. 2 reads `F3_clinvar_20260923.csv`: newly inferred ClinVar scores
and 2,000-draw per-gene intervals. Only its gnomAD comparator is a historical
summary. Fig. 3 reads `clinical_*_20260923.csv`: five-fold same-score OOF
reliability, Brier baselines, gene ECE and paired Brier differences, with
1,000 identity-grouped refits. The old six-gene figure and coverage/model
comparison panels are not current evidence. Fig. 7 and the graphical
abstract use `clinical_macro_20260923.csv`, mean 0.7840632 and the
10,000-draw gene-resampling interval [0.7428235, 0.8237764]. The exporter is
`analysis/revision_20260923/build_clinical_plotdata.py`, called by the current
rebuild entrypoint; old plotdata were not overwritten.

On 30 September 2026, Fig. 7 was reduced to the two source-metric panels:
AUROC and Spearman. The shared null-gap scale and heuristic task-zone counts
were retired; overlapping binding subsets are not independent tasks. Starred
antibody archive nodes keep their unresolved 650M provenance status. No input
CSV or numerical estimate changed. Its current A/B topology is checked by
`scripts/verify_figure_release.py`; the six-panel registry snapshot below is
historical and cannot establish additional current measurements.

Other nonclinical figures preserve their observations. S2 remains a
historical gnomAD sensitivity, and S9 remains the historical predictor table
with missing clinical spread explicitly disclosed. The old calibration
outputs/registries do not define a second current score version. Successful
rendering and hashes do not certify scientific readiness. Visual inspection
is recorded separately in the hash-bound 30 September local review package;
machine rendering and bbox checks do not attest visual or scientific acceptance.
Do not use the historical rebuild/upload directions
below for the current revision; use `scripts/rebuild_current.sh`.

| Field | Value |
|---|---|
| **Date** | 2026-09-09 (visual-contract compositor rewrite) |
| **Carrier** | `manuscript/manuscript.tex` is the TeX source; the checked-in upload authority is `manuscript/submission/flat_Proteins/manuscript.pdf` (43 pages; SHA-256 `45ddad9b9231d5deba8900146c04a13840149ce06fe90eec234a8b67c2b43b4f`) with seven reader-facing main figures plus S1–S9. Optional companion: `a3_ten_figures_floats.tex`. SUPERSEDED: `f0d3fac473976871dbf85f04af0400e11cb407934431f15d1b57ccfc786372ef` (prior Carrier). SUPERSEDED official: `93a4dfafbe4f1c86904d54be007f0b0523b3e750aba24e0d28a0e23f338e324e` (incoming 2026-09-11 ungated face). |
| **Plan** | Historical internal planning record (not part of the submission package) |
| **Option** | Canonical composites regenerate only from committed canonical plotdata (no GPU; no new inference). |
| **Hard gates** | NO fabricated stats · NO unauthorized GPU · NO PanODE merge |

> **Science (2026-09-18 a3-science2).** Follow-up `A3-XFIT-CAL-R2-20260918`
> on the same folds. 7/20 pooled-Platt degradation is prevalence mismatch
> (Spearman $+0.815$), not small $n$; Holm-significant degraders are OTC and
> PIK3CA only; per-gene Platt does not significantly beat pooled. Official face
> `e79ece2efead8019e261adebc333b869ce8ffa9cc650341635ff5ba4ca0e0472`
> (45 pages). Previous official
> `7a52829c35c0813530262c6fcbf5de44b84d1bdf38d4c3d1192ab7739ae22654`
> (44 pages) is the a3-science xfit face. Not a restoration of $0.0207$ / 91\%.
>
> **Science (2026-09-18 a3-science).** New 5-fold cross-fitted calibration
> (`A3-XFIT-CAL-20260918`) printed beside the historical STOP record. Not a
> restoration of $0.0207$ / 91\%. Official face
> `7a52829c35c0813530262c6fcbf5de44b84d1bdf38d4c3d1192ab7739ae22654`
> (44 pages). Previous official
> `45ddad9b9231d5deba8900146c04a13840149ce06fe90eec234a8b67c2b43b4f`
> (43 pages) is the pre-xfit face.

> **Authority (2026-09-11 a3-impl).** Official face rebuilt after deleting the
> named full-20 Platt pair and splitting Methods so printed $0.2443$ is
> uncalibrated (raw, no Platt) only. SHA-256
> `45ddad9b9231d5deba8900146c04a13840149ce06fe90eec234a8b67c2b43b4f`
> (43 pages). Incoming official `93a4dfaf…324e` is SUPERSEDED. Isolated_r2
> 94 still attaches only to `a0f208d2` and is **not** transferred. Fresh
> isolated_r2 of this hash (2026-09-11 protein-visual-ultraqa): official
> 43 pages and 16 included figures all ≥95 (min page 95 / min figure 95);
> standalone graphical abstract 94 (not in the 43-pp face). Scores live in
> `.omc/research/research-20260911-protein-visual-ultraqa/`. `F7_20gene_reliability.csv` was **not** written. `F7_20gene_ece`
> remains unpromoted. Calibrated STOP is held. Fig 2A stays two-column. Fig 7B
> 0.165 axis is still the N4 `antibody_interface` cell (`n=525`).

> **Visual closure (2026-09-04).** The renderer-only candidate was regenerated in
> `/tmp/a3b28_audit/a3_candidate9` and promoted through the external 42-page build at
> `/tmp/a3b28_audit/a3_submission_candidate`. The audit closed the F5 title/tag header
> collision, the F2/F9 compact-card title/value collisions, the F8 and FS4 blocked-card
> text overflow, and the F11 compact-card/x-axis clipping. F5, F8, F11, and FS4 were
> re-inspected as final PDF/PNG pairs; all A--F labels, locked wording, and scalar values
> remain present. No plotdata, registry estimand, statistic, or scientific claim changed;
> only text wrapping, card spacing, and header margins changed. The B28 four-figure set
> was also inspected at native resolution; its labels and leader lines remain readable,
> so no B28 source change was warranted. The promoted authority package and checksum
> now bind the candidate manuscript hash above; the root draft PDF remains a separate,
> non-authoritative artifact by design.

## Historical metrics snapshot (2026-08-23)

The following values are retained as a planning-era provenance snapshot; they do not describe
the current checked-in upload authority above.

| Metric | Before | After |
|---:|---:|---:|
| `\begin{figure}` count | 16 | **16** |
| Candidate full-20 calibration face | retained on disk | **blocked from canonical release** by `F7_20gene_perbin_STATUS.txt` STOP verdict |
| PDF pages | ~28 (pre-campaign) | **35 (historical snapshot)** |
| PDF size | — | **2 523 171 bytes (historical snapshot)** |

**Authority boundary:** `F7_20gene_ece` is candidate-only because Platt-calibrated full-20 per-bin coordinates remain unrestorable (`F7_20gene_perbin_STATUS.txt` STOP). RAW-only bins live in `F7_20gene_raw_bins.csv` and are not that candidate face. `atlas-not-law` is a boundary/reliability synthesis, not an independent data axis or new experiment.

## Rebuild

```bash
cd manuscript
bash scripts/build_release_assets.sh /absolute/empty/staging-directory
```

The compatibility entry point `scripts/build_pdf.sh` delegates to the same
staged release helper. The helper is fail-closed: it requires Tinos, writes
only to the explicit empty staging directory, validates embedded fonts, and
checks that canonical figure hashes did not change. Legacy single/dual-panel
writers are quarantined in `scripts/R/LEGACY-WRITERS.md` and are never called.

## Figure → panel map

| Reader figure | File stem | Panels A–F (views of committed CSV) | CSV lock |
|---:|---|---|---|
| 1 | `F1_scale_plateau` | A AUROC vs size · B n discriminating · C ΔAUROC step · D fraction discriminating · E paired 3B−650M · F scoring-convention sensitivity | `F1_scale.csv`; E/F `F1_EF.csv` |
| 2 | `F3_clinvar_forest` | A forest · B ClinVar vs gnomAD · C AUROC vs n_path · D threshold fractions · E length–AUROC · F not estimable (no per-variant coordinates) | `F3_clinvar.csv`; E/F `F3_EF.csv` |
| 3 | `F4_clinical_calibration` | A reliability · B abstention · C ESM-2 vs ESM-1v · D n pathogenic calls · E P53 contrast · F coverage-range | `F7_reliability.csv`, `F7_coverage.csv`, `F7_headtohead.csv`; E/F `F4_EF.csv` |
| 4 | `F5_stability_ddg` | A discrimination Spearman · B calib RMSE · C inflation · D split sizes · E two-variant ΔSpearman · F two-variant ΔRMSE | `F10_discrimination.csv`, `F10_calibration.csv`, `F10_meta.csv`, `F10_reference.csv`; E/F `F5_EF.csv` |
| 5 | `F7_binding_boundary` | A Spearman forest · B n · C CI width · D key contrast · E not estimable (no binding-row coordinates) · F locked Ab–Ag CI width | `F4_binding.csv`; E/F `F7_EF.csv` |
| 6 | `F8_developability_map` | A frac BH-sig · B n_sig vs n_ds · C liability pooled · D multi-germline n · E/F not estimable (no dataset-level effects/variance) | `F5_developability.csv`; E/F `F8_EF.csv` |
| 7 | `F11_boundary_atlas` (`atlas-not-law`) | A AUROC atlas · B Spearman atlas · C gap from null · D zone counts · E atlas-separation · F provenance inventory | `F11_atlas_auroc.csv`, `F11_atlas_spearman.csv` (from `results-snapshot/a3_meta/`); E/F `F11_EF.csv`; synthesis only |
| S1 | `FS1_proteingym_atlas` | A forest · B histogram · C n_test vs AUROC · D pooled summary · E assay-CI coverage · F not estimable (no exposure labels) | `FS1_proteingym_atlas.csv`, `FS1_proteingym_pooled.csv`; E/F `S1_EF.csv` |
| S2 | `FS2_gnomad_control` | A linked sources · B scatter · C Δ histogram · D class means · E gene-level · F class-level | `FS2_clinvar_vs_gnomad.csv`, `FS2_class_means.csv`; E/F `S2_EF.csv` |
| S3 | `FS3_antibody_robustness` | A grouped forest · B overall · C antibody · D 3B−650M Δρ · E paired size contrast · F not estimable (no influence rows) | `FS3_antibody_robustness.csv`; E/F `S3_EF.csv` |
| S4 | `FS4_abstention_tiers` | A tiers · B coverage curve · C full−keep-hi · D shuffle−conf · E risk card · F not estimable (no protein resample) | `FS4_abstention_tiers.csv`, `FS4_abstention_coverage.csv`; E/F `S4_EF.csv` |
| S5 | `FS5_vep_comparison` | A median \|ρ\| · B class mean · C rank · D n proteins · E range · F missingness | `FS5_vep_comparison.csv`; E/F `S5_EF.csv` |
| S6 | `F2_cross_model` | A model AUROC · B agreement · C Δ vs ESM · D family mean · E pair · F direction | `F2_crossmodel.csv`, `F2_agreement.csv`; E/F `S6_EF.csv` |
| S7 | `F6_protein_uq` | A mean ρ_unc · B BO · C summary ρ · D UCB beats greedy · E shift · F negative-control | `F9_rho_unc.csv`, `F9_bo_acquisition.csv`, `F9_summary.csv`; E/F `S7_EF.csv` |
| S8 | `F9_antibody_lm` | A liability frac · B Ab–Ag ρ · C BH counts · D CI width · E effect · F FDR | `F8_antibody_n9.csv`, `F8_antibody_n4.csv`; E/F `S8_EF.csv` |
| S9 | `F10_reliability_predictor` | A score-spread · B length · C class balance · D cross-axis ρ · E partial · F missing | `F10_reliability_predictor.csv`; E/F `S9_EF.csv` |

## Honest gaps / blockers (documented, not filled)

| Item | Handling |
|---|---|
| **7-face main spine vs 10×4** | Spine kept (F1, F3–5, F7–8, F11 + S1–S9). **16 floats ≥ 10**; campaign standard met without promoting SI to main. |
| **Page growth** | ~28 → **35** pages under Wiley double-space / line-number draft. Journal page budget not enforced here; flag for human before portal. |
| **Panel C/D are alternate views** | Extra panels are deterministic transforms of the same committed tidy tables (deltas, counts, ranks, group means) — **not** new experiments or invented CIs. |
| **F11 / FS5 CSVs** | Canonical locked plotdata are read directly; the release helper never silently refreshes or copies inputs. |
| **Legacy plotdata name skew** | `F4_binding`→Fig7, `F5_developability`→Fig8, `F7_*`→Fig4, `F10_*`→Fig5/S9, `F8_*`→S8, `F9_*`→S7 — unchanged; captions name the CSV locks. |
| **Fig 2 A 20-gene forest** | Twenty ClinVar names need ~137 pt in a ~70–82 pt 3×2 cell. Split into two sub-columns of ten (identity preserved; measured 10.54 pt pitch, 3.90 pt gaps). A single-column 20-label forest remains structurally blocked in an equal cell; do not paper over. |
| **Atlas 0.165 antibody-interface axis** | Locked display cell in `F11_atlas_spearman.csv`; parent `results-snapshot/N4_antibody_ddg/results.json` key `antibody_interface` (n=525, ρ=0.16455 rounds to 0.165). Not a Table 2 / `F4_binding.csv` partition. Cited in Fig 7 caption, Discussion, Methods, and `scope_atlas_not_law.tex`. |
| **FS1 ProteinGym 19-id forest** | Equal 3×2 cell overlapped by 1.22 bp/pair (19×6.64 pt glyphs in ~104 pt). Two-sub-column split collapses a 156 pt cell (18-char ticks ~62 pt). Resolved by `a3_compose_tall_a`: panel A spans two rows and two columns of the same 6.5×4.3 in canvas. Identifiers and type scale unchanged. |

## Self-refutation record (unchanged)

See `figure_ledger.md`: clamped-confidence toolkit gap; long-protein scoring artifact; inapplicable germline confound.
