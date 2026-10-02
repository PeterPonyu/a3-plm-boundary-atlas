# Reproduction without private-workspace layout

## Default verification

`python3 -B tools/verify_public_release.py` needs only Python >=3.10. It checks
saved evidence and public checksums. It does not regenerate weights, scores,
bootstrap samples or fitted calibration maps. It writes nothing. Moving or
copying the release directory must not require `labs/active`, a capsule sibling,
the original absolute path or a network connection.

`python3 -B -m pytest -q -p no:cacheprovider tests/test_public_release.py reproducibility/tests/test_a3_method_integrity.py`
runs public contracts and the original 200-draw identity-isolation tests without
new paper fits. The tests require an existing pytest installation. Additional
original numerical/figure tests may require the recorded dependencies. No
`test.skip`, placeholder test or missing dependency is a verification success.

Original `test_current_contract.py`'s no-log-warning test needs compiler logs,
which are excluded as temporary files; it becomes runnable after a new compile.
Original private-submission/package/font-layout tests need excluded submission
files and are historical tests, not the default public suite. Original scoring
tests need PyTorch and upstream clinical/reference CSVs. Original GA export tests
create PDFs, and are not run by the no-PDF-edit public export workflow.

## Rebuild from saved tables and assets (explicit)

Run `bash tools/build_public.sh` only in a working copy and only when PDF
regeneration is intended. This verifies public inputs before writes, exports
clinical plotdata from saved outputs, renders the current R figures and native
graphical abstract, rebuilds the binding illustrations, compiles Article/SI,
and refreshes the public inventory. It does **not** call `validate_clinical.py`
(which refits OOF), `analyze_clinical.py`, any model scorer, a private submission
packager, or the original hash-accepting `--lock-inputs` operation.

Requirements already installed on the build host must include the Python/R
packages named in `environments/DEPENDENCIES.md`, `latexmk`/pdfLaTeX and document
packages, `pdffonts`, and Arial regular/bold. Fonts are not bundled. The original
renderers refuse fallback. No dependencies are installed by this build. PyMOL
is optional: default structure plotting consumes bundled existing PNGs.

The coordinator performs the final DOI/DAS injection and compilation for this
release. This sidecar did not execute the PDF-producing builder or modify PDFs.
Source/PDF binding after DOI edits is a release-coordinator responsibility.

## Rerun current clinical analysis (new scientific computation)

The exact completed implementation is
`analysis/revision_20260923/analyze_clinical.py`, with imports of
`analysis/xfit_calibration_20260918/run_xfit_calibration.py` and `calibrators.py`.
It is exported unchanged and refuses to overwrite an existing output directory.
An intentional new analysis is not checksum verification:

```bash
python3 analysis/revision_20260923/analyze_clinical.py --scores analysis/revision_20260923/row-scores/scores.csv --out /absolute/new-analysis-directory
```

Do not execute this merely to release the existing results. The original
`validate_clinical.py` also calls `crossfit`; the public read-only validator checks
stored OOF values and refit records instead, and does not claim it reran fits.

## Full model scoring

`analysis/revision_20260922/rescore_clinvar_candidate.py` is the actual completed
scorer, preserved byte-for-byte, including its historical workspace constants.
Use the portable adapter `tools/score_clinvar.py` to provide explicit upstream
locations without altering that scientific implementation:

```bash
python3 -B tools/score_clinvar.py --clinical-csvs /absolute/upstream/clinical_csvs --n6c-summary /absolute/upstream/n6c_consistent.json --n6c-scorer /absolute/upstream/n6c_consistent_rescore.py --preflight-only
python3 -B tools/score_clinvar.py --clinical-csvs /absolute/upstream/clinical_csvs --n6c-summary /absolute/upstream/n6c_consistent.json --n6c-scorer /absolute/upstream/n6c_consistent_rescore.py --checkpoint /absolute/official/esm2_t33_650M_UR50D.pt --device cuda --out /absolute/new-row-score-directory
```

The adapter requires exact recorded upstream clinical CSV, N6c summary, scorer,
manifest and reference identities before scoring. `--preflight-only` imports no
model runtime and performs no model forwards. Full scoring still requires the
official completed checkpoint and its matching `download-provenance.json` beside
it, the original PyTorch/fair-esm environment and adequate CPU/GPU resources.
Official weights may be acquired using the retained download helper only when
network retrieval is intentionally authorized. No helper is called by default.

No complete reference FASTA exists in this source-only export. The exact 20
hashes and lengths are in `row-scores/reference-identities.csv`; the matching
source sequences must be obtained upstream. There is no public raw-data locator
recorded for these prepared gene CSV bytes, which limits exact inference
reproduction until those upstream inputs are supplied.

## Refresh inventory after intentional publication edits

```bash
python3 -B tools/update_public_manifest.py
python3 -B tools/verify_public_release.py
```

The update command uses `apply_patch` for manifest edits. Before replacement it
checks that the protected scientific score/analysis bytes and all unmodified
source dependencies retain their source hashes. Only manuscript text/PDFs,
source-version release-manifest metadata and new public wrappers/docs may change
under this publication refresh; a changed result or scientific implementation
is a hard failure, not silently accepted. Original logs, caches and compile
temporary files are not part of the public inventory.
