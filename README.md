# A3: protein language model reliability-boundary atlas

Public research-release export of **Reliability boundaries of protein language
model zero-shot variant scores**, by Zeyu Fu (2026).

Repository: https://github.com/PeterPonyu/a3-plm-boundary-atlas

Version archive: https://doi.org/10.5281/zenodo.23101189 . The existing concept
series is preserved. The older `10.5281/zenodo.21869995` record is historical,
not this release. Current manuscript availability text binds the code and
archive identifiers to this author-review version. Anonymous published access
and downloaded ZIP identity are verified separately in the release receipts;
a DOI written in a file alone is not proof of publication.

Scientific status remains
`AUTHOR_REVIEW_SAME_SCORE_INTEGRATED_NOT_SUBMISSION_READY` (**HOLD**).
Byte integrity and reproducibility checks are not independent scientific approval,
external clinical validation, clinical-use certification or visual acceptance.

## What is included

- Actual current scoring/analysis code, **4,677 unique variant scores across 20
  genes**, benign/pathogenic labels, mutation identities, exact window bounds,
  20 reference-sequence hashes/lengths, and immutable model-run provenance.
- The current `analysis/revision_20260923/clinical-unified/` results: 4,677 OOF
  predictions, per-gene discrimination/calibration, reliability bins, selective
  prediction and decision-curve tables, and **1,000 recorded identity-grouped
  refitting draws for each of five maps** (5,000 saved metric records).
- Current Article (29 pages) and SI (10 pages), LaTeX sources, consumed figure
  producers and assets, graphical-abstract vector masters, structural illustration
  dependencies, tests and environment records. The current Article was rebuilt
  after version-specific DOI/DAS integration; the original source-capsule hashes
  remain provenance, not the new public PDF identity.
- Locked historical nonclinical result/plot snapshots and separately labelled
  historical clinical/gnomAD/predictor/calibration records. See
  [EVIDENCE_BOUNDARIES.md](EVIDENCE_BOUNDARIES.md); historical tables are not a
  second current clinical score version.

The explicit source allowlist is `provenance/export-allowlist.json`. All 256
method-locked files and all 237 current manuscript build inputs were included.
No whole-analysis-directory copy was used. Official model weights, checkpoint
parts, raw upstream datasets, private submission bundles/cover letters/admin,
backups, old compiled PDFs, caches and compilation temporary files are excluded.

## Portable offline verification (default; no model inference or fitting)

From the repository root, with Python 3.10 or later:

```bash
python3 -B tools/verify_public_release.py
```

This uses only Python's standard library. It checks the current public file
inventory and SHA-256 values, the 4,677-row identity/label/reference-hash contract,
per-gene and pooled AUROC from saved scores, saved OOF metrics and bins, the 5,000
recorded refit metrics/percentiles, paired Brier differences and current plotting
values. It neither imports PyTorch nor calls the calibrator fitting routines.
It performs no downloads, model forwards, new bootstrap/refit runs or writes.

The current public authority is `PUBLIC_RELEASE_MANIFEST.json` plus
`SHA256SUMS.txt`. Source-capsule `reproducibility/method-contract.json`,
`manuscript/release-current.json` and `manuscript/release-spec.json` are preserved
byte-for-byte as **source-version provenance**, not a substitute for the fresh
public inventory. After final DOI/DAS/PDF edits, the coordinator must refresh the
public manifest and rerun verification; stale checksums fail closed.

Offline exported tests (requires the already available pytest/analysis packages;
no install command is run by the verifier):

```bash
python3 -B -m pytest -q -p no:cacheprovider tests/test_public_release.py reproducibility/tests/test_a3_method_integrity.py
```

Additional original tests are retained. Tests that require excluded historical
submission files, compiler logs, private upstream inputs, a model runtime or
create PDFs are **not** part of the default no-inference/no-fit route. Their
availability is documented in [REPRODUCING.md](REPRODUCING.md), not hidden by skips.

## Reproduction tiers

1. **Verify saved evidence:** the default command above is portable and offline.
2. **Inspect current PDFs:** `manuscript/manuscript.pdf` and `manuscript/SI.pdf`.
3. **Rebuild from saved tables/assets:** see `REPRODUCING.md` and
   `tools/build_public.sh`. Requires existing R/Python/LaTeX dependencies and
   regular/bold Arial; no model scoring is needed. The public builder does not
   invoke historical submission packaging or refit clinical estimates.
4. **Rerun clinical analysis / full model scoring:** an explicit separate
   scientific action, not verification. Full inference requires the official
   weights and separately acquired upstream clinical/reference inputs. See
   `tools/score_clinvar.py` and `DATA_SOURCES.md`.

**Full reference sequences are not present in the permitted source capsule.**
This export includes their exact hashes and lengths, not fabricated FASTA files.
The 20 source clinical CSV hashes are retained for upstream input matching.
Obtaining the full matching reference/clinical input bytes is a genuine
prerequisite for model-inference reproduction.

Current within-panel results include mean per-gene AUROC 0.7840632 and
gene-resampling interval [0.7428235, 0.8237764]. The five-fold OOF and refitting
protocol is conditional on this panel's gene/class composition; it is not
gene-disjoint, position-disjoint or external validation. Historical antibody
interval identities remain unresolved, and gnomAD was not newly rescored.

## Rights and attribution

Author-written code: MIT, copyright **Zeyu Fu 2026** (`LICENSE`). Author-created
documentation, manuscript text, figures and derived tables: CC BY 4.0
(`LICENSE-DOCUMENTATION.md`). These grants cover only the author's contributions;
they do not relicense upstream labels, reference sequences, structures, datasets,
model weights, third-party software or fonts. Read `THIRD_PARTY_NOTICES.md` and
`DATA_SOURCES.md` before retrieving or redistributing upstream inputs.
