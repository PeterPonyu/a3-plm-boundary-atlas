# Data-source locators, versions and access boundaries

This is a transcription of **recorded** source provenance, not a new network
audit. Dataset fetch dates, release IDs, accession numbers and license assertions
not present in the capsule are not invented. A publication DOI identifies a source
study, not necessarily an exact downloadable dataset version. Obtain raw inputs
from the source custodians and follow their applicable terms.

## ProteinGym

Recorded source citation: Notin et al., ProteinGym (NeurIPS 2023), in
`manuscript/manuscript.tex`, reference `Notin2023`:

https://proceedings.neurips.cc/paper_files/paper/2023/hash/cac723e5ff29f65e3fcbb0739ae91bee-Abstract-Datasets_and_Benchmarks.html

Recorded sampling frames: 217 DMS assays attempted, 201 scored after QC, 173
proteins and 79,220 variants for the full panel; a separate 69-protein scale
subset. Exact raw benchmark archive version/download date is **not recorded in
the permitted capsule**. Included outputs are saved derived plotdata and
snapshots, not raw assay archives or a promise of complete row-level DMS inputs.

The old clinical manifests have schema `b01_b02_manifest_v1`, source
`proteingym+esm2_650M_fairesm`, source version
`ESM-2 650M wt-marginal (fair-esm)`, and source URL
https://github.com/facebookresearch/esm. That URL locates the model, not the
clinical dataset. Generic DMS labels in those historical metadata fields do not
override the actual current scorer's ClinVar label interpretation.

## ClinVar

Recorded source citation: Landrum et al. (2018), DOI
https://doi.org/10.1093/nar/gkx1153 (`Landrum2018`). The capsule does **not** record
an exact ClinVar release date/archive accession for this panel. The methods select
missense variants with exact pathogenic/benign significance terms, exclude
conflicting identities and isoform WT-residue mismatches, and use benign=1 for
LLR discrimination and pathogenic=1 for calibration.

Included: 4,677 labelled variants, 20 gene identities, all score/window fields,
clinical source-file hashes and 20 WT reference hashes/lengths. Provenance:
`analysis/revision_20260923/row-scores/provenance.json`, `run-binding.json`, and
`reference-identities.csv`. Original upstream inputs were gene-named CSVs under
`research/results/N6_clinvar/gof_test/clinical_csvs/`; that recorded historical
locator is not an accessible public dataset URL. **Full reference sequences and
source clinical CSV bytes are absent from this allowed capsule** and must be
retrieved separately. Hash matching, not a guessed isoform, determines identity.

## gnomAD

Methods record **gnomAD v4**, missense variants with allele count at least two,
excluding sites labelled pathogenic in the same gene. Exact v4 subrelease,
download date and raw archive identity are **not recorded**. Source citation:
Karczewski et al. (2020), https://doi.org/10.1038/s41586-020-2308-7
(`Karczewski2020`); this older study DOI is not asserted to be a v4 archive ID.

Only `FS2_clinvar_vs_gnomad.csv`, `FS2_class_means.csv`, the old `F3_clinvar.csv`
and the current figure's explicitly historical gnomAD column are retained.
They are historical summary-level sensitivity evidence. No gnomAD row scores
were newly generated, and no raw population variants are distributed here.

## MegaScale and ThermoMPNN

Recorded source studies: Tsuboyama et al. (2023),
https://doi.org/10.1038/s41586-023-06328-6; Dieckhaus et al. (2024),
https://doi.org/10.1073/pnas.2314853121 (`Tsuboyama2023`, `Dieckhaus2024`).
The methods use the MegaScale reliable single-mutation subset and ThermoMPNN
MMseqs2 clusters at 25% identity. Exact archive/revision and access dates are
**not recorded**. The response is dataset-native `ddG_ML`, not a universal
kcal/mol conversion. Only locked derived stability summaries/plotdata are
included (`F10_*`, `F5_EF`, related figures); raw measurements, split archives and
reusable fit coefficients are not supplied by these aggregate snapshots.

## SKEMPI v2 and illustrative structure

Recorded benchmark version: **SKEMPI v2 / SKEMPI 2.0**. Source study:
Jankauskaite et al. (2019), https://doi.org/10.1093/bioinformatics/bty635
(`Jankauskaite2019`). Exact raw file revision/access date is **not recorded**.
The retained analysis concerns 3,424 validated single mutations. Archived
650M subset summaries, 3B sensitivity scores, robustness outputs and author
plotdata are included. They are not a complete raw SKEMPI distribution.

The consumed illustrative structure is **1A22**, recorded locator
https://files.rcsb.org/download/1A22.cif. Exact bytes and source-table binding are
in `manuscript/objects/provenance.json`. Illustration is human growth hormone
G120R chain A / receptor chain B; existing I58A and I58L records are contextual
examples, not new mutant-coordinate generation or an independent cohort.

Antibody archived intervals differ: the N4 summary retains approximately
[0.002, 0.157], whereas the separate antibody-language-model ESM-2 snapshot
retains [-0.0022, 0.1597]. Row/scorer/resampling identities are unresolved;
they are not pooled or used to establish subgroup significance.

## FLAb

Recorded source study: Chungyoun, Ruffolo and Gray (2024),
https://doi.org/10.1101/2024.01.13.575504 (`Chungyoun2024FLAb`). Exact dataset
archive version and access date are **not recorded**. The panel contains 63
datasets and seven properties. Only author-derived discoveries, applicability
counts and plotting snapshots are included. Acquire original experiments from
FLAb and their cited providers under their own terms. Single-parent designs
with constant germline identity do not have an applicable partial-correlation
control; inapplicability is not a failed discovery or a deconfounding result.

## Model checkpoint and other retained summaries

Official recorded weight locator:
https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt

The immutable download record has length 2,604,537,549 bytes and locally computed
SHA-256 `ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c`.
It is an integrity binding, **not** an independently published upstream signature.
Official weights and all 156 parts remain upstream/not bundled.

The historical external VEP comparison cites Livesey and Marsh,
https://doi.org/10.15252/msb.20199380, with the old capsule artifact index naming
Figshare `38248524`. The exact source file/download date is unrecorded here;
only `results-snapshot/a3_meta/vep_comparison.csv` and derived plotting values
are shipped. Historical UQ plotdata cite Greenman et al.,
https://doi.org/10.1371/journal.pcbi.1012639, through the manuscript bibliography.
No additional source archive is inferred from these citations.
