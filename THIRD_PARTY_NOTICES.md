# Third-party rights and research-input boundaries

The MIT grant in `LICENSE` applies to author-written code only. CC BY 4.0 in
`LICENSE-DOCUMENTATION.md` applies to the author's research-material contributions
only. Neither notice is an assertion that all inputs are owned by the author.

- **ProteinGym / ClinVar-derived labels:** retained only as the identified
  4,677-variant analysis panel and historical derived manifests. The old manifest
  field says `ProteinGym MIT; ESM-2 MIT`; that is a copied provenance statement,
  not a newly verified blanket license for all original assays, clinical data or
  reference sequences. Consult the specific upstream release and original studies.
- **gnomAD:** only archived summary-level comparators are shipped. No raw
  population dataset, row-level score recovery or new gnomAD scoring is included.
- **MegaScale / ThermoMPNN:** aggregate derived tables only. The raw experimental
  dataset and cluster assignments must be obtained from their source records.
- **SKEMPI v2 / FLAb:** derived research snapshots, plotting tables and the
  pre-existing SKEMPI 3B sensitivity table are included; no complete raw upstream
  benchmark is claimed. Source observations remain under upstream terms.
- **PDB / RCSB structure 1A22:** the previously included `manuscript/objects/1a22.cif`
  and author-generated illustrative assets are retained as a consumed figure
  dependency. The source locator and exact hash are in
  `manuscript/objects/provenance.json`. The structure is not newly authored or
  relicensed by this release. No mutant conformation is claimed.
- **ESM / fair-esm and other models:** official weight binaries, model checkpoint
  parts and caches are excluded. The fair-esm retrieval helper and original
  download integrity metadata are author-code/provenance, not model distribution.
  No external or upstream signature is claimed for a locally computed weight hash.
- **Software and fonts:** no third-party package implementation or Arial font
  files are bundled. Python/R/LaTeX/PyMOL/RDKit and fonts require their own
  installations and licenses. Embedded PDF fonts do not grant permission to
  redistribute standalone font files. Arial regular/bold are required by the
  original rendering contract; do not silently substitute a different font.

`DATA_SOURCES.md` records only locators and version statements actually present
in capsule provenance/manuscript citations. Missing release IDs or access dates
remain explicitly unrecorded. Public availability alone is not a license grant.

This package contains no patient-level private submission material, cover letters,
administrator records, credentials, official weights or raw upstream archives.
