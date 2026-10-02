# Legacy A3 renderers — quarantine

The canonical release path is only:

```bash
bash scripts/build_release_assets.sh /absolute/empty/staging-directory
```

It invokes `make_four_panel_all.R`, `make_f7_20gene_ece.R`, and
`make_a3_graphical_abstract.R`. The other `make_*.R` files in this directory
are retained as historical provenance only. They are not release writers and
must not be invoked by a release helper: several predate the outside-legend,
explicit-staging, and strict-Tinos contract.
