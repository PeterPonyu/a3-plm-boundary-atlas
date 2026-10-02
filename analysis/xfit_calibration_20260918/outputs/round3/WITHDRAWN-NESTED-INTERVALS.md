# Nested-refit uncertainty withdrawn

The nested-bootstrap intervals in this historical directory are invalid as held-out refitting uncertainty: bootstrap copies of original examples crossed fitting/evaluation folds. All 200 original draws exhibit the problem. The numerical ledger is preserved for traceability, not for reuse as a corrected result.

Current code groups copies by original example identity and writes any future estimates to a separate candidate directory. Identity regression tests pass; no replacement nested-refit intervals have been calculated or validated. Per-gene tables based on saved unique-row OOF predictions are not withdrawn by this specific duplicate-copy finding.

See ../../../../reproducibility/METHOD-AUDIT.md. Historical STATUS/PASS statements do not supersede this withdrawal.
