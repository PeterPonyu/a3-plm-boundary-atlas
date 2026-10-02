#!/usr/bin/env python3
"""Export tidy plot-ready CSVs for A3 supplementary figure FS1 — the ProteinGym
(B05) variant-effect reliability atlas. The pooled AUROC is quoted in the A3 text
but the per-assay atlas is not visualized anywhere. Reshapes the committed result
artifact research/results/C_proteingym/panel.json into long CSVs under
figures/plotdata/. base csv only; no network/GPU; deterministic.

Source schema (panel.json):
  per_target[]  : {target, verdict, n_test, auroc, auroc_ci[lo,hi], ...}  (19 assays)
  pooled        : {n, auroc_calibrated, ...}
  pooled_ci_95  : {auroc_calibrated[lo,hi], ...}
The pooled AUROC reported here is the held-out *calibrated* AUROC the A3 text
quotes as 0.781 [0.757, 0.803] on n=1551 variants.
"""
import json, os, csv

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
RES = os.path.join(ROOT, "research", "results")
OUT = os.path.join(os.path.dirname(__file__), "..", "figures", "plotdata")
os.makedirs(OUT, exist_ok=True)


def write(name, header, rows):
    p = os.path.join(OUT, name)
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    print("wrote", os.path.relpath(p, ROOT), f"({len(rows)} rows)")


panel = json.load(open(os.path.join(RES, "C_proteingym", "panel.json")))

# ---- FS1 per-assay atlas (19 ProteinGym DMS assays) ----
# A CI that excludes 0.5 == the assay discriminates (chance-beating with 95% conf).
rows = []
for t in panel["per_target"]:
    lo, hi = t["auroc_ci"]
    # strip the trailing DMS-study suffix for a compact, readable assay label
    label = t["target"]
    discriminates = lo > 0.5
    # NB: the upstream result's per-target "verdict" is a stale label from an
    # unrelated pipeline (it reads KILL_OR_REDESIGN even for AUROC>0.9 rows) and
    # is never rendered; drop it so the exported atlas carries no orphan column.
    rows.append([label, t["auroc"], lo, hi, t["n_test"],
                 int(discriminates)])
# sort ascending by AUROC for a clean forest/atlas
rows.sort(key=lambda r: r[1])
write("FS1_proteingym_atlas.csv",
      ["assay", "auroc", "ci_lo", "ci_hi", "n_test", "ci_excludes_chance"],
      rows)

# ---- FS1 pooled reference band (held-out calibrated AUROC) ----
pooled = panel["pooled"]
pci = panel["pooled_ci_95"]
n_discriminating = sum(1 for r in rows if r[5] == 1)
write("FS1_proteingym_pooled.csv",
      ["metric", "value", "ci_lo", "ci_hi", "n_variants", "n_assays",
       "n_discriminating"],
      [["auroc_calibrated", pooled["auroc_calibrated"],
        pci["auroc_calibrated"][0], pci["auroc_calibrated"][1],
        pooled["n"], len(rows), n_discriminating]])

print("A3 export_fs1_proteingym: done")
