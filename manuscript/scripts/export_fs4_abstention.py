#!/usr/bin/env python3
"""Export plot-ready CSVs for A3 supplementary figure FS4 (confidence is anti-helpful
for selective abstention across ESM-2 tiers). Reshapes the committed BIO abstention
artifacts into long CSVs under figures/plotdata/. No network/GPU. Deterministic.

Sources (research/results/BIO_abstention/):
  bio_tiers_result.json    per ESM tier (35M/150M/650M/ESM-1v): auroc_full,
                           auroc_keephi50, auroc_shuffle50, confidence_anti_helpful
  bio_abstention_result.json  650M coverage-risk curve (protein_mean_auroc vs
                           shuffled_conf_auroc at coverage 1.0..0.5) for the inset panel

Key message preserved verbatim from the JSONs:
  Keeping the highest-confidence 50% of variants DROPS protein-mean AUROC below both
  the full set and the shuffled-confidence control in 4/4 tiers (keep-hi approx<= shuffle):
  ESM-2's clamped confidence does NOT select the reliable predictions. Honest negative.
"""
import json, os, csv

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
RES = os.path.join(ROOT, "research", "results", "BIO_abstention")
OUT = os.path.join(os.path.dirname(__file__), "..", "figures", "plotdata")
os.makedirs(OUT, exist_ok=True)


def load(fname):
    return json.load(open(os.path.join(RES, fname)))


def write(name, header, rows):
    p = os.path.join(OUT, name)
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    print("wrote", os.path.relpath(p, ROOT), f"({len(rows)} rows)")


# ---- Panel A: per-tier full vs keep-hi-50% vs shuffle-50% ----
tiers = load("bio_tiers_result.json")
DISPLAY = {"35M": "ESM-2 35M", "150M": "ESM-2 150M", "650M": "ESM-2 650M", "1v": "ESM-1v"}
COND = [("auroc_full", "full (100%)"),
        ("auroc_keephi50", "keep-hi 50%"),
        ("auroc_shuffle50", "shuffle 50%")]
rows = []
for i, r in enumerate(tiers["results"]):
    for key, cond in COND:
        rows.append([i, DISPLAY[r["tier"]], cond, r[key], r["confidence_anti_helpful"]])
write("FS4_abstention_tiers.csv",
      ["order", "tier", "condition", "auroc", "anti_helpful"], rows)

# ---- Panel B: 650M coverage-risk curve (confidence vs shuffled control) ----
cov = load("bio_abstention_result.json")
crows = []
for d in cov["coverage_risk"]:
    crows.append([d["coverage"], d["protein_mean_auroc"], d["shuffled_conf_auroc"], d["n_proteins"]])
write("FS4_abstention_coverage.csv",
      ["coverage", "conf_auroc", "shuffle_auroc", "n_proteins"], crows)

# Validation echoes for the caption.
print("CHECK n_anti_helpful = %d / %d tiers" % (
      tiers["n_anti_helpful"], len(tiers["results"])))
for r in tiers["results"]:
    print("  %-5s full=%.4f keep-hi-50=%.4f shuffle-50=%.4f anti=%s" % (
        r["tier"], r["auroc_full"], r["auroc_keephi50"],
        r["auroc_shuffle50"], r["confidence_anti_helpful"]))
print("CHECK 650M coverage 1.0 conf=%.4f matches base_auroc=%.4f" % (
      cov["coverage_risk"][0]["protein_mean_auroc"], cov["base_auroc"]))
