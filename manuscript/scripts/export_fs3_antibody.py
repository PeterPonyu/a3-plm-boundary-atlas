#!/usr/bin/env python3
"""Export plot-ready CSV for A3 supplementary figure FS3 (antibody binding-boundary
robustness). Reshapes the committed N4 antibody-ddG artifacts into a long forest CSV
under figures/plotdata/. No network/GPU. Deterministic.

Sources (research/results/N4_antibody_ddg/):
  results.json        ESM-2 650M, WT-marginal  (rho/ci_lo/ci_hi per subset)
  results_3b_wt.json  ESM-2 3B,   WT-marginal  (subsets_3b: spearman_llr_negddg + ci [lo,hi])
  masked_results.json ESM-2 650M, masked-marginal (masked_marginal: rho/ci_lo/ci_hi per subset)

Key messages preserved verbatim from the JSONs:
  - 3B does NOT rescue: 3B overall rho=0.148 < 650M overall rho=0.197.
  - masked approx wt: scoring method is robust (per-subset rho near-identical).
  The binding boundary (weak overall, ~null antibody/non-interface) is not fixed by
  scaling to 3B nor by switching WT-marginal -> masked-marginal scoring.
"""
import json, os, csv

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
RES = os.path.join(ROOT, "research", "results", "N4_antibody_ddg")
OUT = os.path.join(os.path.dirname(__file__), "..", "figures", "plotdata")
os.makedirs(OUT, exist_ok=True)


def load(fname):
    return json.load(open(os.path.join(RES, fname)))


def write(name, header, rows):
    p = os.path.join(OUT, name)
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    print("wrote", os.path.relpath(p, ROOT), f"({len(rows)} rows)")


r650 = load("results.json")            # 650M wt: keys overall/antibody_ABAG/general_nonAB/interface/noninterface
r3b = load("results_3b_wt.json")       # 3B wt:  subsets_3b.{overall,antibody,general_non_ab,interface_COR_SUP_RIM,non_interface_INT_SUR}
rmsk = load("masked_results.json")     # masked: masked_marginal.{overall,antibody_ab,general_nonab,interface,noninterface}

# Subset display order + the matching key in each of the three schemas.
# (display_label, key_650m, key_3b, key_masked)
SUBSETS = [
    ("overall",        "overall",       "overall",                "overall"),
    ("general",        "general_nonAB", "general_non_ab",         "general_nonab"),
    ("antibody",       "antibody_ABAG", "antibody",               "antibody_ab"),
    ("interface",      "interface",     "interface_COR_SUP_RIM",  "interface"),
    ("non-interface",  "noninterface",  "non_interface_INT_SUR",  "noninterface"),
]

rows = []
for i, (label, k650, k3b, kmsk) in enumerate(SUBSETS):
    # ESM-2 650M (WT-marginal)
    s = r650[k650]
    rows.append([i, label, "ESM-2 650M (WT)", s["n"], s["rho"], s["ci_lo"], s["ci_hi"]])
    # ESM-2 3B (WT-marginal) — ci is [lo, hi]
    s = r3b["subsets_3b"][k3b]
    ci = s["ci"]
    rows.append([i, label, "ESM-2 3B (WT)", s["n"], s["spearman_llr_negddg"], ci[0], ci[1]])
    # ESM-2 650M masked-marginal
    s = rmsk["masked_marginal"][kmsk]
    rows.append([i, label, "ESM-2 650M (masked)", s["n"], s["rho"], s["ci_lo"], s["ci_hi"]])

write("FS3_antibody_robustness.csv",
      ["order", "subset", "method", "n", "rho", "ci_lo", "ci_hi"], rows)

# Echo the headline contrasts for the caption / validation.
print("CHECK overall: 650M=%.3f  3B=%.3f  masked=%.3f" % (
    r650["overall"]["rho"],
    r3b["subsets_3b"]["overall"]["spearman_llr_negddg"],
    rmsk["masked_marginal"]["overall"]["rho"]))
print("CHECK 3B-vs-650M overall delta = %.3f (3B does not rescue)" %
      r3b["delta_overall_spearman_3b_minus_650m_overlap"])
print("CHECK wt-vs-masked corr: pearson=%.3f spearman=%.3f" % (
      rmsk["wt_vs_masked_corr"]["pearson"], rmsk["wt_vs_masked_corr"]["spearman"]))
