#!/usr/bin/env python3
"""Export tidy, plot-ready CSVs for A3 figure F9 (protein uncertainty
quantification under distribution shift). Reshapes committed N8 UQ-breadth
artifacts into long CSVs under figures/plotdata/. No network/GPU. Deterministic.

Findings -> figure:
  F9 protein UQ   N8  (research/results/N8_uq/{n8_uq_breadth_results.json,
                       uq_reliability.json, N8_uq_VERDICT.json})

Panel A: rho_unc (predicted-uncertainty vs actual-error rank correlation) by
         UQ method x FLIP split, split into moderate vs extreme shift. Methods:
         BayesianRidge (baseline), Deep Ensemble, GP-RBF, MC-Dropout.
Panel B: retrospective BO acquisition (UCB vs greedy vs random) final-best per
         split; shows UCB never beats greedy (uncertainty-driven acquisition null).
"""
import json
import os
import csv

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
RES = os.path.join(ROOT, "research", "results", "N8_uq")
OUT = os.path.join(os.path.dirname(__file__), "..", "figures", "plotdata")
os.makedirs(OUT, exist_ok=True)


def write(name, header, rows):
    p = os.path.join(OUT, name)
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print("wrote", os.path.relpath(p, ROOT), f"({len(rows)} rows)")


def load(name):
    return json.load(open(os.path.join(RES, name)))


breadth = load("n8_uq_breadth_results.json")
relia = load("uq_reliability.json")  # BayesianRidge baseline, per-split

# Pretty method labels (order = bar/legend order)
METHOD_LABEL = {
    "bayes_ridge": "BayesianRidge",
    "deep_ensemble": "Deep Ensemble",
    "gp_rbf_pca": "GP-RBF",
    "mc_dropout": "MC-Dropout",
}
# 3 added methods carry per-split is_highshift; BayesianRidge per-split comes
# from uq_reliability.json (gb1.* / meltome.*). Highshift set is fixed (VERDICT).
HIGHSHIFT = {"gb1/low_vs_high", "gb1/three_vs_rest", "gb1/two_vs_rest"}
# The single hardest extrapolation split (pred R ~0.26-0.30): all methods ~0.
EXTREME = "gb1/low_vs_high"


def shift_class(split):
    if split == EXTREME:
        return "extreme"
    if split in HIGHSHIFT:
        return "moderate-shift"
    return "in-distribution"


# ---- Panel A: rho_unc by method x split ----
rowsA = []
# BayesianRidge baseline from uq_reliability.json (nested ds -> split)
for ds, splits in relia["per_split"].items():
    for sname, v in splits.items():
        split = f"{ds}/{sname}"
        lo, hi = v["rho_unc_ci"]
        rowsA.append([
            "BayesianRidge", split, shift_class(split),
            v["rho_unc"], lo, hi,
            round(v["test_pred_pearson"], 3),
        ])
# Added methods from breadth per_split
for method, splits in breadth["per_split"].items():
    label = METHOD_LABEL[method]
    for split, v in splits.items():
        lo, hi = v["rho_unc_ci"]
        rowsA.append([
            label, split, shift_class(split),
            v["rho_unc"], lo, hi,
            round(v["test_pred_pearson"], 3),
        ])
write("F9_rho_unc.csv",
      ["method", "split", "shift_class", "rho_unc", "ci_lo", "ci_hi", "pred_pearson"],
      rowsA)

# ---- Panel B: BO acquisition final-best, BayesianRidge baseline per split ----
# uq_reliability.json carries the canonical BayesianRidge BO outcomes.
rowsB = []
for ds, splits in relia["per_split"].items():
    for sname, v in splits.items():
        split = f"{ds}/{sname}"
        bo = v["bo_final_best"]
        rowsB.append([
            split, shift_class(split),
            bo["greedy"], bo["ucb"], bo["random"], bo["thompson"],
        ])
write("F9_bo_acquisition.csv",
      ["split", "shift_class", "greedy", "ucb", "random", "thompson"],
      rowsB)

# ---- aggregate summary (for caption / annotation) ----
rowsS = []
agg = breadth["aggregate"]
# BayesianRidge aggregate from VERDICT (H2 / brr_on_cache_bo) + reliability agg
rel_agg = relia["aggregate"]
rowsS.append([
    "BayesianRidge",
    rel_agg["mean_rho_unc_all"], rel_agg["mean_rho_unc_highshift"],
    0,  # ucb_beats_greedy (brr_on_cache_bo = 0)
    rel_agg["ucb_beats_random"], rel_agg["n_splits"],
])
for method, a in agg.items():
    rowsS.append([
        METHOD_LABEL[method],
        a["mean_rho_unc_all"], a["mean_rho_unc_highshift"],
        a["ucb_beats_greedy"], a["ucb_beats_random"], a["n_splits"],
    ])
write("F9_summary.csv",
      ["method", "mean_rho_unc_all", "mean_rho_unc_highshift",
       "ucb_beats_greedy", "ucb_beats_random", "n_splits"],
      rowsS)
