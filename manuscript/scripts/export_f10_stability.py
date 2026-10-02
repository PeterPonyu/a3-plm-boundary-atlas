#!/usr/bin/env python3
"""Export tidy, plot-ready CSVs for A3 figure F10 (protein stability ddG
zero-shot: a controlled negative). Reshapes committed N10 ddG calibration
artifacts into long CSVs under figures/plotdata/. No network/GPU. Deterministic.

Findings -> figure:
  F10 stability ddG  N10  (research/results/N10_ddg/{ddg_calibration.json,
                           ddg_calibration_v2.json})

Panel A: discrimination Spearman train (in-distribution) vs homology-held-out --
         they transfer ~1:1 (~0.38 vs ~0.34), a pre-registered controlled
         negative (no reliability gap). Dashed reference line at the supervised
         ThermoMPNN literature level (~0.65) shows ESM zero-shot is WEAK but
         stable.
Panel B: calibration RMSE train vs homology-held-out (~equal, inflation ~+0.04)
         = no homology calibration gap.
"""
import json
import os
import csv
import argparse

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--results-dir",
    default=os.path.join(ROOT, "research", "results", "N10_ddg"),
    help="Directory containing ddg_calibration.json and ddg_calibration_v2.json",
)
parser.add_argument(
    "--output-dir",
    default=os.path.join(os.path.dirname(__file__), "..", "figures", "plotdata"),
)
args = parser.parse_args()
RES = os.path.abspath(args.results_dir)
OUT = os.path.abspath(args.output_dir)
os.makedirs(OUT, exist_ok=True)


def write(name, header, rows):
    p = os.path.join(OUT, name)
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    print("wrote", os.path.relpath(p, ROOT), f"({len(rows)} rows)")


def load(name):
    return json.load(open(os.path.join(RES, name)))


v1 = load("ddg_calibration.json")
v2 = load("ddg_calibration_v2.json")  # wt + masked variants, cleaner held-out

# Supervised reference: ThermoMPNN reports per-structure Spearman ~0.65 on the
# Megascale/cDNA-display benchmark (literature). Drawn as a dashed line so the
# reader sees ESM zero-shot is well below the supervised ceiling though stable.
SUPERVISED_REF = 0.65

# ---- Panels A--C: clean v2 split only ----
# v1 is not a cluster-disjoint estimand: four clusters overlap train and
# held-out. Keep it in a separate provenance audit, never in headline panels.
rowsA = [
    ["v2 WT", "train", v2["wt"]["spearman_train"]],
    ["v2 WT", "homology-held-out", v2["wt"]["spearman_heldout"]],
    ["v2 masked", "train", v2["masked"]["spearman_train"]],
    ["v2 masked", "homology-held-out", v2["masked"]["spearman_heldout"]],
]
write("F10_discrimination.csv", ["variant", "split", "spearman"], rowsA)
write("F10_reference.csv", ["label", "spearman"],
      [["supervised reference (ThermoMPNN, literature)", SUPERVISED_REF]])

# ---- Panel B: calibration RMSE, train (held-in) vs homology-held-out ----
rowsB = [
    ["v2 WT", "train", v2["wt"]["calib_rmse_heldin"],
     v2["wt"]["calib_rmse_inflation"]],
    ["v2 WT", "homology-held-out", v2["wt"]["calib_rmse_homology_heldout"],
     v2["wt"]["calib_rmse_inflation"]],
    ["v2 masked", "train", v2["masked"]["calib_rmse_heldin"],
     v2["masked"]["calib_rmse_inflation"]],
    ["v2 masked", "homology-held-out", v2["masked"]["calib_rmse_homology_heldout"],
     v2["masked"]["calib_rmse_inflation"]],
]
write("F10_calibration.csv", ["variant", "split", "calib_rmse", "inflation"], rowsB)

write("F10_v1_overlap_audit.csv",
      ["field", "value"],
      [
          ["status", "superseded; excluded from headline panels"],
          ["n_train", v1["n_train"]],
          ["n_heldout", v1["n_heldout"]],
          ["clusters_train", v1["clusters_train"]],
          ["clusters_heldout", v1["clusters_heldout"]],
          ["cluster_overlap", v1["cluster_overlap"]],
          ["discrimination_spearman_train", v1["discrimination_spearman_train"]],
          ["discrimination_spearman_heldout", v1["discrimination_spearman_heldout"]],
          ["calib_rmse_train_heldin", v1["calib_rmse_train_heldin"]],
          ["calib_rmse_homology_heldout", v1["calib_rmse_homology_heldout"]],
          ["calib_rmse_inflation", v1["calib_rmse_inflation"]],
      ])

# ---- meta (cluster counts, for caption) ----
write("F10_meta.csv",
      ["field", "value"],
      [
          ["n_train", v2["n_train"]],
          ["n_heldout", v2["n_heldout_clean"]],
          ["clusters_train", v2["train_clusters"]],
          ["clusters_heldout", v2["heldout_clusters_clean"]],
          ["cluster_overlap", 0],
      ])
