#!/usr/bin/env python3
"""Export tidy, plot-ready CSVs for the A3 (PLM mutation-effect boundary atlas)
R/ggplot2 figures. Reshapes committed result artifacts into long CSVs under
figures/plotdata/. No network/GPU. Deterministic.

Findings -> figures:
  F1 scale plateau            N1  (research/results/N1_scale/expand80/scale_*_masked.json, n=69)
  F2 cross-MLM-encoder        N5  (research/results/N5_arch/N5_summary.json)
  F3 clinical per-gene forest N6c (research/results/N6_clinvar/n6c_consistent.json)
  F4 binding boundary         N4  (research/results/N4_antibody_ddg/results.json)
  F5 developability map       N9  (research/results/N9_antibody/developability_corrected.json)
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


def load(*parts):
    return json.load(open(os.path.join(RES, *parts)))


# ---- F1 scale plateau (N1, expand80 masked-marginal, n=69) ----
# CANONICAL raw-score AUROC (audit_raw_auroc.py), NOT the scale_*.json values:
# those were computed on the clamped model_confidence_primary, the very post-hoc
# transform the N1 correction flags as distorting rank-based AUROC. Reading off
# the raw score makes 650M (0.725) identical to the N5/Fig 2 value, resolving the
# Fig1/Fig2 inconsistency. 3B = 650M = 0.725 exactly (no gain at 3B).
raw = load("N1_scale", "expand80", "raw_auroc_by_tier.json")
rows = []
for i, lab in enumerate(["35M", "150M", "650M", "3B"]):
    t = raw["tiers"][lab]
    rows.append([lab, i + 1, t["protein_mean_auroc"], raw["n_proteins"], t["n_discriminating"]])
write("F1_scale.csv", ["model_size", "order", "protein_mean_auroc", "n_proteins", "n_discriminating"], rows)

# ---- F2 cross-MLM-encoder agreement (N5) ----
d = load("N5_arch", "N5_summary.json")
rows = [["ProtBERT", "BERT/UniRef100+BFD", d["protbert_protein_mean_auroc"]],
        ["ESM-2 650M (masked)", "ESM/UR50", d["esm2_650M_masked_protein_mean_auroc"]],
        ["ESM-2 650M (WT)", "ESM/UR50", d["esm2_650M_wt_protein_mean_auroc"]]]
write("F2_crossmodel.csv", ["model", "family", "protein_mean_auroc"], rows)
# the decisive cross-model agreement statistic (per-protein Spearman), for the caption/annotation
write("F2_agreement.csv", ["comparison", "spearman"],
      [["ProtBERT vs ESM-2 (masked)", d["spearman_per_protein_vs_esm2_masked"]],
       ["ProtBERT vs ESM-2 (WT)", d["spearman_per_protein_vs_esm2_wt"]]])

# ---- F3 clinical per-gene forest (N6c) ----
d = load("N6_clinvar", "n6c_consistent.json")["per_gene"]
rows = []
for gene, v in d.items():
    auroc = v["A_clinvar_win"] if v.get("windowed") and v.get("A_clinvar_win") is not None else v["A_clinvar_whole"]
    rows.append([gene.replace("_HUMAN", "").replace("_GOF", "").replace("_LOF", ""),
                 auroc, v["n_path"], v["wt_len"], v.get("A_gnomad_win") if v.get("windowed") else v.get("A_gnomad_whole")])
rows.sort(key=lambda r: r[1])  # ascending AUROC for a clean forest
write("F3_clinvar.csv", ["gene", "auroc_clinvar", "n_path", "wt_len", "auroc_gnomad"], rows)

# ---- F4 binding boundary by subset (N4) ----
d = load("N4_antibody_ddg", "results.json")
SUB = [("overall", "overall"), ("general (non-Ab)", "general_nonAB"),
       ("antibody (Ab-Ag)", "antibody_ABAG"), ("interface", "interface"),
       ("non-interface", "noninterface")]
rows = []
for lab, key in SUB:
    v = d[key]
    rows.append([lab, key, v["rho"], v["ci_lo"], v["ci_hi"], v.get("n")])
write("F4_binding.csv", ["subset", "key", "spearman", "ci_lo", "ci_hi", "n"], rows)

# ---- F5 developability map (N9 corrected): per-property significant fraction ----
pp = load("N9_antibody", "developability_corrected.json")["summary"]["per_property"]
PROP_ORDER = ["polyreactivity", "aggregation", "pharmacokinetics", "thermostability",
              "immunogenicity", "expression", "binding"]
rows = []
for prop in PROP_ORDER:
    v = pp.get(prop)
    if not v:
        continue
    n = v["n_datasets"]; sig = v["sig_raw_bh"]
    rows.append([prop, n, sig, round(sig / n, 4) if n else 0.0, v["n_applicable_multigermline"]])
write("F5_developability.csv", ["property", "n_datasets", "n_sig_raw_bh", "frac_sig", "n_applicable_multigermline"], rows)

print("A3 export_plot_data: done")
