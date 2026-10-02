#!/usr/bin/env python3
"""Export FS2 (clinical-axis controls) plot data from the CANONICAL N6c consistent
re-score — NOT the withdrawn N6/N6b artifacts (gof_boundary.json / gnomad_control.json
used the long-protein bf16 scoring that N6c corrected: e.g. MET 0.37->0.80, BRAF
0.39->0.83). Everything here is read off research/results/N6_clinvar/n6c_consistent.json,
where ClinVar-benign and gnomAD-AC>=2-benign were scored against the same pathogenic
set in one deterministic fp32 pass.
"""
import json, os, csv

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
OUT = os.path.join(os.path.dirname(__file__), "..", "figures", "plotdata")
os.makedirs(OUT, exist_ok=True)
d = json.load(open(os.path.join(ROOT, "research", "results", "N6_clinvar", "n6c_consistent.json")))["per_gene"]


def cls(key):
    if key.endswith("_GOF"):
        return "gain-of-function"
    if key.endswith("_LOF"):
        return "loss-of-function"
    return "other"


def name(key):
    return key.replace("_HUMAN", "").replace("_GOF", "").replace("_LOF", "")


rows = []
for g, v in d.items():
    win = v.get("windowed")
    ac = v.get("A_clinvar_win") if win and v.get("A_clinvar_win") is not None else v.get("A_clinvar_whole")
    ag = v.get("A_gnomad_win") if win and v.get("A_gnomad_win") is not None else v.get("A_gnomad_whole")
    rows.append([name(g), cls(g), ac, ag, ag - ac, v.get("n_path"), v.get("wt_len")])

rows.sort(key=lambda r: r[2])
with open(os.path.join(OUT, "FS2_clinvar_vs_gnomad.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["gene", "class", "auroc_clinvar_benign", "auroc_gnomad_benign", "delta_gnomad_minus_clinvar", "n_path", "wt_len"])
    w.writerows(rows)

# class means (N6c finds no LoF/GoF structure)
cm = {}
for r in rows:
    cm.setdefault(r[1], []).append(r[2])
with open(os.path.join(OUT, "FS2_class_means.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["class", "mean_clinvar_auroc", "n_genes"])
    for c, vals in cm.items():
        w.writerow([c, round(sum(vals) / len(vals), 4), len(vals)])

print("wrote FS2_clinvar_vs_gnomad.csv (%d genes) and FS2_class_means.csv" % len(rows))
print("class means:", {c: round(sum(v) / len(v), 3) for c, v in cm.items()})
