#!/usr/bin/env python3
"""Regenerate A3 boundary-atlas plotdata from RAW source artifacts (reproducibility loop).
CRITICAL: DMS AUROC must come from the RAW model_score_primary via audit_raw_auroc (=0.725),
NOT the superseded clamped scale_*_masked.json (0.706) — the clamp distorts rank-based AUROC
(project discipline: AUROC on raw scores, not clamped confidence).
"""
import json, csv, sys
import numpy as np
from pathlib import Path
sys.path.insert(0, "/home/zeyufu/Desktop/labs/active/reliability-program/biochemical/biochem-research/research/scaffolds/protein_ligand/stage0/featurize")
import audit_raw_auroc as ara

ROOT = "/home/zeyufu/Desktop/labs/active/reliability-program/biochemical/biochem-research"
R = f"{ROOT}/research/results"
OUT = f"{ROOT}/research/results/a3_meta"

def boot_ci(vals, seed0=0, n=2000):
    vals = np.asarray(vals, float)
    bs = [np.mean(np.random.default_rng(seed0 + i).choice(vals, len(vals))) for i in range(n)]
    return round(float(np.percentile(bs, 2.5)), 3), round(float(np.percentile(bs, 97.5)), 3)

# ---- Panel A: discrimination AUROC ----
dms = ara.evaluate(Path(f"{R}/N1_scale/expand80/manifest_650M_masked.csv"))  # RAW
dlo, dhi = boot_ci(list(dms["per_protein"].values()))
n5 = json.load(open(f"{R}/N5_arch/N5_summary.json"))
protbert = round(n5.get("protbert_protein_mean_auroc_shared69", n5.get("protbert_protein_mean_auroc")), 3)
pg = json.load(open(f"{R}/N6_clinvar/n6c_consistent.json"))["per_gene"]
ca = [v["A_clinvar_whole"] for v in pg.values() if v.get("A_clinvar_whole") is not None]
clo, chi = boot_ci(ca)
aurocA = [
    ("DMS fitness (ESM-2 650M)", round(dms["protein_mean_auroc"], 3), dlo, dhi, "n=69 proteins (raw masked-marginal)"),
    ("DMS fitness (ProtBERT)", protbert, "", "", "encoder-robust (no committed per-protein CI)"),
    ("ClinVar path/benign (per-gene)", round(float(np.mean(ca)), 3), clo, chi, "n=20 genes"),
]
with open(f"{OUT}/a3_atlas_auroc.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["task", "value", "lo", "hi", "note"]); [w.writerow(r) for r in aurocA]

# ---- Panel B: Spearman physical tasks ----
n4 = json.load(open(f"{R}/N4_antibody_ddg/results.json"))
def g(k):
    o = n4[k]; return round(o.get("rho", o.get("spearman")), 3), round(o["ci_lo"], 3), round(o["ci_hi"], 3)
n10 = json.load(open(f"{R}/N10_ddg/ddg_calibration_v2.json"))["wt"]  # v2, wt.spearman_heldout
spB = [
    ("Stability ΔΔG (MegaScale)", round(n10["spearman_heldout"], 3), "", "", "no committed CI (ddg_calibration_v2.json:wt.spearman_heldout)"),
    ("Binding ΔΔG, general (non-antibody)", *g("general_nonAB"), ""),
    ("Binding ΔΔG, overall", *g("overall"), ""),
    ("Binding ΔΔG, antibody interface", *g("antibody_interface"), ""),
    ("Binding ΔΔG, antibody–antigen", *g("antibody_ABAG"), ""),
    ("Binding ΔΔG, non-interface", *g("noninterface"), ""),
]
with open(f"{OUT}/a3_atlas_spearman.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["task", "value", "lo", "hi", "note"]); [w.writerow(r) for r in spB]

print(f"DMS-ESM2 raw {dms['protein_mean_auroc']} [{dlo},{dhi}] | ProtBERT {protbert} | ClinVar {round(float(np.mean(ca)),3)} [{clo},{chi}]")
print("wrote a3_atlas_auroc.csv, a3_atlas_spearman.csv")
