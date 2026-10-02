#!/usr/bin/env python3
"""A3-XFIT-CAL-20260918: NEW 5-fold cross-fitted clinical calibration.

This is not a restoration of the retired historical full-20 Platt pair
(0.2443 -> 0.0207 / 91%). That pair remains unrestorable; the STOP file is
history. This script fits new calibrators under a pre-specified protocol
so a calibrator is never evaluated on its own fit fold.

Protocol (locked before looking at ECE):
  seed            = 20260918
  folds           = 5, class-stratified (y=pathogenic)
  uncalibrated p  = sigmoid((s - mu_fit) / sd_fit), s = -LLR
  calibrators     = none (uncalibrated), Platt/logistic, isotonic
  ECE             = 10 equal-width bins on [0,1]; last bin includes 1.0
  bootstrap       = 2000 percentile CI on OOF (p, y) pairs; calibrators NOT refit
  six-gene check  = same protocol on ClinVar n=1093 AND on the 20-gene six-subset

Never print 0.021 / 91% as if recovered.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import random
import sys
from collections import defaultdict
from pathlib import Path

from calibrators import ece, fit_isotonic, fit_platt, reliability_points, sigmoid

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
OUT = HERE / "outputs"
SEED = 20260918
K_FOLDS = 5
N_BOOT = 2000
N_BINS = 10
ORIG6 = (
    "AICDA_HUMAN",
    "HEM3_HUMAN",
    "KCNE1_HUMAN",
    "OPSD_HUMAN",
    "OTC_HUMAN",
    "P53_HUMAN",
)
HIST_SPLIT_SEED = 20260621


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest(path: Path) -> list[dict]:
    rows = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            raw = float(row["model_score_primary"])
            label_value = int(float(row["label_value"]))
            rows.append(
                {
                    "example_id": row["example_id"],
                    "gene": row["target_id"],
                    "mut": row["ligand_id"],
                    "raw": raw,
                    "s": -raw,
                    "y": 1 - label_value,
                }
            )
    return rows


def assign_folds(ys: list[int], k: int, seed: int) -> list[int]:
    rng = random.Random(seed)
    by_class: dict[int, list[int]] = {0: [], 1: []}
    for i, y in enumerate(ys):
        by_class[int(y)].append(i)
    folds = [-1] * len(ys)
    for members in by_class.values():
        rng.shuffle(members)
        for j, idx in enumerate(members):
            folds[idx] = j % k
    if any(f < 0 for f in folds):
        raise RuntimeError("fold assignment incomplete")
    return folds


def assert_identity_separation(rows: list[dict], folds: list[int]) -> None:
    """A held-out example must have no bootstrap copy in its fit folds."""
    if len(rows) != len(folds):
        raise ValueError("row/fold length mismatch")
    seen = {}
    for row, fold in zip(rows, folds):
        key = row.get("example_id")
        if key is None or str(key).strip() == "":
            raise ValueError("stable original example_id is required")
        if key in seen and seen[key] != fold:
            raise ValueError(f"example_id {key} crosses fit/evaluation folds")
        seen[key] = fold


def assign_identity_folds(rows: list[dict], k: int, seed: int) -> list[int]:
    """Stratify original identities; keep all resampled copies in one fold.

    For unique rows this is exactly the original assignment. Bootstrap copies
    retain their multiplicity in fitting/scoring, but cannot leak across folds.
    This groups copies of a variant, not genes; it is still in-panel calibration.
    """
    identities, labels, signatures = [], [], {}
    for row in rows:
        key = row.get("example_id")
        if key is None or str(key).strip() == "":
            raise ValueError("stable original example_id is required")
        signature = (row["gene"], row["s"], row["y"])
        if key in signatures and signatures[key] != signature:
            raise ValueError(f"conflicting data for example_id {key}")
        if key not in signatures:
            identities.append(key)
            labels.append(row["y"])
            signatures[key] = signature
    if len(identities) < k:
        raise ValueError("fewer independent example identities than folds")
    by_id = dict(zip(identities, assign_folds(labels, k, seed)))
    folds = [by_id[row["example_id"]] for row in rows]
    if set(folds) != set(range(k)):
        raise ValueError("empty evaluation fold after identity grouping")
    assert_identity_separation(rows, folds)
    return folds


def zscore_sigmoid(scores: list[float], mu: float, sd: float) -> list[float]:
    if sd == 0.0:
        sd = 1.0
    return [sigmoid((x - mu) / sd) for x in scores]


def mean_sd(values: list[float]) -> tuple[float, float]:
    mu = sum(values) / len(values)
    var = sum((x - mu) ** 2 for x in values) / len(values)
    return mu, var ** 0.5


def percentile(values: list[float], q: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    pos = (len(ordered) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(ordered) - 1)
    w = pos - lo
    return ordered[lo] * (1.0 - w) + ordered[hi] * w


def boot_ece(probs: list[float], labels: list[int], seed: int, n_boot: int) -> dict:
    rng = random.Random(seed)
    n = len(probs)
    point = ece(probs, labels, N_BINS)
    draws = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        draws.append(ece([probs[i] for i in idx], [labels[i] for i in idx], N_BINS))
    return {
        "ece": point,
        "ci_lo": percentile(draws, 0.025),
        "ci_hi": percentile(draws, 0.975),
        "n_boot": n_boot,
        "n": n,
        "seed": seed,
    }


def crossfit(rows: list[dict], seed: int = SEED) -> dict:
    ys = [r["y"] for r in rows]
    folds = assign_identity_folds(rows, K_FOLDS, seed)
    fold_counts = []
    oof = {
        "uncalibrated": [0.0] * len(rows),
        "platt": [0.0] * len(rows),
        "isotonic": [0.0] * len(rows),
    }
    platt_coeff = []
    for k in range(K_FOLDS):
        fit_idx = [i for i, f in enumerate(folds) if f != k]
        test_idx = [i for i, f in enumerate(folds) if f == k]
        s_fit = [rows[i]["s"] for i in fit_idx]
        y_fit = [rows[i]["y"] for i in fit_idx]
        s_test = [rows[i]["s"] for i in test_idx]
        mu, sd = mean_sd(s_fit)
        p_fit = zscore_sigmoid(s_fit, mu, sd)
        p_test = zscore_sigmoid(s_test, mu, sd)
        if len(set(y_fit)) < 2:
            raise RuntimeError(f"fit fold {k} is single-class")
        platt = fit_platt(p_fit, y_fit)
        iso = fit_isotonic(p_fit, y_fit)
        q_platt = platt(p_test)
        q_iso = iso(p_test)
        for j, i in enumerate(test_idx):
            oof["uncalibrated"][i] = p_test[j]
            oof["platt"][i] = q_platt[j]
            oof["isotonic"][i] = q_iso[j]
        fold_counts.append(
            {
                "fold": k,
                "n_fit": len(fit_idx),
                "n_eval": len(test_idx),
                "n_eval_pathogenic": sum(rows[i]["y"] for i in test_idx),
                "n_eval_benign": sum(1 - rows[i]["y"] for i in test_idx),
                "n_fit_pathogenic": sum(y_fit),
                "n_fit_benign": len(y_fit) - sum(y_fit),
                "z_mu_fit": mu,
                "z_sd_fit": sd,
                "platt_a": platt.coeff[0],  # type: ignore[attr-defined]
                "platt_b": platt.coeff[1],  # type: ignore[attr-defined]
            }
        )
        platt_coeff.append({"fold": k, "a": platt.coeff[0], "b": platt.coeff[1]})  # type: ignore[attr-defined]
    labels = ys
    summary = {}
    bins = []
    for name, probs in oof.items():
        rec = boot_ece(probs, labels, seed=seed, n_boot=N_BOOT)
        rec["calibrator"] = name
        rec["degrades_vs_uncalibrated"] = None
        summary[name] = rec
        for pt in reliability_points(probs, labels, N_BINS):
            row = {"calibrator": name, **pt}
            bins.append(row)
    raw_ece = summary["uncalibrated"]["ece"]
    for name in ("platt", "isotonic"):
        summary[name]["degrades_vs_uncalibrated"] = summary[name]["ece"] > raw_ece
        if raw_ece > 0:
            summary[name]["delta_vs_uncalibrated"] = summary[name]["ece"] - raw_ece
            summary[name]["reduction_pct"] = 100.0 * (raw_ece - summary[name]["ece"]) / raw_ece
        else:
            summary[name]["delta_vs_uncalibrated"] = 0.0
            summary[name]["reduction_pct"] = float("nan")
    return {
        "n": len(rows),
        "n_pathogenic": sum(ys),
        "n_benign": len(ys) - sum(ys),
        "n_genes": len({r["gene"] for r in rows}),
        "folds": fold_counts,
        "ece": summary,
        "bins": bins,
        "oof": oof,
        "labels": labels,
        "genes": [r["gene"] for r in rows],
        "platt_coeff": platt_coeff,
        "seed": seed,
        "k_folds": K_FOLDS,
    }


def holdout_50_50(rows: list[dict], split_seed: int, z_on: str) -> dict:
    rng = random.Random(split_seed)
    s = [r["s"] for r in rows]
    y = [r["y"] for r in rows]
    pos = [i for i, yi in enumerate(y) if yi == 1]
    neg = [i for i, yi in enumerate(y) if yi == 0]
    rng.shuffle(pos)
    rng.shuffle(neg)
    fit_idx = set(pos[: len(pos) // 2] + neg[: len(neg) // 2])
    test_idx = [i for i in range(len(rows)) if i not in fit_idx]
    if z_on == "all":
        mu, sd = mean_sd(s)
    elif z_on == "fit":
        mu, sd = mean_sd([s[i] for i in fit_idx])
    else:
        raise ValueError(z_on)
    raw_prob = zscore_sigmoid(s, mu, sd)
    p_fit = [raw_prob[i] for i in fit_idx]
    y_fit = [y[i] for i in fit_idx]
    p_test = [raw_prob[i] for i in test_idx]
    y_test = [y[i] for i in test_idx]
    platt = fit_platt(p_fit, y_fit)
    iso = fit_isotonic(p_fit, y_fit)
    q_platt = platt(p_test)
    q_iso = iso(p_test)
    return {
        "split_seed": split_seed,
        "z_on": z_on,
        "n_fit": len(fit_idx),
        "n_test": len(test_idx),
        "ece_uncalibrated": ece(p_test, y_test, N_BINS),
        "ece_platt": ece(q_platt, y_test, N_BINS),
        "ece_isotonic": ece(q_iso, y_test, N_BINS),
        "bins_uncalibrated": reliability_points(p_test, y_test, N_BINS),
        "bins_platt": reliability_points(q_platt, y_test, N_BINS),
        "bins_isotonic": reliability_points(q_iso, y_test, N_BINS),
    }


def per_group_ece(result: dict, key_fn) -> list[dict]:
    labels = result["labels"]
    keys = [key_fn(i) for i in range(result["n"])]
    groups: dict[str, list[int]] = defaultdict(list)
    for i, key in enumerate(keys):
        groups[key].append(i)
    out = []
    for key, idx in sorted(groups.items()):
        rec = {
            "group": key,
            "n": len(idx),
            "n_pathogenic": sum(labels[i] for i in idx),
            "n_benign": sum(1 - labels[i] for i in idx),
        }
        for name, probs in result["oof"].items():
            p = [probs[i] for i in idx]
            y = [labels[i] for i in idx]
            if len(set(y)) < 2 or len(y) < 8:
                rec[f"ece_{name}"] = None
                rec[f"ece_{name}_note"] = "skipped (n<8 or single-class)"
                continue
            rec[f"ece_{name}"] = ece(p, y, N_BINS)
            if len(y) >= 30:
                boot = boot_ece(p, y, seed=SEED, n_boot=N_BOOT)
                rec[f"ece_{name}_ci_lo"] = boot["ci_lo"]
                rec[f"ece_{name}_ci_hi"] = boot["ci_hi"]
        out.append(rec)
    return out


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def round_rec(rec: dict, digits: int = 6) -> dict:
    out = {}
    for key, value in rec.items():
        if isinstance(value, float):
            out[key] = None if math.isnan(value) else round(value, digits)
        else:
            out[key] = value
    return out


def public_ece(rec: dict) -> dict:
    return {
        "calibrator": rec["calibrator"],
        "ece": rec["ece"],
        "ci_lo": rec["ci_lo"],
        "ci_hi": rec["ci_hi"],
        "n": rec["n"],
        "n_boot": rec["n_boot"],
        "seed": rec["seed"],
        "degrades_vs_uncalibrated": rec.get("degrades_vs_uncalibrated"),
        "delta_vs_uncalibrated": rec.get("delta_vs_uncalibrated"),
        "reduction_pct": rec.get("reduction_pct"),
    }


def run_panel(name: str, rows: list[dict], source: dict) -> dict:
    result = crossfit(rows, seed=SEED)
    gene_rows = per_group_ece(result, lambda i: result["genes"][i])
    family_rows = per_group_ece(
        result, lambda i: "ORIG6" if result["genes"][i] in ORIG6 else "OTHER14"
    )
    # original writer used seed=1 for the pooled 50/50 split
    hist_like_all = holdout_50_50(rows, split_seed=1, z_on="all")
    hist_like_fit = holdout_50_50(rows, split_seed=1, z_on="fit")
    panel = {
        "panel": name,
        "protocol": "A3-XFIT-CAL-20260918",
        "not_a_restoration": True,
        "source": source,
        "n": result["n"],
        "n_pathogenic": result["n_pathogenic"],
        "n_benign": result["n_benign"],
        "n_genes": result["n_genes"],
        "k_folds": K_FOLDS,
        "seed": SEED,
        "n_boot": N_BOOT,
        "n_bins": N_BINS,
        "bin_definition": "equal-width 10 bins on [0,1]; bin b=[b/10,(b+1)/10) for b=0..8; bin 9=[0.9,1.0]",
        "bootstrap": "percentile 95% CI by resampling OOF (p,y) pairs; calibrators not refit",
        "fold_counts": result["folds"],
        "ece_xfit": {k: public_ece(v) for k, v in result["ece"].items()},
        "diagnostic_50_50_z_on_all_split_seed_1": {
            "n_fit": hist_like_all["n_fit"],
            "n_test": hist_like_all["n_test"],
            "ece_uncalibrated": hist_like_all["ece_uncalibrated"],
            "ece_platt": hist_like_all["ece_platt"],
            "ece_isotonic": hist_like_all["ece_isotonic"],
            "note": "diagnostic of the historical 50/50 writer, not the printed xfit claim",
        },
        "diagnostic_50_50_z_on_fit_split_seed_1": {
            "n_fit": hist_like_fit["n_fit"],
            "n_test": hist_like_fit["n_test"],
            "ece_uncalibrated": hist_like_fit["ece_uncalibrated"],
            "ece_platt": hist_like_fit["ece_platt"],
            "ece_isotonic": hist_like_fit["ece_isotonic"],
        },
        "per_gene": gene_rows,
        "per_family": family_rows,
        "bins": result["bins"],
    }
    return panel


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gof_path = INPUTS / "manifest_gof_650M.csv"
    clin_path = INPUTS / "manifest_clinical_650M.csv"
    gof = load_manifest(gof_path)
    clin = load_manifest(clin_path)
    if len(gof) != 4677:
        raise SystemExit(f"expected n=4677 gof rows, got {len(gof)}")
    if len(clin) != 1093:
        raise SystemExit(f"expected n=1093 clinical rows, got {len(clin)}")
    gof6 = [r for r in gof if r["gene"] in ORIG6]
    provenance = {
        "protocol_id": "A3-XFIT-CAL-20260918",
        "title": "NEW 5-fold cross-fitted calibration (not a historical restoration)",
        "forbidden_claim": "Do not treat any number here as recovery of 0.0207 or 91%",
        "python": sys.version,
        "platform": platform.platform(),
        "seed": SEED,
        "k_folds": K_FOLDS,
        "n_boot": N_BOOT,
        "inputs": {
            "manifest_gof_650M.csv": {"sha256": sha256(gof_path), "n": len(gof)},
            "manifest_clinical_650M.csv": {"sha256": sha256(clin_path), "n": len(clin)},
        },
    }
    panels = {
        "full20_gof": run_panel(
            "full20_gof",
            gof,
            {"file": "manifest_gof_650M.csv", "label_source": "DMS_score_bin / GOF-LoF 20-gene panel"},
        ),
        "six_from_gof": run_panel(
            "six_from_gof",
            gof6,
            {"file": "manifest_gof_650M.csv", "subset": list(ORIG6), "label_source": "DMS_score_bin six-gene slice"},
        ),
        "six_clinvar": run_panel(
            "six_clinvar",
            clin,
            {"file": "manifest_clinical_650M.csv", "label_source": "ClinVar pathogenic/benign (printed 0.157->0.047 source)"},
        ),
    }
    summary = {"provenance": provenance, "panels": {}}
    for name, panel in panels.items():
        slim = {k: v for k, v in panel.items() if k != "bins"}
        summary["panels"][name] = slim
        write_csv(
            OUT / f"{name}_reliability_bins.csv",
            [round_rec(r, 6) for r in panel["bins"]],
            [
                "calibrator",
                "bin",
                "bin_lo",
                "bin_hi",
                "bin_mid",
                "mean_pred",
                "obs_rate",
                "n",
            ],
        )
        write_csv(
            OUT / f"{name}_per_gene_ece.csv",
            [round_rec(r, 6) for r in panel["per_gene"]],
            [
                "group",
                "n",
                "n_pathogenic",
                "n_benign",
                "ece_uncalibrated",
                "ece_uncalibrated_ci_lo",
                "ece_uncalibrated_ci_hi",
                "ece_platt",
                "ece_platt_ci_lo",
                "ece_platt_ci_hi",
                "ece_isotonic",
                "ece_isotonic_ci_lo",
                "ece_isotonic_ci_hi",
                "ece_uncalibrated_note",
                "ece_platt_note",
                "ece_isotonic_note",
            ],
        )
        write_csv(
            OUT / f"{name}_fold_counts.csv",
            [round_rec(r, 6) for r in panel["fold_counts"]],
            [
                "fold",
                "n_fit",
                "n_eval",
                "n_eval_pathogenic",
                "n_eval_benign",
                "n_fit_pathogenic",
                "n_fit_benign",
                "z_mu_fit",
                "z_sd_fit",
                "platt_a",
                "platt_b",
            ],
        )
    ece_rows = []
    for name, panel in panels.items():
        for cal, rec in panel["ece_xfit"].items():
            ece_rows.append(
                round_rec(
                    {
                        "panel": name,
                        "protocol": "5-fold-xfit",
                        **rec,
                    },
                    6,
                )
            )
    write_csv(
        OUT / "ece_xfit_summary.csv",
        ece_rows,
        [
            "panel",
            "protocol",
            "calibrator",
            "ece",
            "ci_lo",
            "ci_hi",
            "n",
            "n_boot",
            "seed",
            "degrades_vs_uncalibrated",
            "delta_vs_uncalibrated",
            "reduction_pct",
        ],
    )
    (OUT / "xfit_calibration.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("=== A3-XFIT-CAL-20260918 ===")
    for name, panel in panels.items():
        print(f"\n[{name}] n={panel['n']} genes={panel['n_genes']}")
        for cal, rec in panel["ece_xfit"].items():
            flag = " DEGRADES" if rec.get("degrades_vs_uncalibrated") else ""
            print(
                f"  {cal:13s} ECE={rec['ece']:.4f}  95% CI [{rec['ci_lo']:.4f}, {rec['ci_hi']:.4f}]  n={rec['n']}{flag}"
            )
        d = panel["diagnostic_50_50_z_on_all_split_seed_1"]
        print(
            f"  diagnostic 50/50 z-on-all seed=1: raw={d['ece_uncalibrated']:.4f} "
            f"platt={d['ece_platt']:.4f} iso={d['ece_isotonic']:.4f} n_test={d['n_test']}"
        )


if __name__ == "__main__":
    main()
