#!/usr/bin/env python3
"""A3-XFIT-CAL-R3-20260919: per-gene per-bin ECE and nested xfit bootstrap.

NEW protocol on the locked xfit object. This is not a restoration of the
historical full-20 Platt pair (0.2443 -> 0.0207 / 91%). The historical STOP
file manuscript/figures/plotdata/F7_20gene_perbin_STATUS.txt is unchanged
and is not the authority for these bins.

Point-estimate OOF maps reuse A3-XFIT-CAL-20260918 / R2 writers.
Nested bootstrap resamples original identities and re-fits 5-fold calibrators
(B=200). All copies of one example stay in the same fold. The previous
row-copy split was leaky; its archived nested intervals are withdrawn.
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

import numpy as np

from calibrators import ece, fit_isotonic, reliability_points, sigmoid
import run_xfit_calibration as xfit
import run_xfit_round2 as r2

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs" / "round3_identity_grouped_candidate"
STOP = (
    HERE.parents[1]
    / "manuscript"
    / "figures"
    / "plotdata"
    / "F7_20gene_perbin_STATUS.txt"
)
PROTOCOL = "A3-XFIT-NESTED-IDENTITY-GROUPED-CANDIDATE"
PARENT = "A3-XFIT-CAL-20260918"
EXPECTED_UNCAL = 0.24848355100110145
N_NESTED = 200
NESTED_SEED = 20260919
N_BINS = xfit.N_BINS
MAPS = ("uncalibrated", "pooled_platt", "pooled_isotonic", "pergene_platt")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def r6(value):
    if value is None:
        return None
    if isinstance(value, float):
        if math.isnan(value):
            return None
        return round(value, 6)
    return value


def logit_clip(p: float) -> float:
    p = min(1.0 - 1e-6, max(1e-6, float(p)))
    return math.log(p / (1.0 - p))


def fit_platt_vectorized(p_cal: list[float], y_cal: list[int], iters: int = 4000, lr: float = 0.05):
    """Same GD update as ``calibrators.fit_platt``, with array reductions.

    The original candidate used a Python inner loop for every score and every
    gradient step. That made the identity-grouped B=200 protocol take tens of
    minutes while producing no additional audit value. NumPy is used only for
    deterministic scalar reductions; initialization, clipping, step size and
    iteration count remain the locked calibrator protocol.
    """
    x = np.asarray([logit_clip(p) for p in p_cal], dtype=np.float64)
    y = np.asarray(y_cal, dtype=np.float64)
    a, b = 1.0, 0.0
    n = len(x)
    if n == 0:
        raise RuntimeError("empty Platt fit")
    for _ in range(iters):
        z = a * x + b
        q = np.empty_like(z)
        positive = z >= 0.0
        q[positive] = 1.0 / (1.0 + np.exp(-z[positive]))
        ez = np.exp(z[~positive])
        q[~positive] = ez / (1.0 + ez)
        err = q - y
        a -= lr * float(np.dot(err, x)) / n
        b -= lr * float(err.sum()) / n

    def apply(probs: list[float]) -> list[float]:
        return [sigmoid(a * logit_clip(p) + b) for p in probs]

    apply.coeff = (a, b)  # type: ignore[attr-defined]
    return apply


def nested_crossfit(rows: list[dict], seed: int) -> dict[str, float]:
    ys = [r["y"] for r in rows]
    folds = xfit.assign_identity_folds(rows, xfit.K_FOLDS, seed)
    oof_u = [0.0] * len(rows)
    oof_p = [0.0] * len(rows)
    oof_i = [0.0] * len(rows)
    for k in range(xfit.K_FOLDS):
        fit_idx = [i for i, f in enumerate(folds) if f != k]
        test_idx = [i for i, f in enumerate(folds) if f == k]
        s_fit = [rows[i]["s"] for i in fit_idx]
        y_fit = [rows[i]["y"] for i in fit_idx]
        s_test = [rows[i]["s"] for i in test_idx]
        mu, sd = xfit.mean_sd(s_fit)
        p_fit = xfit.zscore_sigmoid(s_fit, mu, sd)
        p_test = xfit.zscore_sigmoid(s_test, mu, sd)
        if len(set(y_fit)) < 2:
            raise ValueError(f"nested fit fold {k} has one class; calibrator cannot be estimated")
        platt = fit_platt_vectorized(p_fit, y_fit)
        iso = fit_isotonic(p_fit, y_fit)
        q_p = platt(p_test)
        q_i = iso(p_test)
        for j, i in enumerate(test_idx):
            oof_u[i] = p_test[j]
            oof_p[i] = q_p[j]
            oof_i[i] = q_i[j]
    return {
        "uncalibrated": ece(oof_u, ys, N_BINS),
        "pooled_platt": ece(oof_p, ys, N_BINS),
        "pooled_isotonic": ece(oof_i, ys, N_BINS),
    }


def percentile(values: list[float], q: float) -> float:
    return xfit.percentile(values, q)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if not STOP.is_file():
        raise SystemExit(f"historical STOP missing: {STOP}")
    stop_hash = sha256(STOP)
    manifest = HERE / "inputs" / "manifest_gof_650M.csv"
    rows = xfit.load_manifest(manifest)
    if len(rows) != 4677:
        raise SystemExit(f"unexpected n={len(rows)}")

    parent = xfit.crossfit(rows, seed=xfit.SEED)
    uncal = parent["ece"]["uncalibrated"]["ece"]
    if abs(uncal - EXPECTED_UNCAL) > 1e-12:
        raise SystemExit(f"replay failed: uncal ECE {uncal} != {EXPECTED_UNCAL}")

    genes = parent["genes"]
    labels = parent["labels"]
    maps = {
        "uncalibrated": parent["oof"]["uncalibrated"],
        "pooled_platt": parent["oof"]["platt"],
        "pooled_isotonic": parent["oof"]["isotonic"],
    }
    folds = xfit.assign_folds(labels, xfit.K_FOLDS, xfit.SEED)
    pergene = r2.xfit_per_gene(maps["uncalibrated"], labels, genes, folds)
    maps["pergene_platt"] = pergene["platt"]

    gene_idx: dict[str, list[int]] = defaultdict(list)
    for i, gene in enumerate(genes):
        gene_idx[gene].append(i)

    bin_rows = []
    summary_rows = []
    for gene, idx in sorted(gene_idx.items()):
        y_g = [labels[i] for i in idx]
        n = len(idx)
        n_path = sum(y_g)
        for map_name in MAPS:
            p_g = [maps[map_name][i] for i in idx]
            gene_ece = ece(p_g, y_g, N_BINS)
            points = reliability_points(p_g, y_g, N_BINS)
            occupied = {int(pt["bin"]) for pt in points}
            for b in range(N_BINS):
                rec = {
                    "gene": gene,
                    "map": map_name,
                    "bin": b,
                    "bin_lo": b / N_BINS,
                    "bin_hi": (b + 1) / N_BINS,
                    "bin_mid": (b + 0.5) / N_BINS,
                    "n_gene": n,
                    "occupied": b in occupied,
                }
                match = next((pt for pt in points if int(pt["bin"]) == b), None)
                if match:
                    rec["mean_pred"] = r6(match["mean_pred"])
                    rec["obs_rate"] = r6(match["obs_rate"])
                    rec["n"] = int(match["n"])
                    rec["abs_gap"] = r6(abs(match["mean_pred"] - match["obs_rate"]))
                    rec["ece_share"] = r6((match["n"] / n) * abs(match["mean_pred"] - match["obs_rate"]))
                else:
                    rec["mean_pred"] = None
                    rec["obs_rate"] = None
                    rec["n"] = 0
                    rec["abs_gap"] = None
                    rec["ece_share"] = 0.0
                bin_rows.append(rec)
            summary_rows.append(
                {
                    "gene": gene,
                    "map": map_name,
                    "n": n,
                    "n_pathogenic": n_path,
                    "base_rate": r6(n_path / n if n else float("nan")),
                    "ece": r6(gene_ece),
                    "n_occupied_bins": len(occupied),
                    "low_n_flag": n < 40,
                }
            )

    write_csv(
        OUT / "per_gene_reliability_bins.csv",
        bin_rows,
        [
            "gene",
            "map",
            "bin",
            "bin_lo",
            "bin_hi",
            "bin_mid",
            "mean_pred",
            "obs_rate",
            "n",
            "n_gene",
            "abs_gap",
            "ece_share",
            "occupied",
        ],
    )
    write_csv(
        OUT / "per_gene_ece_by_map.csv",
        summary_rows,
        ["gene", "map", "n", "n_pathogenic", "base_rate", "ece", "n_occupied_bins", "low_n_flag"],
    )

    pooled_bins = []
    for map_name in MAPS:
        for pt in reliability_points(maps[map_name], labels, N_BINS):
            pooled_bins.append(
                {
                    "map": map_name,
                    **{k: r6(v) if k not in {"bin", "n"} else v for k, v in pt.items()},
                }
            )
    write_csv(
        OUT / "pooled_reliability_bins.csv",
        pooled_bins,
        ["map", "bin", "bin_lo", "bin_hi", "bin_mid", "mean_pred", "obs_rate", "n"],
    )

    print(f"nested bootstrap B={N_NESTED} seed={NESTED_SEED}", flush=True)
    rng = random.Random(NESTED_SEED)
    nested = {k: [] for k in ("uncalibrated", "pooled_platt", "pooled_isotonic")}
    n = len(rows)
    for b in range(N_NESTED):
        idx = [rng.randrange(n) for _ in range(n)]
        boot_rows = [rows[i] for i in idx]
        rec = nested_crossfit(boot_rows, seed=NESTED_SEED + b)
        for key, value in rec.items():
            nested[key].append(value)
        if (b + 1) % 20 == 0:
            print(f"  nested {b + 1}/{N_NESTED}", flush=True)

    nested_summary = {}
    for key, draws in nested.items():
        nested_summary[key] = {
            "n_boot": len(draws),
            "mean": r6(sum(draws) / len(draws)),
            "ci_lo": r6(percentile(draws, 0.025)),
            "ci_hi": r6(percentile(draws, 0.975)),
            "point_oof_pair_ece": r6(
                {
                    "uncalibrated": parent["ece"]["uncalibrated"]["ece"],
                    "pooled_platt": parent["ece"]["platt"]["ece"],
                    "pooled_isotonic": parent["ece"]["isotonic"]["ece"],
                }[key]
            ),
            "oof_pair_ci_lo": r6(
                {
                    "uncalibrated": parent["ece"]["uncalibrated"]["ci_lo"],
                    "pooled_platt": parent["ece"]["platt"]["ci_lo"],
                    "pooled_isotonic": parent["ece"]["isotonic"]["ci_lo"],
                }[key]
            ),
            "oof_pair_ci_hi": r6(
                {
                    "uncalibrated": parent["ece"]["uncalibrated"]["ci_hi"],
                    "pooled_platt": parent["ece"]["platt"]["ci_hi"],
                    "pooled_isotonic": parent["ece"]["isotonic"]["ci_hi"],
                }[key]
            ),
        }

    payload = {
        "protocol_id": PROTOCOL,
        "parent_protocol": PARENT,
        "not_a_restoration": True,
        "forbidden_claim": "Do not treat any number here as recovery of 0.0207 or 91%",
        "historical_stop_file": str(STOP),
        "historical_stop_sha256": stop_hash,
        "historical_stop_unchanged": True,
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "seed": xfit.SEED,
        "nested_seed": NESTED_SEED,
        "k_folds": xfit.K_FOLDS,
        "n_bins": N_BINS,
        "n_nested": N_NESTED,
        "n": len(rows),
        "n_genes": len(gene_idx),
        "input_sha256": sha256(manifest),
        "replayed_uncal_ece": uncal,
        "bin_definition": "equal-width 10 bins on [0,1]; bin b=[b/10,(b+1)/10) for b=0..8; bin 9=[0.9,1.0]",
        "nested_bootstrap": (
            "resample n rows with replacement; stratify original example identities into 5 folds; "
            "keep all copies of an original example in one fold; "
            "re-fit z-score sigmoid + Platt + isotonic on complementary folds; "
            "evaluate ECE on the bootstrap OOF map. Calibrators ARE refit."
        ),
        "oof_pair_bootstrap_note": "R1/R2 CIs resample OOF (p,y) pairs and do not refit calibrators",
        "supersedes_only": "The archived nested bootstrap; original unique-row OOF maps are unchanged.",
        "nested": nested_summary,
        "ready": "NO",
        "printed": False,
    }
    (OUT / "round3_summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    (OUT / "STATUS.txt").write_text(
        f"{PROTOCOL}\n"
        "Per-gene 10-bin reliability under the locked xfit maps, plus nested "
        "B=200 bootstrap that re-fits calibrators.\n"
        "NOT a restoration of 0.2443/0.0207. STOP file unchanged.\n"
        f"replayed_uncal_ece={uncal}\n"
        f"stop_sha256={stop_hash}\n"
        "Ready: not claimed. Official PDF not rebuilt.\n"
    )
    print("wrote", OUT, flush=True)
    print("nested", json.dumps(nested_summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
