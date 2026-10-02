#!/usr/bin/env python3
"""A3-XFIT-CAL-R2-20260918: second science-hardening round.

Explains the 7/20 pooled-Platt degradation, compares per-gene / grouped /
gene-intercept calibrators to the pooled map, reports Cox slope+intercept,
decision-curve net benefit, and Holm/BH-corrected per-gene inference.

Cross-fitted throughout. Same seed/folds as A3-XFIT-CAL-20260918.
Not a restoration of the retired historical pair (0.2443 -> 0.0207 / 91%).
"""
from __future__ import annotations

import csv
import json
import math
import platform
import random
import sys
from collections import defaultdict
from pathlib import Path

from calibrators import ece, fit_isotonic, fit_platt, logit, sigmoid
import run_xfit_calibration as xfit

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs" / "round2"
SEED = xfit.SEED
K_FOLDS = xfit.K_FOLDS
N_BOOT = xfit.N_BOOT
N_BINS = xfit.N_BINS
MIN_FIT = 8
THRESHOLDS = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
PROTOCOL = "A3-XFIT-CAL-R2-20260918"
EXPECTED_UNCAL = 0.24848355100110145


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


def round_row(row: dict) -> dict:
    return {k: r6(v) for k, v in row.items()}


def mean(xs: list[float]) -> float:
    return sum(xs) / len(xs) if xs else float("nan")


def pop_sd(xs: list[float]) -> float:
    if not xs:
        return float("nan")
    mu = mean(xs)
    return (sum((x - mu) ** 2 for x in xs) / len(xs)) ** 0.5


def rankdata(xs: list[float]) -> list[float]:
    n = len(xs)
    order = sorted(range(n), key=lambda i: xs[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = mean(xs), mean(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    dy = math.sqrt(sum((b - my) ** 2 for b in ys))
    if dx == 0.0 or dy == 0.0:
        return float("nan")
    return num / (dx * dy)


def spearman(xs: list[float], ys: list[float]) -> float:
    return pearson(rankdata(xs), rankdata(ys))


def ks_two_sample(xs: list[float], ys: list[float]) -> float:
    a = sorted(xs)
    b = sorted(ys)
    i = j = 0
    d = 0.0
    na, nb = len(a), len(b)
    while i < na and j < nb:
        if a[i] <= b[j]:
            i += 1
        else:
            j += 1
        d = max(d, abs(i / na - j / nb))
    return d


def solve_3x3(A: list[list[float]], rhs: list[float]) -> list[float]:
    m = [row[:] + [rhs[i]] for i, row in enumerate(A)]
    n = 3
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(m[r][col]))
        m[col], m[pivot] = m[pivot], m[col]
        if abs(m[col][col]) < 1e-12:
            return [float("nan")] * n
        scale = m[col][col]
        for k in range(col, n + 1):
            m[col][k] /= scale
        for r in range(n):
            if r == col:
                continue
            fac = m[r][col]
            for k in range(col, n + 1):
                m[r][k] -= fac * m[col][k]
    return [m[i][n] for i in range(n)]


def ols_two(x1: list[float], x2: list[float], y: list[float]) -> dict:
    n = len(y)
    xtx = [[0.0] * 3 for _ in range(3)]
    xty = [0.0, 0.0, 0.0]
    for a, b, yi in zip(x1, x2, y):
        row = (1.0, a, b)
        for i in range(3):
            xty[i] += row[i] * yi
            for j in range(3):
                xtx[i][j] += row[i] * row[j]
    beta = solve_3x3(xtx, xty)
    yhat = [beta[0] + beta[1] * a + beta[2] * b for a, b in zip(x1, x2)]
    resid = [yi - yh for yi, yh in zip(y, yhat)]
    ss_res = sum(r * r for r in resid)
    ym = mean(y)
    ss_tot = sum((yi - ym) ** 2 for yi in y)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {"b0": beta[0], "b1": beta[1], "b2": beta[2], "r2": r2, "n": n}


def holm_adjust(pvalues: list[float], alpha: float = 0.05) -> tuple[list[float], list[bool]]:
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    adj = [1.0] * m
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, pvalues[i] * (m - rank)))
        adj[i] = running
    return adj, [a <= alpha for a in adj]


def bh_qvalues(pvalues: list[float]) -> list[float]:
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    q = [1.0] * m
    running = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        running = min(running, pvalues[i] * m / (rank + 1))
        q[i] = min(1.0, running)
    return q


def two_sided_boot_p(draws: list[float]) -> float:
    b = len(draws)
    if b == 0:
        return float("nan")
    frac_pos = (1 + sum(1 for d in draws if d >= 0.0)) / (b + 1)
    frac_neg = (1 + sum(1 for d in draws if d <= 0.0)) / (b + 1)
    return min(1.0, 2.0 * min(frac_pos, frac_neg))


def fit_cox(probs: list[float], labels: list[int], max_iter: int = 30) -> tuple[float, float]:
    x = [logit(p) for p in probs]
    a, b = 1.0, 0.0
    n = len(x)
    if n < 4 or len(set(labels)) < 2:
        return float("nan"), float("nan")
    for _ in range(max_iter):
        ga = gb = 0.0
        haa = hbb = hab = 0.0
        for xi, yi in zip(x, labels):
            qi = sigmoid(a * xi + b)
            w = qi * (1.0 - qi)
            err = qi - yi
            ga += err * xi
            gb += err
            haa += w * xi * xi
            hbb += w
            hab += w * xi
        det = haa * hbb - hab * hab
        if abs(det) < 1e-12:
            break
        da = (hbb * ga - hab * gb) / det
        db = (haa * gb - hab * ga) / det
        a -= da
        b -= db
        if abs(da) + abs(db) < 1e-12:
            break
    return a, b


def fit_platt_gene_offset(
    p_cal: list[float],
    y_cal: list[int],
    g_cal: list[str],
    genes: list[str],
    iters: int = 4000,
    lr: float = 0.05,
):
    x = [logit(p) for p in p_cal]
    index = {g: i for i, g in enumerate(genes)}
    a, b = 1.0, 0.0
    c = [0.0] * len(genes)
    n = len(x)
    counts = [0] * len(genes)
    for g in g_cal:
        counts[index[g]] += 1
    for _ in range(iters):
        ga = gb = 0.0
        gc = [0.0] * len(genes)
        for xi, yi, g in zip(x, y_cal, g_cal):
            k = index[g]
            qi = sigmoid(a * xi + b + c[k])
            err = qi - yi
            ga += err * xi
            gb += err
            gc[k] += err
        a -= lr * ga / n
        b -= lr * gb / n
        for k in range(len(genes)):
            if counts[k] == 0:
                continue
            c[k] -= lr * gc[k] / n
        mu_c = sum(c) / len(c)
        b += mu_c
        c = [ck - mu_c for ck in c]

    def apply(probs: list[float], gene_ids: list[str]) -> list[float]:
        out = []
        for p, g in zip(probs, gene_ids):
            k = index.get(g)
            off = 0.0 if k is None else c[k]
            out.append(sigmoid(a * logit(p) + b + off))
        return out

    apply.coeff = (a, b, {g: c[index[g]] for g in genes})  # type: ignore[attr-defined]
    return apply


def net_benefit(probs: list[float], labels: list[int], t: float) -> float:
    n = len(labels)
    if n == 0 or t <= 0.0 or t >= 1.0:
        return float("nan")
    tp = sum(1 for p, y in zip(probs, labels) if p >= t and y == 1)
    fp = sum(1 for p, y in zip(probs, labels) if p >= t and y == 0)
    return tp / n - (fp / n) * (t / (1.0 - t))


def treat_all_nb(labels: list[int], t: float) -> float:
    prev = mean(labels)
    return prev - (1.0 - prev) * (t / (1.0 - t))


def boot_pairs(n: int, seed: int, n_boot: int):
    rng = random.Random(seed)
    for _ in range(n_boot):
        yield [rng.randrange(n) for _ in range(n)]


def can_fit(y: list[int]) -> bool:
    return len(y) >= MIN_FIT and len(set(y)) >= 2


def group_indices(genes: list[str]) -> dict[str, list[int]]:
    out: dict[str, list[int]] = defaultdict(list)
    for i, g in enumerate(genes):
        out[g].append(i)
    return dict(out)


def apply_identity(probs: list[float]) -> list[float]:
    return list(probs)


def xfit_restricted(
    uncal: list[float],
    labels: list[int],
    genes: list[str],
    folds: list[int],
    allowed,
    scheme: str,
) -> dict:
    """Fit calibrators on complementary folds restricted to `allowed` rows.

    `allowed(i)` is a predicate on the *fit* index. Test rows that fail the
    matching test-side predicate are left as None and filled with identity.
    """
    n = len(uncal)
    oof_platt = [None] * n
    oof_iso = [None] * n
    fold_notes = []
    for k in range(K_FOLDS):
        fit_idx = [i for i in range(n) if folds[i] != k and allowed(i)]
        test_idx = [i for i in range(n) if folds[i] == k and allowed(i)]
        y_fit = [labels[i] for i in fit_idx]
        note = {
            "scheme": scheme,
            "fold": k,
            "n_fit": len(fit_idx),
            "n_eval": len(test_idx),
            "n_fit_pathogenic": sum(y_fit),
            "skipped": False,
        }
        if not can_fit(y_fit) or not test_idx:
            for i in test_idx:
                oof_platt[i] = uncal[i]
                oof_iso[i] = uncal[i]
            note["skipped"] = True
            fold_notes.append(note)
            continue
        p_fit = [uncal[i] for i in fit_idx]
        p_test = [uncal[i] for i in test_idx]
        platt = fit_platt(p_fit, y_fit)
        iso = fit_isotonic(p_fit, y_fit)
        q_p = platt(p_test)
        q_i = iso(p_test)
        for j, i in enumerate(test_idx):
            oof_platt[i] = q_p[j]
            oof_iso[i] = q_i[j]
        fold_notes.append(note)
    for i in range(n):
        if oof_platt[i] is None:
            oof_platt[i] = uncal[i]
            oof_iso[i] = uncal[i]
    return {
        "platt": oof_platt,
        "isotonic": oof_iso,
        "fold_notes": fold_notes,
        "n_skipped_folds": sum(1 for r in fold_notes if r["skipped"]),
    }


def xfit_per_gene(uncal, labels, genes, folds) -> dict:
    n = len(uncal)
    oof_p = [0.0] * n
    oof_i = [0.0] * n
    notes = []
    skipped = {}
    for gene, idx in group_indices(genes).items():
        allowed = set(idx)
        rec = xfit_restricted(
            uncal, labels, genes, folds, lambda i, a=allowed: i in a, f"pergene:{gene}"
        )
        for i in idx:
            oof_p[i] = rec["platt"][i]
            oof_i[i] = rec["isotonic"][i]
        skipped[gene] = rec["n_skipped_folds"]
        notes.extend(rec["fold_notes"])
    return {
        "platt": oof_p,
        "isotonic": oof_i,
        "skipped_folds": skipped,
        "fold_notes": notes,
    }


def xfit_grouped(uncal, labels, genes, folds) -> dict:
    n = len(uncal)
    oof_p = [0.0] * n
    oof_i = [0.0] * n
    notes = []
    for family, members in (
        ("ORIG6", set(xfit.ORIG6)),
        ("OTHER14", None),
    ):
        if family == "OTHER14":
            pred = lambda i, m=set(xfit.ORIG6): genes[i] not in m
        else:
            pred = lambda i, m=members: genes[i] in m
        rec = xfit_restricted(uncal, labels, genes, folds, pred, f"grouped:{family}")
        for i in range(n):
            if pred(i):
                oof_p[i] = rec["platt"][i]
                oof_i[i] = rec["isotonic"][i]
        notes.extend(rec["fold_notes"])
    return {"platt": oof_p, "isotonic": oof_i, "fold_notes": notes}


def xfit_gene_offset(uncal, labels, genes, folds) -> dict:
    n = len(uncal)
    gene_list = sorted(set(genes))
    oof = [0.0] * n
    coeff = []
    for k in range(K_FOLDS):
        fit_idx = [i for i in range(n) if folds[i] != k]
        test_idx = [i for i in range(n) if folds[i] == k]
        model = fit_platt_gene_offset(
            [uncal[i] for i in fit_idx],
            [labels[i] for i in fit_idx],
            [genes[i] for i in fit_idx],
            gene_list,
        )
        q = model([uncal[i] for i in test_idx], [genes[i] for i in test_idx])
        for j, i in enumerate(test_idx):
            oof[i] = q[j]
        a, b, c = model.coeff  # type: ignore[attr-defined]
        coeff.append({"fold": k, "a": a, "b": b, "offsets": c})
    return {"platt": oof, "coeff": coeff}


def logo_calibrate(uncal, labels, genes) -> dict:
    """Leave-one-gene-out: fit on all other genes, apply to the held-out gene.

    No target-gene labels enter the calibrator. Not fold-cross-fit; the held-out
    unit is the gene.
    """
    n = len(uncal)
    oof_p = [0.0] * n
    oof_i = [0.0] * n
    notes = []
    by_gene = group_indices(genes)
    for gene, idx in by_gene.items():
        fit_idx = [i for i in range(n) if genes[i] != gene]
        y_fit = [labels[i] for i in fit_idx]
        if not can_fit(y_fit):
            for i in idx:
                oof_p[i] = uncal[i]
                oof_i[i] = uncal[i]
            notes.append({"gene": gene, "skipped": True, "n_fit": len(fit_idx)})
            continue
        p_fit = [uncal[i] for i in fit_idx]
        p_test = [uncal[i] for i in idx]
        platt = fit_platt(p_fit, y_fit)
        iso = fit_isotonic(p_fit, y_fit)
        q_p = platt(p_test)
        q_i = iso(p_test)
        for j, i in enumerate(idx):
            oof_p[i] = q_p[j]
            oof_i[i] = q_i[j]
        notes.append({"gene": gene, "skipped": False, "n_fit": len(fit_idx)})
    return {"platt": oof_p, "isotonic": oof_i, "notes": notes}


def boot_ece_ci(probs, labels, seed, n_boot=N_BOOT) -> dict:
    return xfit.boot_ece(probs, labels, seed=seed, n_boot=n_boot)


def boot_delta(p_a, p_b, labels, seed, n_boot=N_BOOT) -> dict:
    n = len(labels)
    point = ece(p_b, labels, N_BINS) - ece(p_a, labels, N_BINS)
    draws = []
    for idx in boot_pairs(n, seed, n_boot):
        ya = [labels[i] for i in idx]
        draws.append(
            ece([p_b[i] for i in idx], ya, N_BINS) - ece([p_a[i] for i in idx], ya, N_BINS)
        )
    return {
        "delta": point,
        "ci_lo": xfit.percentile(draws, 0.025),
        "ci_hi": xfit.percentile(draws, 0.975),
        "p": two_sided_boot_p(draws),
        "draws": draws,
    }


def boot_cox(probs, labels, seed, n_boot=N_BOOT) -> dict:
    slope, intercept = fit_cox(probs, labels)
    n = len(labels)
    s_draws, i_draws = [], []
    for idx in boot_pairs(n, seed, n_boot):
        a, b = fit_cox([probs[i] for i in idx], [labels[i] for i in idx])
        if math.isnan(a):
            continue
        s_draws.append(a)
        i_draws.append(b)
    return {
        "slope": slope,
        "slope_ci_lo": xfit.percentile(s_draws, 0.025) if s_draws else float("nan"),
        "slope_ci_hi": xfit.percentile(s_draws, 0.975) if s_draws else float("nan"),
        "intercept": intercept,
        "intercept_ci_lo": xfit.percentile(i_draws, 0.025) if i_draws else float("nan"),
        "intercept_ci_hi": xfit.percentile(i_draws, 0.975) if i_draws else float("nan"),
        "n_boot_ok": len(s_draws),
    }


def boot_nb(probs, labels, t, seed, n_boot=N_BOOT) -> dict:
    point = net_benefit(probs, labels, t)
    n = len(labels)
    draws = [net_benefit([probs[i] for i in idx], [labels[i] for i in idx], t) for idx in boot_pairs(n, seed, n_boot)]
    return {
        "nb": point,
        "ci_lo": xfit.percentile(draws, 0.025),
        "ci_hi": xfit.percentile(draws, 0.975),
    }


def reconstruct_folds(rows: list[dict]) -> list[int]:
    return xfit.assign_folds([r["y"] for r in rows], K_FOLDS, SEED)


def perm_spearman(xs, ys, seed, n_perm=2000) -> dict:
    obs = spearman(xs, ys)
    rng = random.Random(seed)
    y = list(ys)
    count = 0
    for _ in range(n_perm):
        rng.shuffle(y)
        val = spearman(xs, y)
        if abs(val) >= abs(obs):
            count += 1
    return {"rho": obs, "p": (1 + count) / (n_perm + 1), "n_perm": n_perm}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gof_path = xfit.INPUTS / "manifest_gof_650M.csv"
    rows = xfit.load_manifest(gof_path)
    if len(rows) != 4677:
        raise SystemExit(f"expected n=4677, got {len(rows)}")
    pooled = xfit.crossfit(rows, seed=SEED)
    uncal = pooled["oof"]["uncalibrated"]
    platt = pooled["oof"]["platt"]
    iso = pooled["oof"]["isotonic"]
    labels = pooled["labels"]
    genes = pooled["genes"]
    scores = [r["s"] for r in rows]
    folds = reconstruct_folds(rows)
    n = len(rows)
    if abs(pooled["ece"]["uncalibrated"]["ece"] - EXPECTED_UNCAL) > 1e-12:
        raise SystemExit("pooled uncalibrated ECE does not replay A3-XFIT-CAL-20260918")

    print("=== fitting restricted calibrators ===", flush=True)
    pergene = xfit_per_gene(uncal, labels, genes, folds)
    grouped = xfit_grouped(uncal, labels, genes, folds)
    geneoff = xfit_gene_offset(uncal, labels, genes, folds)
    logo = logo_calibrate(uncal, labels, genes)

    schemes = {
        "uncalibrated": uncal,
        "pooled_platt": platt,
        "pooled_isotonic": iso,
        "pergene_platt": pergene["platt"],
        "pergene_isotonic": pergene["isotonic"],
        "grouped_platt": grouped["platt"],
        "grouped_isotonic": grouped["isotonic"],
        "geneoffset_platt": geneoff["platt"],
        "logo_platt": logo["platt"],
        "logo_isotonic": logo["isotonic"],
    }

    print("=== pooled scheme ECE ===", flush=True)
    scheme_rows = []
    scheme_summary = {}
    for name, probs in schemes.items():
        rec = boot_ece_ci(probs, labels, seed=SEED)
        if name != "uncalibrated":
            dlt = boot_delta(uncal, probs, labels, seed=SEED + 1)
        else:
            dlt = {"delta": 0.0, "ci_lo": 0.0, "ci_hi": 0.0, "p": 1.0}
        vs_pooled = None
        if name not in ("uncalibrated", "pooled_platt", "pooled_isotonic"):
            if name.endswith("_isotonic"):
                vs_pooled = boot_delta(iso, probs, labels, seed=SEED + 2)
            else:
                vs_pooled = boot_delta(platt, probs, labels, seed=SEED + 2)
        row = {
            "scheme": name,
            "ece": rec["ece"],
            "ece_ci_lo": rec["ci_lo"],
            "ece_ci_hi": rec["ci_hi"],
            "delta_vs_uncal": dlt["delta"],
            "delta_vs_uncal_ci_lo": dlt["ci_lo"],
            "delta_vs_uncal_ci_hi": dlt["ci_hi"],
            "delta_vs_uncal_p": dlt["p"],
            "delta_vs_pooled": None if vs_pooled is None else vs_pooled["delta"],
            "delta_vs_pooled_ci_lo": None if vs_pooled is None else vs_pooled["ci_lo"],
            "delta_vs_pooled_ci_hi": None if vs_pooled is None else vs_pooled["ci_hi"],
            "delta_vs_pooled_p": None if vs_pooled is None else vs_pooled["p"],
            "n": n,
            "n_boot": N_BOOT,
            "seed": SEED,
        }
        scheme_rows.append(round_row(row))
        scheme_summary[name] = row
        print(
            f"  {name:20s} ECE={rec['ece']:.4f} [{rec['ci_lo']:.4f},{rec['ci_hi']:.4f}] "
            f"dUncal={dlt['delta']:+.4f}",
            flush=True,
        )

    print("=== per-gene mechanism / slope / multiplicity ===", flush=True)
    by_gene = group_indices(genes)
    pooled_rate = mean(labels)
    pooled_mean_s = mean(scores)
    pooled_sd_s = pop_sd(scores)
    gene_rows = []
    pvals = []
    for gene in sorted(by_gene):
        idx = by_gene[gene]
        y = [labels[i] for i in idx]
        s = [scores[i] for i in idx]
        p_u = [uncal[i] for i in idx]
        p_p = [platt[i] for i in idx]
        p_i = [iso[i] for i in idx]
        p_pg = [pergene["platt"][i] for i in idx]
        p_pgi = [pergene["isotonic"][i] for i in idx]
        p_go = [geneoff["platt"][i] for i in idx]
        p_lg = [logo["platt"][i] for i in idx]
        n_g = len(idx)
        n_pos = sum(y)
        base = n_pos / n_g
        e_u = ece(p_u, y, N_BINS)
        e_p = ece(p_p, y, N_BINS)
        e_i = ece(p_i, y, N_BINS)
        e_pg = ece(p_pg, y, N_BINS)
        e_pgi = ece(p_pgi, y, N_BINS)
        e_go = ece(p_go, y, N_BINS)
        e_lg = ece(p_lg, y, N_BINS)
        d_pooled = boot_delta(p_u, p_p, y, seed=SEED)
        d_iso = boot_delta(p_u, p_i, y, seed=SEED + 3)
        d_pg = boot_delta(p_u, p_pg, y, seed=SEED + 4)
        cox_u = boot_cox(p_u, y, seed=SEED)
        cox_p = boot_cox(p_p, y, seed=SEED + 5)
        gene_rows.append(
            {
                "gene": gene,
                "n": n_g,
                "n_pathogenic": n_pos,
                "n_benign": n_g - n_pos,
                "base_rate": base,
                "base_minus_pooled": base - pooled_rate,
                "mean_s": mean(s),
                "sd_s": pop_sd(s),
                "s_shift_vs_pooled": mean(s) - pooled_mean_s,
                "sd_s_ratio": pop_sd(s) / pooled_sd_s if pooled_sd_s else float("nan"),
                "mean_p_uncal": mean(p_u),
                "mean_p_pooled_platt": mean(p_p),
                "mean_p_pergene_platt": mean(p_pg),
                "citl_uncal": mean(p_u) - base,
                "citl_pooled_platt": mean(p_p) - base,
                "citl_pergene_platt": mean(p_pg) - base,
                "ks_uncal_vs_pooled": ks_two_sample(p_u, uncal),
                "ece_uncal": e_u,
                "ece_uncal_ci_lo": d_pooled["ci_lo"] and None,
                "ece_pooled_platt": e_p,
                "ece_pooled_isotonic": e_i,
                "ece_pergene_platt": e_pg,
                "ece_pergene_isotonic": e_pgi,
                "ece_geneoffset_platt": e_go,
                "ece_logo_platt": e_lg,
                "delta_pooled_platt": d_pooled["delta"],
                "delta_pooled_platt_ci_lo": d_pooled["ci_lo"],
                "delta_pooled_platt_ci_hi": d_pooled["ci_hi"],
                "delta_pooled_platt_p": d_pooled["p"],
                "delta_pooled_isotonic": d_iso["delta"],
                "delta_pooled_isotonic_ci_lo": d_iso["ci_lo"],
                "delta_pooled_isotonic_ci_hi": d_iso["ci_hi"],
                "delta_pergene_platt": d_pg["delta"],
                "delta_pergene_platt_ci_lo": d_pg["ci_lo"],
                "delta_pergene_platt_ci_hi": d_pg["ci_hi"],
                "degrades_pooled_platt": e_p > e_u,
                "degrades_pergene_platt": e_pg > e_u,
                "pergene_skipped_folds": pergene["skipped_folds"].get(gene, 0),
                "slope_uncal": cox_u["slope"],
                "slope_uncal_ci_lo": cox_u["slope_ci_lo"],
                "slope_uncal_ci_hi": cox_u["slope_ci_hi"],
                "intercept_uncal": cox_u["intercept"],
                "intercept_uncal_ci_lo": cox_u["intercept_ci_lo"],
                "intercept_uncal_ci_hi": cox_u["intercept_ci_hi"],
                "slope_pooled_platt": cox_p["slope"],
                "slope_pooled_platt_ci_lo": cox_p["slope_ci_lo"],
                "slope_pooled_platt_ci_hi": cox_p["slope_ci_hi"],
                "intercept_pooled_platt": cox_p["intercept"],
                "intercept_pooled_platt_ci_lo": cox_p["intercept_ci_lo"],
                "intercept_pooled_platt_ci_hi": cox_p["intercept_ci_hi"],
                "low_n_flag": n_g < 30,
            }
        )
        # fill uncal ECE CI separately (cheap)
        u_ci = boot_ece_ci(p_u, y, seed=SEED)
        p_ci = boot_ece_ci(p_p, y, seed=SEED)
        i_ci = boot_ece_ci(p_i, y, seed=SEED)
        gene_rows[-1]["ece_uncal_ci_lo"] = u_ci["ci_lo"]
        gene_rows[-1]["ece_uncal_ci_hi"] = u_ci["ci_hi"]
        gene_rows[-1]["ece_pooled_platt_ci_lo"] = p_ci["ci_lo"]
        gene_rows[-1]["ece_pooled_platt_ci_hi"] = p_ci["ci_hi"]
        gene_rows[-1]["ece_pooled_isotonic_ci_lo"] = i_ci["ci_lo"]
        gene_rows[-1]["ece_pooled_isotonic_ci_hi"] = i_ci["ci_hi"]
        pvals.append(d_pooled["p"])
        print(
            f"  {gene:16s} n={n_g:4d} rate={base:.3f} "
            f"u={e_u:.3f} P={e_p:.3f} d={d_pooled['delta']:+.3f} p={d_pooled['p']:.4f}",
            flush=True,
        )

    holm_p, holm_rej = holm_adjust(pvals, alpha=0.05)
    qvals = bh_qvalues(pvals)
    for rec, hp, hr, q in zip(gene_rows, holm_p, holm_rej, qvals):
        rec["p_holm"] = hp
        rec["reject_holm_05"] = hr
        rec["q_bh"] = q
        rec["reject_bh_05"] = q <= 0.05

    print("=== mechanism correlations ===", flush=True)
    deltas = [r["delta_pooled_platt"] for r in gene_rows]
    corr_specs = [
        ("n", [r["n"] for r in gene_rows]),
        ("log_n", [math.log(r["n"]) for r in gene_rows]),
        ("base_rate", [r["base_rate"] for r in gene_rows]),
        ("abs_base_minus_pooled", [abs(r["base_minus_pooled"]) for r in gene_rows]),
        ("base_minus_pooled", [r["base_minus_pooled"] for r in gene_rows]),
        ("ece_uncal", [r["ece_uncal"] for r in gene_rows]),
        ("citl_uncal", [r["citl_uncal"] for r in gene_rows]),
        ("citl_pooled_platt", [r["citl_pooled_platt"] for r in gene_rows]),
        ("s_shift_vs_pooled", [r["s_shift_vs_pooled"] for r in gene_rows]),
        ("sd_s", [r["sd_s"] for r in gene_rows]),
        ("ks_uncal_vs_pooled", [r["ks_uncal_vs_pooled"] for r in gene_rows]),
        ("mean_p_uncal", [r["mean_p_uncal"] for r in gene_rows]),
    ]
    corr_rows = []
    for name, xs in corr_specs:
        rec = perm_spearman(xs, deltas, seed=SEED)
        rec["predictor"] = name
        rec["outcome"] = "delta_pooled_platt"
        corr_rows.append(rec)
        print(f"  Spearman({name}, dECE)={rec['rho']:+.3f} p={rec['p']:.4f}", flush=True)
    ols = ols_two(
        [r["base_rate"] for r in gene_rows],
        [math.log(r["n"]) for r in gene_rows],
        deltas,
    )
    print(
        f"  OLS dECE ~ base_rate + log(n): b_rate={ols['b1']:+.4f} "
        f"b_logn={ols['b2']:+.4f} R2={ols['r2']:.3f}",
        flush=True,
    )

    print("=== Cox slope/intercept (pooled OOF) ===", flush=True)
    slope_rows = []
    for name, probs in (
        ("uncalibrated", uncal),
        ("pooled_platt", platt),
        ("pooled_isotonic", iso),
        ("pergene_platt", pergene["platt"]),
        ("geneoffset_platt", geneoff["platt"]),
    ):
        rec = boot_cox(probs, labels, seed=SEED)
        rec["scheme"] = name
        rec["n"] = n
        slope_rows.append(rec)
        print(
            f"  {name:20s} slope={rec['slope']:.3f} [{rec['slope_ci_lo']:.3f},{rec['slope_ci_hi']:.3f}] "
            f"int={rec['intercept']:.3f} [{rec['intercept_ci_lo']:.3f},{rec['intercept_ci_hi']:.3f}]",
            flush=True,
        )

    print("=== decision-curve net benefit ===", flush=True)
    dca_rows = []
    for t in THRESHOLDS:
        ta = treat_all_nb(labels, t)
        block = {"threshold": t, "treat_all": ta, "treat_none": 0.0, "prevalence": pooled_rate}
        for name, probs in (
            ("uncalibrated", uncal),
            ("pooled_platt", platt),
            ("pooled_isotonic", iso),
            ("pergene_platt", pergene["platt"]),
        ):
            rec = boot_nb(probs, labels, t, seed=SEED)
            calls = sum(1 for p in probs if p >= t)
            flips = sum(
                1
                for a, b in zip(uncal, probs)
                if (a >= t) != (b >= t)
            )
            block[f"nb_{name}"] = rec["nb"]
            block[f"nb_{name}_ci_lo"] = rec["ci_lo"]
            block[f"nb_{name}_ci_hi"] = rec["ci_hi"]
            block[f"n_called_{name}"] = calls
            block[f"frac_flip_vs_uncal_{name}"] = flips / n
        d_platt = boot_delta_nb(uncal, platt, labels, t, seed=SEED + 6)
        block["nb_platt_minus_uncal"] = d_platt["delta"]
        block["nb_platt_minus_uncal_ci_lo"] = d_platt["ci_lo"]
        block["nb_platt_minus_uncal_ci_hi"] = d_platt["ci_hi"]
        dca_rows.append(block)
        print(
            f"  t={t:.2f} uncal={block['nb_uncalibrated']:.4f} "
            f"platt={block['nb_pooled_platt']:.4f} iso={block['nb_pooled_isotonic']:.4f} "
            f"treat_all={ta:.4f} dNB={d_platt['delta']:+.4f}",
            flush=True,
        )

    degraders = [r["gene"] for r in gene_rows if r["degrades_pooled_platt"]]
    holm_deg = [
        r["gene"]
        for r in gene_rows
        if r["degrades_pooled_platt"] and r["reject_holm_05"]
    ]
    holm_imp = [
        r["gene"]
        for r in gene_rows
        if (not r["degrades_pooled_platt"]) and r["reject_holm_05"]
    ]
    pergene_beats = scheme_summary["pergene_platt"]["ece"] < scheme_summary["pooled_platt"]["ece"]
    grouped_beats = scheme_summary["grouped_platt"]["ece"] < scheme_summary["pooled_platt"]["ece"]
    offset_beats = scheme_summary["geneoffset_platt"]["ece"] < scheme_summary["pooled_platt"]["ece"]

    summary = {
        "protocol_id": PROTOCOL,
        "parent_protocol": "A3-XFIT-CAL-20260918",
        "not_a_restoration": True,
        "forbidden_claim": "Do not treat any number here as recovery of 0.0207 or 91%",
        "python": sys.version,
        "platform": platform.platform(),
        "seed": SEED,
        "k_folds": K_FOLDS,
        "n_boot": N_BOOT,
        "n_bins": N_BINS,
        "min_fit": MIN_FIT,
        "thresholds": list(THRESHOLDS),
        "n": n,
        "n_pathogenic": sum(labels),
        "pooled_base_rate": pooled_rate,
        "input_sha256": xfit.sha256(gof_path),
        "replayed_uncal_ece": pooled["ece"]["uncalibrated"]["ece"],
        "schemes": {
            k: {
                "ece": v["ece"],
                "ci_lo": v["ece_ci_lo"],
                "ci_hi": v["ece_ci_hi"],
                "delta_vs_uncal": v["delta_vs_uncal"],
                "delta_vs_pooled": v["delta_vs_pooled"],
            }
            for k, v in scheme_summary.items()
        },
        "n_degrade_pooled_platt": len(degraders),
        "degraders_pooled_platt": degraders,
        "holm_significant_degraders": holm_deg,
        "holm_significant_improvers": holm_imp,
        "pergene_platt_beats_pooled": pergene_beats,
        "grouped_platt_beats_pooled": grouped_beats,
        "geneoffset_platt_beats_pooled": offset_beats,
        "ols_delta_on_base_and_logn": ols,
        "correlations": corr_rows,
        "pergene_skipped_folds": pergene["skipped_folds"],
        "geneoffset_fold_coeff": [
            {"fold": c["fold"], "a": c["a"], "b": c["b"], "offsets": c["offsets"]}
            for c in geneoff["coeff"]
        ],
    }

    write_csv(
        OUT / "scheme_ece.csv",
        scheme_rows,
        [
            "scheme",
            "ece",
            "ece_ci_lo",
            "ece_ci_hi",
            "delta_vs_uncal",
            "delta_vs_uncal_ci_lo",
            "delta_vs_uncal_ci_hi",
            "delta_vs_uncal_p",
            "delta_vs_pooled",
            "delta_vs_pooled_ci_lo",
            "delta_vs_pooled_ci_hi",
            "delta_vs_pooled_p",
            "n",
            "n_boot",
            "seed",
        ],
    )
    gene_fields = list(gene_rows[0].keys())
    write_csv(OUT / "per_gene_mechanism.csv", [round_row(r) for r in gene_rows], gene_fields)
    write_csv(
        OUT / "mechanism_correlations.csv",
        [round_row(r) for r in corr_rows],
        ["predictor", "outcome", "rho", "p", "n_perm"],
    )
    write_csv(
        OUT / "slope_intercept.csv",
        [round_row(r) for r in slope_rows],
        [
            "scheme",
            "slope",
            "slope_ci_lo",
            "slope_ci_hi",
            "intercept",
            "intercept_ci_lo",
            "intercept_ci_hi",
            "n",
            "n_boot_ok",
        ],
    )
    dca_fields = list(dca_rows[0].keys())
    write_csv(OUT / "decision_curve.csv", [round_row(r) for r in dca_rows], dca_fields)
    (OUT / "round2_summary.json").write_text(json.dumps(json_safe(summary), indent=2) + "\n")
    print("=== wrote", OUT, "===", flush=True)


def boot_delta_nb(p_a, p_b, labels, t, seed, n_boot=N_BOOT) -> dict:
    n = len(labels)
    point = net_benefit(p_b, labels, t) - net_benefit(p_a, labels, t)
    draws = []
    for idx in boot_pairs(n, seed, n_boot):
        ya = [labels[i] for i in idx]
        draws.append(
            net_benefit([p_b[i] for i in idx], ya, t) - net_benefit([p_a[i] for i in idx], ya, t)
        )
    return {
        "delta": point,
        "ci_lo": xfit.percentile(draws, 0.025),
        "ci_hi": xfit.percentile(draws, 0.975),
    }


def json_safe(obj):
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_safe(v) for v in obj]
    if isinstance(obj, float):
        return None if math.isnan(obj) else obj
    if isinstance(obj, bool):
        return obj
    return obj


if __name__ == "__main__":
    main()
