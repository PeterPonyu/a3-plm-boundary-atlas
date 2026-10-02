#!/usr/bin/env python3
"""One row-level scorer for A3 discrimination, calibration and uncertainty.

No historic scores or AUROC-derived scores are accepted. A successful N6c
recovery provenance is required. Outputs are new observations, not an attempt
to reproduce the invalid old calibration numbers.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
CAP = HERE.parents[1]
sys.path.insert(0, str(CAP / "analysis/xfit_calibration_20260918"))
from run_xfit_calibration import assign_folds

SEED = 20260923
K = 5
MAPS = ("sigmoid_z", "platt", "isotonic", "pooled_prevalence", "gene_prevalence")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(path, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def isotonic_fit(x, y, weight):
    """Weighted PAV on distinct scores; step map, constant at outer bounds."""
    unique, inverse = np.unique(x, return_inverse=True)
    mass = np.bincount(inverse, weights=weight)
    target = np.bincount(inverse, weights=weight*y)
    blocks = []
    for j, (w, s) in enumerate(zip(mass, target)):
        blocks.append([j, j+1, w, s])
        while len(blocks) > 1 and blocks[-2][3]/blocks[-2][2] >= blocks[-1][3]/blocks[-1][2]:
            b = blocks.pop()
            a = blocks.pop()
            blocks.append([a[0], b[1], a[2]+b[2], a[3]+b[3]])
    values = np.empty(len(unique))
    for start, stop, w, s in blocks:
        values[start:stop] = s/w

    def predict(q):
        index = np.clip(np.searchsorted(unique, q, side="right")-1, 0, len(unique)-1)
        return np.clip(values[index], 1e-6, 1-1e-6)
    return predict


def platt_fit(z, y, weight):
    design = np.column_stack((np.ones(len(z)), z))
    def objective(beta):
        eta = design @ beta
        value = np.sum(weight*(np.logaddexp(0, eta)-y*eta))/weight.sum()
        gradient = design.T @ (weight*(expit(eta)-y))/weight.sum()
        return value, gradient
    fit = minimize(objective, [0., 1.], jac=True, method="BFGS",
                   options={"gtol": 1e-8, "maxiter": 200})
    if np.max(np.abs(objective(fit.x)[1])) > 1e-6:
        raise ValueError(f"Logistic fit did not converge: {fit.message}")
    return lambda q: expit(fit.x[0]+fit.x[1]*q)


def crossfit(rows, weight=None, seed=SEED):
    """Each array entry is one original identity; weights retain boot copies."""
    n = len(rows)
    if len({r["example_id"] for r in rows}) != n:
        raise ValueError("One row per original identity is required")
    score = -np.array([float(r["primary_llr_fp32"]) for r in rows])
    y = 1-np.array([int(r["label_benign"]) for r in rows])
    genes = np.array([r["gene"] for r in rows])
    weight = np.ones(n) if weight is None else np.asarray(weight, dtype=float)
    active = np.flatnonzero(weight > 0)
    folds = np.full(n, -1, dtype=int)
    folds[active] = assign_folds(y[active].tolist(), K, seed)
    if set(folds[active]) != set(range(K)):
        raise ValueError("Every evaluation fold must be nonempty")
    predictions = {name: np.full(n, np.nan) for name in MAPS}
    for fold in range(K):
        train = (folds >= 0) & (folds != fold)
        test = folds == fold
        fit_ids = {rows[i]["example_id"] for i in np.flatnonzero(train)}
        test_ids = {rows[i]["example_id"] for i in np.flatnonzero(test)}
        if fit_ids & test_ids:
            raise AssertionError("Original identities cross fitting/evaluation")
        wt = weight[train]
        yt = y[train]
        if len(set(yt)) != 2:
            raise ValueError("Single-class calibrator fit is not estimable")
        mu = np.average(score[train], weights=wt)
        sd = np.sqrt(np.average((score[train]-mu)**2, weights=wt))
        if not sd > 0:
            raise ValueError("Constant fitting scores")
        zfit = (score[train]-mu)/sd
        ztest = (score[test]-mu)/sd
        predictions["sigmoid_z"][test] = expit(ztest)
        predictions["platt"][test] = platt_fit(zfit, yt, wt)(ztest)
        predictions["isotonic"][test] = isotonic_fit(score[train], yt, wt)(score[test])
        prevalence = np.average(yt, weights=wt)
        predictions["pooled_prevalence"][test] = prevalence
        for gene in np.unique(genes[test]):
            member = train & (genes == gene)
            if not member.any():
                raise ValueError(f"No training examples for gene baseline: {gene}")
            predictions["gene_prevalence"][test & (genes == gene)] = np.average(y[member], weights=weight[member])
    return predictions, y, folds


def metrics(p, y, weight=None):
    weight = np.ones(len(p)) if weight is None else np.asarray(weight)
    active = weight > 0
    p, y, weight = np.asarray(p)[active], np.asarray(y)[active], weight[active]
    if not np.all(np.isfinite(p)):
        raise ValueError("Missing out-of-fold predictions")
    bins = np.minimum(9, (p*10).astype(int))
    mass = np.bincount(bins, weights=weight, minlength=10)
    pred = np.bincount(bins, weights=weight*p, minlength=10)
    obs = np.bincount(bins, weights=weight*y, minlength=10)
    pclip = np.clip(p, 1e-6, 1-1e-6)
    return dict(ece=float(np.abs(pred-obs).sum()/weight.sum()),
                brier=float(np.average((p-y)**2, weights=weight)),
                log_loss=float(np.average(-y*np.log(pclip)-(1-y)*np.log1p(-pclip), weights=weight)))


def auc_interval(y, score, rng, draws):
    positive = score[y == 1]
    negative = score[y == 0]
    estimates = []
    for _ in range(draws):
        pos = rng.choice(positive, len(positive), replace=True)
        neg = np.sort(rng.choice(negative, len(negative), replace=True))
        credit = (np.searchsorted(neg, pos, side="left") + np.searchsorted(neg, pos, side="right"))/2
        estimates.append(float(np.mean(credit)/len(neg)))
    return np.quantile(estimates, [.025, .975]).tolist()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=HERE / "clinical-unified")
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--auc-bootstrap", type=int, default=2000)
    args = ap.parse_args()
    provenance = json.loads((args.scores.parent / "provenance.json").read_text())
    if provenance["n_rows"] != 4677 or len(provenance["auc_gates"]) != 24:
        raise ValueError("N6c score recovery gates not complete")
    if provenance.get("scores_sha256") != sha(args.scores):
        raise ValueError("Row scores do not match validated inference provenance")
    rows = sorted(list(csv.DictReader(args.scores.open(newline=""))), key=lambda r: r["example_id"])
    if len(rows) != 4677 or len({r["example_id"] for r in rows}) != 4677:
        raise ValueError("Expected exactly 4677 uniquely scored variants")
    if len({r["gene"] for r in rows}) != 20:
        raise ValueError("Expected exactly 20 genes")
    if args.bootstrap < 2 or args.auc_bootstrap < 2:
        raise ValueError("At least two bootstrap draws are required")
    if args.out.exists():
        raise SystemExit(f"Refuse overwrite: {args.out}")
    args.out.mkdir(parents=True)
    predictions, y, folds = crossfit(rows)
    per_gene = []
    rng = np.random.default_rng(SEED)
    for gene in sorted({r["gene"] for r in rows}):
        index = np.array([r["gene"] == gene for r in rows])
        llr = np.array([float(r["primary_llr_fp32"]) for r in rows])[index]
        label = 1-y[index]
        lo, hi = auc_interval(label, llr, rng, args.auc_bootstrap)
        per_gene.append(dict(gene=gene, n=int(index.sum()), n_path=int(y[index].sum()),
                             auroc=float(roc_auc_score(label, llr)), ci_lo=lo, ci_hi=hi,
                             score_mean=float(np.mean(llr)), score_std=float(np.std(llr)),
                             prevalence_path=float(y[index].mean())))
    aucs = np.array([r["auroc"] for r in per_gene])
    gene_mean_draws = rng.choice(aucs, (10000, len(aucs)), replace=True).mean(axis=1)
    oof = [dict(example_id=r["example_id"], gene=r["gene"], mutation=r["mutation"],
                label_pathogenic=int(y[i]), fold=int(folds[i]),
                **{name: float(p[i]) for name, p in predictions.items()}) for i, r in enumerate(rows)]
    write_csv(args.out / "oof_predictions.csv", oof)
    write_csv(args.out / "per_gene_discrimination.csv", per_gene)
    bins = []
    gene_metrics = []
    for name, p in predictions.items():
        for gene in ["ALL"] + sorted({r["gene"] for r in rows}):
            index = np.ones(len(rows), dtype=bool) if gene == "ALL" else np.array([r["gene"] == gene for r in rows])
            gene_metrics.append(dict(gene=gene, map=name, n=int(index.sum()), **metrics(p[index], y[index])))
            for b in range(10):
                mask = index & (np.minimum(9, (p*10).astype(int)) == b)
                bins.append(dict(gene=gene, map=name, bin=b, n=int(mask.sum()),
                                 mean_pred=float(p[mask].mean()) if mask.any() else None,
                                 obs_rate=float(y[mask].mean()) if mask.any() else None))
    write_csv(args.out / "reliability_bins.csv", bins)
    write_csv(args.out / "per_gene_calibration.csv", gene_metrics)
    coverage = []
    for name in ("platt", "isotonic"):
        p = predictions[name]
        order = np.argsort(-np.abs(p-.5), kind="stable")
        for c in np.arange(.1, 1.01, .1):
            selected = order[:max(1, int(np.ceil(c*len(rows))))]
            called = p[selected] > .5
            coverage.append(dict(map=name, coverage=float(c), n=len(selected),
                                 accuracy=float(np.mean(called == y[selected])),
                                 pathogenic_calls=int(called.sum()),
                                 pathogenic_precision=float(y[selected][called].mean()) if called.any() else None))
    write_csv(args.out / "selective_prediction.csv", coverage)
    dca = []
    for t in (.05,.10,.15,.20,.25,.30,.40,.50):
        for name, p in predictions.items():
            call = p >= t
            nb = np.sum(call*y)/len(y)-np.sum(call*(1-y))/len(y)*t/(1-t)
            dca.append(dict(threshold=t, map=name, net_benefit=float(nb)))
        dca.append(dict(threshold=t, map="treat_all", net_benefit=float(y.mean()-(1-y.mean())*t/(1-t))))
        dca.append(dict(threshold=t, map="treat_none", net_benefit=0.))
    write_csv(args.out / "decision_curve.csv", dca)
    bootstrap = []
    strata = [np.flatnonzero(np.array([r["gene"] == gene for r in rows]) & (y == label))
              for gene in sorted({r["gene"] for r in rows}) for label in (0,1)]
    for draw in range(args.bootstrap):
        sampled = np.concatenate([rng.choice(index, len(index), replace=True) for index in strata])
        weight = np.bincount(sampled, minlength=len(rows))
        pb, yb, fb = crossfit(rows, weight, seed=SEED+draw+1)
        for name, p in pb.items():
            bootstrap.append(dict(draw=draw, map=name, **metrics(p, yb, weight)))
        if (draw+1) % 50 == 0:
            print(f"identity-grouped refits {draw+1}/{args.bootstrap}", flush=True)
    write_csv(args.out / "identity_bootstrap.csv", bootstrap)
    summary = dict(protocol="A3-CLINICAL-UNIFIED-20260923", seed=SEED, folds=K,
                   n_rows=len(rows), n_genes=len(per_gene), n_bootstrap=args.bootstrap,
                   auc_bootstrap=args.auc_bootstrap,
                   macro_auroc=float(aucs.mean()), gene_mean_ci=np.quantile(gene_mean_draws,[.025,.975]).tolist(),
                   pooled_raw_auroc=float(roc_auc_score(1-y,[float(r["primary_llr_fp32"]) for r in rows])),
                   n_nominal_gene_ci_above_chance=sum(r["ci_lo"] > .5 for r in per_gene),
                   prevalence_pathogenic=float(y.mean()), maps={}, paired_brier_differences={},
                   input_sha256=sha(args.scores), scorer_provenance_sha256=sha(args.scores.parent / "provenance.json"),
                   analysis_sha256=sha(Path(__file__)),
                   interval_scope="Percentile refitting bootstrap, within gene and class; each original identity confined to one fold, multiplicities preserved. Not gene/position-disjoint or external validation.",
                   calibration="Unpenalized logistic MLE (BFGS); weighted distinct-score PAV with step predictions; fitting-fold standardization only.")
    for name, p in predictions.items():
        records = [r for r in bootstrap if r["map"] == name]
        summary["maps"][name] = dict(point=metrics(p,y), refit={metric: dict(
            ci=np.quantile([r[metric] for r in records],[.025,.975]).tolist(),
            bootstrap_mean=float(np.mean([r[metric] for r in records])),
            first_half_ci=np.quantile([r[metric] for r in records[:len(records)//2]],[.025,.975]).tolist(),
            second_half_ci=np.quantile([r[metric] for r in records[len(records)//2:]],[.025,.975]).tolist())
            for metric in ("ece","brier","log_loss")})
    for name in ("platt", "isotonic"):
        a = [r["brier"] for r in bootstrap if r["map"] == name]
        b = [r["brier"] for r in bootstrap if r["map"] == "gene_prevalence"]
        summary["paired_brier_differences"][name+"_minus_gene_prevalence"] = dict(
            point=metrics(predictions[name],y)["brier"]-metrics(predictions["gene_prevalence"],y)["brier"],
            ci=np.quantile(np.array(a)-np.array(b),[.025,.975]).tolist())
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
