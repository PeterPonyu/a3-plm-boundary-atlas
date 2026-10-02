#!/usr/bin/env python3
"""Independently reconcile new row scores, OOF metrics, bins and intervals."""
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from sklearn.metrics import roc_auc_score, brier_score_loss, log_loss

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from analyze_clinical import crossfit


def load(path):
    return list(csv.DictReader(path.open(newline="")))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    data=HERE/"clinical-unified"
    scores=sorted(load(HERE/"row-scores/scores.csv"),key=lambda r:r["example_id"])
    prov=json.loads((HERE/"row-scores/provenance.json").read_text())
    summary=json.loads((data/"summary.json").read_text())
    oof=load(data/"oof_predictions.csv")
    draws=load(data/"identity_bootstrap.csv")
    bins=load(data/"reliability_bins.csv")
    checks={}
    checks["scores_bound_to_model_run_and_analysis"]=(sha(HERE/"row-scores/scores.csv")==prov["scores_sha256"]==summary["input_sha256"])
    checks["source_code_bound_to_analysis"]=sha(HERE/"analyze_clinical.py")==summary["analysis_sha256"]
    checks["4677_unique_gene_mutation_identities"]=(len(scores)==len({(r["gene"],r["mutation"]) for r in scores})==len({r["example_id"] for r in scores})==4677)
    checks["oof_rows_labels_same_as_scores"]=(all(a["example_id"]==b["example_id"] and int(b["label_pathogenic"])==1-int(a["label_benign"]) for a,b in zip(scores,oof)) and len(oof)==len(scores))
    y=np.array([int(r["label_pathogenic"]) for r in oof])
    checks["clinical_counts"]=(int(y.sum())==1108 and len(y)-int(y.sum())==3569)
    aurocs=[]
    genes={r["gene"] for r in scores}
    for gene in genes:
        rr=[r for r in scores if r["gene"]==gene]
        labels=[int(r["label_benign"]) for r in rr]
        primary=roc_auc_score(labels,[float(r["primary_llr_fp32"]) for r in rr])
        aurocs.append(primary)
        for mode,column in (("whole","llr_whole_fp32"),("win","llr_window_fp32")):
            key=f"{gene}:{mode}"
            if key not in prov["auc_gates"]:
                continue
            got=roc_auc_score(labels,[float(r[column]) for r in rr])
            checks[f"N6c_{key}"]=(abs(got-prov["auc_gates"][key]["N6c_rounded"])<.000051)
    checks["macro_auroc_recomputed"]=abs(np.mean(aurocs)-summary["macro_auroc"])<1e-14
    checks["twenty_repeat_forwards_exact"]=len(prov["diagnostics"])==20 and all(r["repeat_max_abs_logp_diff"]==0 for r in prov["diagnostics"].values())
    reconstructed,labels,folds=crossfit(scores)
    checks["oof_folds_and_labels_reproduced"]=np.array_equal(labels,y) and np.array_equal(folds,[int(r["fold"]) for r in oof])
    max_end_shift={}
    for name in summary["maps"]:
        p=np.array([float(r[name]) for r in oof])
        checks[f"{name}_oof_reproduced"]=np.max(np.abs(p-reconstructed[name]))<1e-12
        got={"brier":brier_score_loss(y,p),"log_loss":log_loss(y,np.clip(p,1e-6,1-1e-6))}
        for metric,value in got.items():
            checks[f"{name}_{metric}_independent"]=abs(value-summary["maps"][name]["point"][metric])<1e-12
        group=np.minimum((p*10).astype(int),9)
        ece=0
        for b in range(10):
            mask=group==b
            saved=next(r for r in bins if r["gene"]=="ALL" and r["map"]==name and int(r["bin"])==b)
            if int(mask.sum())!=int(saved["n"]):
                raise AssertionError(f"bin count differs {name} {b}")
            if mask.any():
                assert abs(float(saved["mean_pred"])-np.mean(p[mask]))<1e-12
                assert abs(float(saved["obs_rate"])-np.mean(y[mask]))<1e-12
                ece+=mask.sum()/len(y)*abs(p[mask].mean()-y[mask].mean())
        checks[f"{name}_ece_bins_independent"]=abs(ece-summary["maps"][name]["point"]["ece"])<1e-12
        dd=[r for r in draws if r["map"]==name]
        checks[f"{name}_1000_refits"]=len(dd)==len({r["draw"] for r in dd})==1000
        for metric in ("ece","brier","log_loss"):
            ci=np.quantile([float(r[metric]) for r in dd],[.025,.975])
            checks[f"{name}_{metric}_percentiles"]=np.allclose(ci,summary["maps"][name]["refit"][metric]["ci"],atol=1e-14,rtol=0)
        max_end_shift[name]=max(abs(a-b) for metric in ("ece","brier","log_loss") for a,b in zip(summary["maps"][name]["refit"][metric]["first_half_ci"],summary["maps"][name]["refit"][metric]["second_half_ci"]))
    for name in ("platt","isotonic"):
        dd={int(r["draw"]):float(r["brier"]) for r in draws if r["map"]==name}
        baseline={int(r["draw"]):float(r["brier"]) for r in draws if r["map"]=="gene_prevalence"}
        ci=np.quantile([dd[k]-baseline[k] for k in sorted(dd)],[.025,.975])
        checks[name+"_paired_brier_draws"]=np.allclose(ci,summary["paired_brier_differences"][name+"_minus_gene_prevalence"]["ci"],atol=1e-14,rtol=0)
    aicda=next(r for r in load(data/"per_gene_discrimination.csv") if r["gene"]=="AICDA_HUMAN")
    checks["AICDA_perfect_bootstrap_1_1"]=(float(aicda["ci_lo"])==float(aicda["ci_hi"])==1)
    checks={name:bool(value) for name,value in checks.items()}
    report=dict(all_pass=all(checks.values()),checks=checks,
                half_split_max_endpoint_shift_across_three_metrics=max_end_shift,
                within_panel_only=True,independent_clinical_validation=False,
                files={str(p.relative_to(HERE)):sha(p) for p in sorted(data.iterdir()) if p.is_file()})
    (HERE/"clinical-validation.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(dict(all_pass=report["all_pass"],n_checks=len(checks),failed=[k for k,v in checks.items() if not v],half_split_max_endpoint_shift=max_end_shift),indent=2))
    if not report["all_pass"]:
        raise SystemExit(1)


if __name__=="__main__":
    main()
