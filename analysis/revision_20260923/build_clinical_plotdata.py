#!/usr/bin/env python3
"""Export version-named figure inputs and TeX values from the unified analysis.

Historical observations stay untouched. The gnomAD column is explicitly a
saved-summary comparator, not a newly rescored population-control cohort.
"""
import csv
import hashlib
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
CAP = HERE.parents[1]
MAN = CAP / "manuscript"
DATA = HERE / "clinical-unified"
PD = MAN / "figures/plotdata"


def read(path):
    return list(csv.DictReader(path.open(newline="")))


def write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    summary = json.loads((DATA / "summary.json").read_text())
    scores = HERE / "row-scores/scores.csv"
    if hashlib.sha256(scores.read_bytes()).hexdigest() != summary["input_sha256"]:
        raise ValueError("Unified analysis is not bound to current row scores")
    old = {r["gene"]: r for r in read(PD / "F3_clinvar.csv")}
    lengths = {r["gene"]: int(r["wt_len"]) for r in read(scores)}
    per_gene = read(DATA / "per_gene_discrimination.csv")
    rows = []
    for row in sorted(per_gene, key=lambda r: float(r["auroc"])):
        gene = row["gene"].removesuffix("_HUMAN").removesuffix("_GOF")
        rows.append(dict(gene=gene, gene_id=row["gene"], auroc_clinvar=row["auroc"],
                         n_path=row["n_path"], wt_len=lengths[row["gene"]],
                         historical_auroc_gnomad=old[gene]["auroc_gnomad"],
                         ci_lo=row["ci_lo"], ci_hi=row["ci_hi"],
                         evidence="row-rescored-20260923; gnomAD historical-summary-only"))
    write(PD / "F3_clinvar_20260923.csv", rows)
    for name in ("per_gene_calibration", "reliability_bins"):
        shutil.copy2(DATA / f"{name}.csv", PD / f"clinical_{name}_20260923.csv")
    metrics = []
    for name, values in summary["maps"].items():
        metrics.append(dict(map=name, **values["point"],
                            **{f"{metric}_{bound}": values["refit"][metric]["ci"][i]
                               for metric in ("ece", "brier", "log_loss")
                               for i, bound in enumerate(("lo", "hi"))}))
    write(PD / "clinical_metrics_20260923.csv", metrics)
    differences = [dict(map=name.removesuffix("_minus_gene_prevalence"),
                        value=values["point"], lo=values["ci"][0], hi=values["ci"][1])
                   for name, values in summary["paired_brier_differences"].items()]
    write(PD / "clinical_brier_difference_20260923.csv", differences)
    write(PD / "clinical_macro_20260923.csv", [dict(value=summary["macro_auroc"],
               lo=summary["gene_mean_ci"][0], hi=summary["gene_mean_ci"][1])])
    vals = {
        "ClinicalMacroAUC": f'{summary["macro_auroc"]:.4f}',
        "ClinicalMacroLow": f'{summary["gene_mean_ci"][0]:.4f}',
        "ClinicalMacroHigh": f'{summary["gene_mean_ci"][1]:.4f}',
    }
    for name, prefix in (("sigmoid_z", "Sigmoid"), ("platt", "Platt"),
                         ("isotonic", "Iso"), ("pooled_prevalence", "PooledPrev"),
                         ("gene_prevalence", "GenePrev")):
        for metric, suffix in (("ece", "ECE"), ("brier", "Brier"), ("log_loss", "LogLoss")):
            vals[prefix+suffix] = f'{summary["maps"][name]["point"][metric]:.4f}'
    (MAN / "text/clinical_values_20260923.tex").write_text(
        "% Generated from clinical-unified/summary.json; do not hand edit.\n"+
        "".join(f"\\newcommand{{\\{name}}}{{{value}}}\n" for name, value in vals.items()))
    table = [r"\begin{tabular}{lrrrr}", r"\toprule",
             r"Map & ECE & Brier & Log loss & Brier refit range \\", r"\midrule"]
    names = {"sigmoid_z": r"Sigmoid($z$)", "platt": "Platt", "isotonic": "Isotonic",
             "pooled_prevalence": "Pooled prevalence", "gene_prevalence": "Gene prevalence"}
    for m in metrics:
        table.append(f'{names[m["map"]]} & {m["ece"]:.4f} & {m["brier"]:.4f} & '
                     f'{m["log_loss"]:.4f} & [{m["brier_lo"]:.4f}, {m["brier_hi"]:.4f}] '+r"\\")
    table += [r"\bottomrule", r"\end{tabular}"]
    (MAN / "text/clinical_table_20260923.tex").write_text("\n".join(table)+"\n")
    print("Exported unified clinical figure inputs and source-derived TeX values")


if __name__ == "__main__":
    main()
