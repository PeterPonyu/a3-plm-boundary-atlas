#!/usr/bin/env python3
"""Read-only portable validation of saved A3 evidence; no inference, fits or writes."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys

from update_public_manifest import CHECKSUMS, MANIFEST, STATUS, inventory, publication_editable

ROOT = Path(__file__).resolve().parents[1]
CLINICAL = "analysis/revision_20260923/"
MAPS = ("sigmoid_z", "platt", "isotonic", "pooled_prevalence", "gene_prevalence")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(root: Path, relative: str) -> list[dict[str, str]]:
    with (root / relative).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def document(root: Path, relative: str) -> dict:
    return json.loads((root / relative).read_text(encoding="utf-8"))


def close(a: float, b: float, tolerance: float = 1e-12) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tolerance


def percentile(values: list[float], q: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * q
    lo = math.floor(position)
    hi = min(lo + 1, len(ordered) - 1)
    return ordered[lo] * (1 - position + lo) + ordered[hi] * (position - lo)


def auc(labels: list[int], values: list[float]) -> float:
    """Mann-Whitney AUROC, half credit for score ties; no fitted estimator."""
    ordered = sorted(zip(values, labels))
    positive = sum(labels)
    negative = len(labels) - positive
    if not positive or not negative:
        raise ValueError("AUROC requires both classes")
    rank_sum = 0.0
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end][0] == ordered[start][0]:
            end += 1
        rank_sum += ((start + 1 + end) / 2) * sum(label for _, label in ordered[start:end])
        start = end
    return (rank_sum - positive * (positive + 1) / 2) / (positive * negative)


def metrics(probabilities: list[float], labels: list[int]) -> dict[str, float]:
    bins = defaultdict(list)
    for p, y in zip(probabilities, labels):
        if not math.isfinite(p) or not 0 <= p <= 1:
            raise ValueError("Invalid saved OOF probability")
        bins[min(9, int(p * 10))].append((p, y))
    n = len(labels)
    return {
        "brier": math.fsum((p - y) ** 2 for p, y in zip(probabilities, labels)) / n,
        "log_loss": math.fsum(-y * math.log(min(1 - 1e-6, max(1e-6, p))) - (1 - y) * math.log1p(-min(1 - 1e-6, max(1e-6, p))) for p, y in zip(probabilities, labels)) / n,
        "ece": math.fsum(abs(math.fsum(p - y for p, y in group)) for group in bins.values()) / n,
    }


def verify(root: Path = ROOT) -> dict:
    root = root.resolve(strict=True)
    checks: dict[str, bool] = {}

    def require(name: str, passed: bool) -> None:
        checks[name] = bool(passed)
        if not passed:
            raise ValueError(f"Public release check failed: {name}")

    manifest = document(root, MANIFEST)
    require("hold_not_scientific_ready", manifest["status"] == STATUS and manifest["scientific_ready_attested"] is False)
    require("version_specific_doi_binding", manifest["zenodo"]["version_doi"] == "10.5281/zenodo.23101189" and manifest["zenodo"]["concept_doi"] == "10.5281/zenodo.21869994")
    by_path = {item["path"]: item for item in manifest["files"]}
    require("unique_manifest_paths", len(by_path) == len(manifest["files"]))
    for relative, item in by_path.items():
        path = Path(relative)
        require(f"safe_path:{relative}", not path.is_absolute() and ".." not in path.parts and (root / path).is_file() and not (root / path).is_symlink())
        require(f"hash:{relative}", sha(root / path) == item["sha256"] and (root / path).stat().st_size == item["bytes"])
    actual = {path.relative_to(root).as_posix() for path in inventory(root)}
    require("complete_clean_inventory", actual == set(by_path) | {MANIFEST, CHECKSUMS})
    require("manifest_cardinality_and_bytes", len(by_path) == manifest["file_count_excluding_manifest_and_checksums"] and sum(r["bytes"] for r in by_path.values()) == manifest["bytes_excluding_manifest_and_checksums"])
    checksum_rows = []
    for line in (root / CHECKSUMS).read_text().splitlines():
        digest, relative = line.split("  ", 1)
        checksum_rows.append((relative, digest))
    require("checksum_complete_and_unique", len(checksum_rows) == len({r for r, _ in checksum_rows}) == len(by_path) + 1 and {r for r, _ in checksum_rows} == set(by_path) | {MANIFEST})
    require("checksum_inventory_bytes_match", all(sha(root / relative) == digest for relative, digest in checksum_rows))
    allowlist = document(root, "provenance/export-allowlist.json")
    source = {r["path"]: r for r in allowlist["source_files"]}
    require("365_explicit_source_files_present", len(source) == allowlist["source_file_count"] == 365 and set(source) <= set(by_path))
    require("source_immutable_except_publication_edits", all(publication_editable(p) or sha(root / p) == r["source_sha256"] for p, r in source.items()))
    require("original_locked_build_dependency_closure_included", allowlist["locked_input_count"] == 256 and allowlist["current_build_input_count"] == 237 and allowlist["all_locked_and_current_build_inputs_included"] is True)

    scores = sorted(rows(root, CLINICAL + "row-scores/scores.csv"), key=lambda r: r["example_id"])
    provenance = document(root, CLINICAL + "row-scores/provenance.json")
    binding = document(root, CLINICAL + "row-scores/run-binding.json")
    summary = document(root, CLINICAL + "clinical-unified/summary.json")
    score_hash = sha(root / CLINICAL / "row-scores/scores.csv")
    require("current_scores_model_run_analysis_same_hash", score_hash == provenance["scores_sha256"] == summary["input_sha256"])
    require("analysis_and_provenance_bound", sha(root / CLINICAL / "analyze_clinical.py") == summary["analysis_sha256"] and sha(root / CLINICAL / "row-scores/provenance.json") == summary["scorer_provenance_sha256"])
    require("actual_scorer_run_binding", sha(root / "analysis/revision_20260922/rescore_clinvar_candidate.py") == binding["scorer_sha256"])
    require("4677_unique_identity_rows", len(scores) == len({r["example_id"] for r in scores}) == len({(r["gene"], r["mutation"]) for r in scores}) == 4677)
    gene_rows = defaultdict(list)
    for row in scores:
        gene_rows[row["gene"]].append(row)
        mutation = re.fullmatch(r"([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])", row["mutation"])
        require(f"variant_contract:{row['example_id']}", mutation is not None and 1 <= int(mutation.group(2)) <= int(row["wt_len"]) and row["source_file"] == row["gene"] + ".csv" and row["label_benign"] in {"0", "1"} and row["wt_sha256"] == provenance["wt_sequence_sha256"][row["gene"]])
        start, end = int(row["window_start_zero_based"]), int(row["window_end_exclusive"])
        windowed = int(row["wt_len"]) > provenance["context"]
        selected = "llr_window_fp32" if windowed else "llr_whole_fp32"
        require(f"saved_window_score:{row['example_id']}", row["primary_mode"] == ("windowed" if windowed else "whole") and 0 <= start < end <= int(row["wt_len"]) and (not windowed or end - start <= 1024) and start <= int(mutation.group(2)) - 1 < end and close(float(row["primary_llr_fp32"]), float(row[selected])))
    require("20_reference_and_gene_identities", len(gene_rows) == 20 and set(gene_rows) == set(provenance["wt_sequence_sha256"]) == set(provenance["clinical_source_sha256"]))
    references = rows(root, CLINICAL + "row-scores/reference-identities.csv")
    require("reference_hash_length_source_export", len(references) == 20 and {r["gene"] for r in references} == set(gene_rows) and all(r["wt_sha256"] == provenance["wt_sequence_sha256"][r["gene"]] and r["clinical_source_sha256"] == provenance["clinical_source_sha256"][r["gene"]] and all(int(s["wt_len"]) == int(r["wt_len"]) for s in gene_rows[r["gene"]]) for r in references))
    require("no_invented_full_reference_sequences", allowlist["reference_sequence_payload"]["full_sequences_present_in_permitted_source"] is False)
    label_manifest = rows(root, "analysis/xfit_calibration_20260918/inputs/manifest_gof_650M.csv")
    labels = {r["example_id"]: r for r in label_manifest}
    require("source_label_manifest_bound", sha(root / "analysis/xfit_calibration_20260918/inputs/manifest_gof_650M.csv") == provenance["manifest_sha256"] == binding["manifest_sha256"] and len(labels) == len(label_manifest) == 4677)
    require("labels_match_historical_source_identities_not_old_scores", set(labels) == {r["example_id"] for r in scores} and all(labels[r["example_id"]]["target_id"] == r["gene"] and labels[r["example_id"]]["ligand_id"] == r["mutation"] and int(float(labels[r["example_id"]]["label_value"])) == int(r["label_benign"]) for r in scores))
    by_gene = {r["gene"]: r for r in rows(root, CLINICAL + "clinical-unified/per_gene_discrimination.csv")}
    gene_aucs = []
    for gene, group in gene_rows.items():
        value = auc([int(r["label_benign"]) for r in group], [float(r["primary_llr_fp32"]) for r in group])
        gene_aucs.append(value)
        require(f"gene_auc:{gene}", close(value, float(by_gene[gene]["auroc"])) and len(group) == int(by_gene[gene]["n"]) and sum(1 - int(r["label_benign"]) for r in group) == int(by_gene[gene]["n_path"]))
        saved_gene = document(root, CLINICAL + f"row-scores/{gene}.json")
        require(f"gene_cache_rows:{gene}", len(saved_gene["rows"]) == len(group) and {r["example_id"] for r in saved_gene["rows"]} == {r["example_id"] for r in group} and all(saved_gene["diagnostics"][k] == provenance["diagnostics"][gene][k] for k in saved_gene["diagnostics"]))
    require("24_recorded_n6c_auroc_gates", len(provenance["auc_gates"]) == 24)
    for key, gate in provenance["auc_gates"].items():
        gene, mode = key.split(":")
        group = gene_rows[gene]
        column = "llr_window_fp32" if mode == "win" else "llr_whole_fp32"
        value = auc([int(r["label_benign"]) for r in group], [float(r[column]) for r in group])
        require(f"n6c:{key}", close(value, gate["auroc"]) and abs(value - gate["N6c_rounded"]) < 0.000051)
    require("macro_pooled_auc_saved_scores", close(statistics.mean(gene_aucs), summary["macro_auroc"]) and close(auc([int(r["label_benign"]) for r in scores], [float(r["primary_llr_fp32"]) for r in scores]), summary["pooled_raw_auroc"]))
    require("run_diagnostics_not_new_forwards", len(provenance["diagnostics"]) == 20 and all(r["repeat_max_abs_logp_diff"] == 0 for r in provenance["diagnostics"].values()))
    oof = rows(root, CLINICAL + "clinical-unified/oof_predictions.csv")
    require("4677_oof_labels_and_identities", len(oof) == 4677 and all(a["example_id"] == b["example_id"] and a["gene"] == b["gene"] and a["mutation"] == b["mutation"] and 1 - int(a["label_benign"]) == int(b["label_pathogenic"]) for a, b in zip(scores, oof)))
    y = [int(r["label_pathogenic"]) for r in oof]
    require("1108_pathogenic_3569_benign_five_folds", sum(y) == 1108 and len(y) - sum(y) == 3569 and {int(r["fold"]) for r in oof} == set(range(5)))
    refits = rows(root, CLINICAL + "clinical-unified/identity_bootstrap.csv")
    require("5000_saved_refit_metric_records", len(refits) == 5000 and {(int(r["draw"]), r["map"]) for r in refits} == {(draw, name) for draw in range(1000) for name in MAPS})
    saved_bins = {(r["gene"], r["map"], int(r["bin"])): r for r in rows(root, CLINICAL + "clinical-unified/reliability_bins.csv")}
    saved_metrics = {(r["gene"], r["map"]): r for r in rows(root, CLINICAL + "clinical-unified/per_gene_calibration.csv")}
    for name in MAPS:
        probabilities = [float(r[name]) for r in oof]
        measured = metrics(probabilities, y)
        require(f"stored_oof_metrics:{name}", all(close(value, summary["maps"][name]["point"][metric]) for metric, value in measured.items()))
        for gene in ["ALL"] + sorted(gene_rows):
            selected = [i for i, r in enumerate(oof) if gene == "ALL" or r["gene"] == gene]
            vals, targets = [probabilities[i] for i in selected], [y[i] for i in selected]
            require(f"per_gene_metrics:{gene}:{name}", int(saved_metrics[gene, name]["n"]) == len(selected) and all(close(value, float(saved_metrics[gene, name][metric])) for metric, value in metrics(vals, targets).items()))
            for b in range(10):
                selected_bin = [(p, target) for p, target in zip(vals, targets) if min(9, int(p * 10)) == b]
                saved = saved_bins[gene, name, b]
                require(f"saved_bin:{gene}:{name}:{b}", len(selected_bin) == int(saved["n"]) and (not selected_bin or close(statistics.mean(p for p, _ in selected_bin), float(saved["mean_pred"])) and close(statistics.mean(target for _, target in selected_bin), float(saved["obs_rate"]))))
        stored_draws = sorted([r for r in refits if r["map"] == name], key=lambda r: int(r["draw"]))
        for metric in ("ece", "brier", "log_loss"):
            values = [float(r[metric]) for r in stored_draws]
            saved = summary["maps"][name]["refit"][metric]
            require(f"recorded_refit_percentiles:{name}:{metric}", all(close(percentile(values, q), saved["ci"][i]) for i, q in enumerate((0.025, 0.975))) and close(statistics.mean(values), saved["bootstrap_mean"]) and all(close(percentile(values[:500], q), saved["first_half_ci"][i]) and close(percentile(values[500:], q), saved["second_half_ci"][i]) for i, q in enumerate((0.025, 0.975))))
    for name in ("platt", "isotonic"):
        by_draw = {(int(r["draw"]), r["map"]): float(r["brier"]) for r in refits}
        differences = [by_draw[draw, name] - by_draw[draw, "gene_prevalence"] for draw in range(1000)]
        saved = summary["paired_brier_differences"][name + "_minus_gene_prevalence"]
        require(f"paired_brier_records:{name}", all(close(percentile(differences, q), saved["ci"][i]) for i, q in enumerate((0.025, 0.975))) and close(summary["maps"][name]["point"]["brier"] - summary["maps"]["gene_prevalence"]["point"]["brier"], saved["point"]))
    plotted = {r["gene_id"]: r for r in rows(root, "manuscript/figures/plotdata/F3_clinvar_20260923.csv")}
    old = {r["gene"]: r for r in rows(root, "manuscript/figures/plotdata/F3_clinvar.csv")}
    require("current_clinical_plot_historical_gnomad_boundary", set(plotted) == set(gene_rows) and all(close(float(r["auroc_clinvar"]), float(by_gene[g]["auroc"])) and r["evidence"] == "row-rescored-20260923; gnomAD historical-summary-only" and close(float(r["historical_auroc_gnomad"]), float(old[r["gene"]]["auroc_gnomad"])) for g, r in plotted.items()))
    macro = rows(root, "manuscript/figures/plotdata/clinical_macro_20260923.csv")[0]
    require("macro_plot_same_score_version", close(float(macro["value"]), summary["macro_auroc"]) and close(float(macro["lo"]), summary["gene_mean_ci"][0]) and close(float(macro["hi"]), summary["gene_mean_ci"][1]))
    plotted_metrics = rows(root, "manuscript/figures/plotdata/clinical_metrics_20260923.csv")
    require("metrics_plot_same_saved_points", len(plotted_metrics) == 5 and all(close(float(r[m]), summary["maps"][r["map"]]["point"][m]) for r in plotted_metrics for m in ("brier", "ece", "log_loss")))
    for name in ("per_gene_calibration", "reliability_bins"):
        require(f"clinical_plot_copy:{name}", sha(root / CLINICAL / f"clinical-unified/{name}.csv") == sha(root / f"manuscript/figures/plotdata/clinical_{name}_20260923.csv"))
    require("scope_remains_within_panel", "Not gene/position-disjoint or external validation" in summary["interval_scope"] and "unresolved" in (root / "EVIDENCE_BOUNDARIES.md").read_text())
    return {
        "all_pass": True, "n_checks": len(checks),
        "files": len(actual), "bytes": sum(path.stat().st_size for path in inventory(root)),
        "rows": len(scores), "genes": len(gene_rows), "oof_rows": len(oof),
        "refit_metric_records": len(refits), "n6c_gates": len(provenance["auc_gates"]),
        "macro_auroc": summary["macro_auroc"],
        "scientific_ready_attested": False, "model_inference_performed": False,
        "calibrator_fits_performed": False, "network_required": False,
        "full_reference_sequences_bundled": False,
        "source_capsule_binding": "historical/source-version provenance; public inventory is current",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    try:
        print(json.dumps(verify(args.root), indent=2))
    except (OSError, ValueError, KeyError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
