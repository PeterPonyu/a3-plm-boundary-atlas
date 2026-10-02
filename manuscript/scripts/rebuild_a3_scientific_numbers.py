#!/usr/bin/env python3
"""Rebuild the A3 scale-axis number artifacts from raw-score manifests.

This CPU-only script closes two provenance gaps in the manuscript:

1. the paired 3B-minus-650M protein-mean AUROC bootstrap interval; and
2. the explicit WT- versus masked-marginal scoring labels for the full panel.

AUROC is computed per assay on ``model_score_primary`` and then averaged
within protein, using the first two underscore-delimited assay-name tokens as
the protein key. This matches the canonical ``audit_raw_auroc.py`` estimator.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence


def protein_key(assay: str) -> str:
    parts = assay.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else assay


def auroc(scores: Sequence[float], labels: Sequence[int]) -> float | None:
    pairs = sorted(zip(scores, labels), key=lambda pair: pair[0])
    ranks = [0.0] * len(pairs)
    start = 0
    while start < len(pairs):
        end = start
        while end + 1 < len(pairs) and pairs[end + 1][0] == pairs[start][0]:
            end += 1
        average_rank = (start + end) / 2.0 + 1.0
        for index in range(start, end + 1):
            ranks[index] = average_rank
        start = end + 1

    n_positive = sum(label == 1 for _, label in pairs)
    n_negative = len(pairs) - n_positive
    if n_positive == 0 or n_negative == 0:
        return None
    positive_rank_sum = sum(
        ranks[index] for index, (_, label) in enumerate(pairs) if label == 1
    )
    return (
        positive_rank_sum - n_positive * (n_positive + 1) / 2.0
    ) / (n_positive * n_negative)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def evaluate_manifest(path: Path, min_per_assay: int = 20) -> dict:
    by_assay: dict[str, list[tuple[float, int]]] = defaultdict(list)
    labels_by_row: dict[tuple[str, str], int] = {}
    n_rows = 0
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            try:
                assay = row["split_group_target"]
                example_id = row["example_id"]
                score = float(row["model_score_primary"])
                label = int(float(row["label_value"]))
            except (KeyError, TypeError, ValueError):
                continue
            row_key = (assay, example_id)
            if row_key in labels_by_row:
                raise ValueError(f"duplicate manifest row key: {row_key!r}")
            labels_by_row[row_key] = label
            by_assay[assay].append((score, label))
            n_rows += 1

    assay_aurocs: dict[str, float] = {}
    for assay, values in by_assay.items():
        if len(values) < min_per_assay:
            continue
        value = auroc(
            [score for score, _ in values],
            [label for _, label in values],
        )
        if value is not None:
            assay_aurocs[assay] = value

    by_protein: dict[str, list[float]] = defaultdict(list)
    for assay, value in assay_aurocs.items():
        by_protein[protein_key(assay)].append(value)
    protein_aurocs = {
        protein: statistics.mean(values)
        for protein, values in sorted(by_protein.items())
    }
    return {
        "n_rows": n_rows,
        "n_scored_assays": len(assay_aurocs),
        "n_proteins": len(protein_aurocs),
        "n_discriminating_proteins": sum(
            value > 0.5 for value in protein_aurocs.values()
        ),
        "protein_mean_auroc": statistics.mean(protein_aurocs.values()),
        "assay_aurocs": assay_aurocs,
        "protein_aurocs": protein_aurocs,
        "labels_by_row": labels_by_row,
    }


def percentile(values: Iterable[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("percentile requires at least one value")
    position = (len(ordered) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def paired_bootstrap(
    deltas: Sequence[float], *, replicates: int, seed: int
) -> dict[str, float | int | bool]:
    if not deltas:
        raise ValueError("paired bootstrap requires at least one protein")
    generator = random.Random(seed)
    bootstrap_means = []
    for _ in range(replicates):
        bootstrap_means.append(
            statistics.mean(deltas[generator.randrange(len(deltas))] for _ in deltas)
        )
    ci_low = percentile(bootstrap_means, 0.025)
    ci_high = percentile(bootstrap_means, 0.975)
    return {
        "replicates": replicates,
        "seed": seed,
        "estimate": statistics.mean(deltas),
        "ci_low": ci_low,
        "ci_high": ci_high,
        "spans_zero": ci_low <= 0.0 <= ci_high,
    }


def public_summary(result: dict) -> dict:
    return {
        key: result[key]
        for key in (
            "n_rows",
            "n_scored_assays",
            "n_proteins",
            "n_discriminating_proteins",
            "protein_mean_auroc",
        )
    }


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subset-650m", required=True, type=Path)
    parser.add_argument("--subset-3b", required=True, type=Path)
    parser.add_argument("--full-650m-wt", required=True, type=Path)
    parser.add_argument("--full-650m-masked", required=True, type=Path)
    parser.add_argument("--full-attempted-assays", type=int, default=217)
    parser.add_argument("--bootstrap-replicates", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260823)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "figures" / "plotdata",
    )
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    subset_650m = evaluate_manifest(args.subset_650m)
    subset_3b = evaluate_manifest(args.subset_3b)
    if subset_650m["labels_by_row"] != subset_3b["labels_by_row"]:
        raise ValueError("650M and 3B subset manifests do not contain identical labelled rows")
    if subset_650m["protein_aurocs"].keys() != subset_3b["protein_aurocs"].keys():
        raise ValueError("650M and 3B subset manifests do not contain identical proteins")

    paired_rows = []
    for protein in subset_650m["protein_aurocs"]:
        value_650m = subset_650m["protein_aurocs"][protein]
        value_3b = subset_3b["protein_aurocs"][protein]
        paired_rows.append(
            {
                "protein": protein,
                "auroc_650m_masked": f"{value_650m:.12g}",
                "auroc_3b_masked": f"{value_3b:.12g}",
                "delta_3b_minus_650m": f"{value_3b - value_650m:.12g}",
            }
        )
    bootstrap = paired_bootstrap(
        [float(row["delta_3b_minus_650m"]) for row in paired_rows],
        replicates=args.bootstrap_replicates,
        seed=args.seed,
    )
    write_csv(
        args.output_dir / "F1_scale_paired_650M_3B.csv",
        ["protein", "auroc_650m_masked", "auroc_3b_masked", "delta_3b_minus_650m"],
        paired_rows,
    )
    paired_json = {
        "analysis": "paired protein bootstrap of 3B minus 650M raw-score AUROC",
        "estimator": "per-assay AUROC, averaged within protein, then equal-protein mean",
        "scoring": "masked-marginal raw model_score_primary",
        "subset_650m": public_summary(subset_650m),
        "subset_3b": public_summary(subset_3b),
        "paired_n_proteins": len(paired_rows),
        "paired_bootstrap": bootstrap,
        "sources": {
            "subset_650m": {
                "file": args.subset_650m.name,
                "sha256": sha256(args.subset_650m),
            },
            "subset_3b": {
                "file": args.subset_3b.name,
                "sha256": sha256(args.subset_3b),
            },
        },
    }
    (args.output_dir / "F1_scale_paired_650M_3B.json").write_text(
        json.dumps(paired_json, indent=2, sort_keys=True) + "\n"
    )

    full_results = {
        "WT-marginal": evaluate_manifest(args.full_650m_wt),
        "masked-marginal": evaluate_manifest(args.full_650m_masked),
    }
    if full_results["WT-marginal"]["labels_by_row"] != full_results["masked-marginal"]["labels_by_row"]:
        raise ValueError("full-panel WT and masked manifests do not contain identical labelled rows")
    full_rows = []
    for scoring, result in full_results.items():
        full_rows.append(
            {
                "model": "ESM-2 650M",
                "scoring": scoring,
                "score_field": "raw model_score_primary",
                "n_attempted_assays": args.full_attempted_assays,
                "n_scored_assays": result["n_scored_assays"],
                "n_proteins": result["n_proteins"],
                "n_variants": result["n_rows"],
                "protein_mean_auroc": f'{result["protein_mean_auroc"]:.12g}',
                "n_discriminating_proteins": result["n_discriminating_proteins"],
                "fraction_discriminating_proteins": f'{result["n_discriminating_proteins"] / result["n_proteins"]:.12g}',
            }
        )
    write_csv(
        args.output_dir / "F1_full217_scoring.csv",
        list(full_rows[0].keys()),
        full_rows,
    )
    full_json = {
        "analysis": "full ProteinGym ESM-2 650M raw-score comparison",
        "estimator": "per-assay AUROC, averaged within protein, then equal-protein mean",
        "n_attempted_assays": args.full_attempted_assays,
        "results": {
            scoring: public_summary(result) for scoring, result in full_results.items()
        },
        "sources": {
            "WT-marginal": {
                "file": args.full_650m_wt.name,
                "sha256": sha256(args.full_650m_wt),
            },
            "masked-marginal": {
                "file": args.full_650m_masked.name,
                "sha256": sha256(args.full_650m_masked),
            },
        },
    }
    (args.output_dir / "F1_full217_scoring.json").write_text(
        json.dumps(full_json, indent=2, sort_keys=True) + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
