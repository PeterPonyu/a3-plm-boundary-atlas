#!/usr/bin/env python3
"""Export deterministic, CPU-only A3 E/F scientific payloads from locked tables.

This exporter deliberately operates on summary tables only.  It never invokes a
model, reconstructs unavailable variant coordinates, or promotes a payload into
the canonical figure registry.  A ``BLOCKED_*`` row is an auditable result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import os
import statistics
import tempfile
from pathlib import Path


SEED = "20260824"
COHORT = "A3-frozen"
STATUSES = ("READY", "BLOCKED_NO_COORDINATES", "BLOCKED_NO_ROW_LEVEL", "BLOCKED_NO_RETAINED_COEFFICIENTS")
COLUMNS = (
    "logical_figure", "target_panel", "estimand_id", "status", "value", "interval_lo", "interval_hi",
    "units", "cohort_id", "input_id", "source_path", "source_sha256", "method", "uncertainty_method",
    "cluster_unit", "seed", "tie_handling", "multiplicity", "claim", "nonclaim", "blocked_reason", "downstream_action", "downstream_target",
)
DERIVED_COLUMNS = (
    "logical_figure", "target_panel", "estimand_id", "status", "group_1", "group_2", "value",
    "units", "cohort_id", "n", "source_path", "source_sha256", "method", "uncertainty_method",
    "cluster_unit", "seed", "sort_order", "claim", "nonclaim",
)
F8C_COLUMNS = ("logical_figure", "target_panel", "estimand_id", "status", "group", "n_sig", "n_datasets", "frac", "units", "cohort_id", "source_path", "source_sha256", "method", "uncertainty_method", "cluster_unit", "seed", "sort_order", "claim", "nonclaim")
PLOTDATA = "manuscript/figures/plotdata"
SOURCE_CONTRACT = Path(__file__).with_name("a3_expansion_sources.json")
SOURCE_CONTRACT_DATA = json.loads(SOURCE_CONTRACT.read_text(encoding="utf-8"))
BLOCKED_ACTIONS = {
    ("F3", "F"): ("replace_with", "F3 gene-summary sensitivity"),
    ("F5", "F"): ("replace_with", "F5 locked split-calibration sensitivity"),
    ("F7", "E"): ("replace_with", "F7 locked CI-width diagnostic"),
    ("F8", "E"): ("merge_into", "F8 developability provenance/family audit"),
    ("F8", "F"): ("retire", "F8 power/MDE until dataset-level locked results exist"),
    ("S1", "F"): ("merge_into", "S1 assay CI coverage"),
    ("S3", "F"): ("retire", "S3 influence analysis until binding rows exist"),
    ("S4", "F"): ("replace_with", "S4 locked coverage sensitivity"),
}
DERIVED_TABLES = {
    "S5B_supervised_group_mean_median_abs_spearman": {
        "figure": "S5", "panel": "B", "estimand": "EQ-S5", "source": "FS5_vep_comparison.csv",
        "units": "mean median absolute Spearman rho", "cluster": "predictor summary",
        "method": "arithmetic mean of locked median_abs_spearman values grouped by supervised",
        "group_keys": ("supervised",), "order": (("no",), ("yes",)),
        "claim": "Locked predictor summaries are aggregated descriptively by supervision group.",
        "nonclaim": "This mean of predictor medians is not a paired benchmark comparison or causal supervision effect.",
    },
    "S6D_family_mean_protein_auroc": {
        "figure": "S6", "panel": "D", "estimand": "EQ-S6", "source": "F2_crossmodel.csv",
        "units": "mean protein-level AUROC", "cluster": "model summary",
        "method": "arithmetic mean of locked protein_mean_auroc values grouped by encoder family",
        "group_keys": ("family",), "order": (("BERT/UniRef100+BFD",), ("ESM/UR50",)),
        "claim": "Locked model summaries are aggregated descriptively by encoder family.",
        "nonclaim": "Family means do not establish an encoder-family effect or independent model replication.",
    },
    "S7A_method_shift_mean_rho_unc": {
        "figure": "S7", "panel": "A", "estimand": "EQ-S7", "source": "F9_rho_unc.csv",
        "units": "mean uncertainty-signal Spearman rho", "cluster": "locked shift-summary row",
        "method": "arithmetic mean of locked rho_unc values grouped by method and shift_class",
        "group_keys": ("method", "shift_class"),
        "order": tuple((method, shift) for method in ("BayesianRidge", "Deep Ensemble", "GP-RBF", "MC-Dropout") for shift in ("extreme", "moderate-shift", "in-distribution")),
        "claim": "Locked uncertainty summaries are aggregated descriptively within method and shift class.",
        "nonclaim": "These grouped means are not independent experimental replications or a transport claim.",
    },
    "F8C_liability_pooled_counts": {
        "figure": "F8", "panel": "C", "estimand": "A3-EQ-F8-POOL-01", "source": "F5_developability.csv",
        "units": "pooled BH-significant fraction", "cluster": "prespecified property group",
        "method": "sum locked BH-significant counts and dataset counts within prespecified liability group, then divide counts",
        "group_keys": ("liability_group",), "order": (("liability",), ("non-liability",)), "aggregate": "sum_counts",
        "claim": "Locked BH counts are pooled descriptively within prespecified liability groups.",
        "nonclaim": "This is not a re-computed BH analysis, interval, or causal liability effect.",
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(root: Path, name: str) -> list[dict[str, str]]:
    path = root / PLOTDATA / name
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or not rows[0]:
        raise ValueError(f"malformed or empty source: {path}")
    return rows


SOURCE_NAMES = tuple(SOURCE_CONTRACT_DATA["sources"])


def _columns(spec: dict[str, object]) -> list[tuple[str, str]]:
    return [tuple(item.split(":", 1)) for item in str(spec["columns"]).split(",")]


def _validate_cell(value: str, kind: str, *, nullable: bool, name: str, column: str) -> None:
    if value in {"", "NA", "None", "null"}:
        if nullable:
            return
        raise ValueError(f"{name}:{column} is unexpectedly null")
    if kind == "str":
        return
    if kind == "bool":
        if value not in {"True", "False"}:
            raise ValueError(f"{name}:{column} is not boolean")
        return
    try:
        parsed = float(value)
    except ValueError as error:
        raise ValueError(f"{name}:{column} is not numeric") from error
    if not math.isfinite(parsed):
        raise ValueError(f"{name}:{column} is not finite")
    if kind == "int" and not parsed.is_integer():
        raise ValueError(f"{name}:{column} is not an integer")


def validate_source(root: Path, name: str) -> list[dict[str, str]]:
    spec = SOURCE_CONTRACT_DATA["sources"].get(name)
    if not isinstance(spec, dict):
        raise ValueError(f"uncontracted source: {name}")
    path = root / PLOTDATA / name
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = [column for column, _ in _columns(spec)]
        if reader.fieldnames != expected:
            raise ValueError(f"{name}: malformed header")
        rows = list(reader)
    if len(rows) != spec["rows"]:
        raise ValueError(f"{name}: row count mismatch")
    nullable = set(spec.get("nullable", []))
    for row in rows:
        if set(row) != set(expected) or any(value is None for value in row.values()):
            raise ValueError(f"{name}: malformed row")
        for column, kind in _columns(spec):
            _validate_cell(row[column], kind, nullable=column in nullable, name=name, column=column)
    for unique in spec.get("unique", []):
        if len({row[unique] for row in rows}) != len(rows):
            raise ValueError(f"{name}: duplicate {unique}")
    for column, values in spec.get("required_values", {}).items():
        if {row[column] for row in rows} != set(values):
            raise ValueError(f"{name}: required values mismatch")
    for pair in spec.get("required_pairs", []):
        if not any(row[pair[0]] == pair[1] and row[pair[2]] == pair[3] for row in rows):
            raise ValueError(f"{name}: required named row missing")
    if name == "F1_scale_paired_650M_3B.csv":
        for row in rows:
            expected_delta = float(row["auroc_3b_masked"]) - float(row["auroc_650m_masked"])
            if not math.isclose(float(row["delta_3b_minus_650m"]), expected_delta, abs_tol=1e-12, rel_tol=0.0):
                raise ValueError("F1 paired delta is inconsistent")
    if _sha256(path) != spec["sha256"]:
        raise ValueError(f"{name}: SHA-256 mismatch")
    return rows


def source_manifest(root: Path) -> list[dict[str, object]]:
    manifest: list[dict[str, object]] = []
    for name in SOURCE_NAMES:
        path = root / PLOTDATA / name
        rows = validate_source(root, name)
        manifest.append({"path": f"{PLOTDATA}/{name}", "sha256": _sha256(path), "rows": len(rows), "columns": list(rows[0])})
    return manifest


def _float(value: str) -> float:
    if value in {"", "NA", "None", "null"}:
        raise ValueError(f"missing numeric value: {value!r}")
    return float(value)


def _fmt(value: float | None) -> str:
    return "" if value is None or not math.isfinite(value) else f"{value:.9g}"


def _fmt_derived(value: float) -> str:
    """Retain enough decimal precision for locked mean aggregation replay."""
    if not math.isfinite(value):
        raise ValueError("derived output value must be finite")
    return f"{value:.17g}"


def _spearman(xs: list[float], ys: list[float]) -> float:
    def ranks(values: list[float]) -> list[float]:
        ordered = sorted(enumerate(values), key=lambda item: (item[1], item[0]))
        output = [0.0] * len(values)
        start = 0
        while start < len(ordered):
            stop = start + 1
            while stop < len(ordered) and ordered[stop][1] == ordered[start][1]:
                stop += 1
            rank = (start + 1 + stop) / 2.0
            for index, _ in ordered[start:stop]:
                output[index] = rank
            start = stop
        return output
    rx, ry = ranks(xs), ranks(ys)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    numerator = sum((x - mx) * (y - my) for x, y in zip(rx, ry))
    denom = math.sqrt(sum((x - mx) ** 2 for x in rx) * sum((y - my) ** 2 for y in ry))
    return numerator / denom if denom else float("nan")


def _percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    low, high = math.floor(position), math.ceil(position)
    return ordered[low] if low == high else ordered[low] + (ordered[high] - ordered[low]) * (position - low)


def _protein_bootstrap_interval(values: list[float], *, replicates: int = 2000) -> tuple[float, float]:
    """Deterministic percentile interval; each complete protein is one cluster."""
    if len(values) < 2:
        raise ValueError("protein bootstrap requires at least two complete proteins")
    rng = random.Random(int(SEED))
    n = len(values)
    means = [sum(values[rng.randrange(n)] for _ in range(n)) / n for _ in range(replicates)]
    return _percentile(means, 0.025), _percentile(means, 0.975)


def f1_paired_statistics(rows: list[dict[str, str]]) -> tuple[float, tuple[float, float]]:
    """Order-invariant paired statistic after identity/delta validation."""
    if len(rows) != 69 or len({row["protein"] for row in rows}) != 69:
        raise ValueError("F1 paired table must contain exactly 69 unique protein IDs")
    values = []
    for row in sorted(rows, key=lambda item: item["protein"]):
        delta = _float(row["delta_3b_minus_650m"])
        expected = _float(row["auroc_3b_masked"]) - _float(row["auroc_650m_masked"])
        if not math.isclose(delta, expected, abs_tol=1e-12, rel_tol=0.0):
            raise ValueError("F1 paired delta is inconsistent")
        values.append(delta)
    return statistics.mean(values), _protein_bootstrap_interval(values)


def _source(root: Path, name: str) -> tuple[str, str]:
    path = root / PLOTDATA / name
    _rows(root, name)  # validates before a row is emitted
    return f"{PLOTDATA}/{name}", _sha256(path)


def _derived_rows(root: Path, table_name: str, spec: dict[str, object]) -> list[dict[str, str]]:
    """Compute one renderer aggregate from a contracted source table only."""
    source_name = str(spec["source"])
    rows = _rows(root, source_name)
    source_path, source_sha = _source(root, source_name)
    group_keys = tuple(spec["group_keys"])
    grouped: dict[tuple[str, ...], list[float]] = {}
    if spec.get("aggregate") == "sum_counts":
        rows = _rows(root, str(spec["source"]))
        for row in rows:
            group = "liability" if row["property"] in {"polyreactivity", "aggregation"} else "non-liability"
            grouped.setdefault((group,), []).append((_float(row["n_sig_raw_bh"]), _float(row["n_datasets"])))
        if set(grouped) != set(spec["order"]):
            raise ValueError(f"{table_name}: unexpected grouping cardinality")
        source_path, source_sha = _source(root, str(spec["source"]))
        output = []
        for sort_order, key in enumerate(spec["order"], start=1):
            n_sig = sum(value[0] for value in grouped[key]); n_ds = sum(value[1] for value in grouped[key])
            if n_ds <= 0: raise ValueError(f"{table_name}: non-positive pooled denominator")
            output.append({"logical_figure": str(spec["figure"]), "target_panel": str(spec["panel"]), "estimand_id": str(spec["estimand"]), "status": "READY", "group": key[0], "n_sig": str(int(n_sig)), "n_datasets": str(int(n_ds)), "frac": _fmt_derived(n_sig / n_ds), "units": str(spec["units"]), "cohort_id": COHORT, "source_path": source_path, "source_sha256": source_sha, "method": str(spec["method"]), "uncertainty_method": "descriptive; not resampled", "cluster_unit": str(spec["cluster"]), "seed": SEED, "sort_order": str(sort_order), "claim": str(spec["claim"]), "nonclaim": str(spec["nonclaim"])})
        return output
    value_column = "median_abs_spearman" if source_name == "FS5_vep_comparison.csv" else "protein_mean_auroc" if source_name == "F2_crossmodel.csv" else "rho_unc"
    for row in rows:
        key = tuple(row[column] for column in group_keys)
        value = _float(row[value_column])
        if not math.isfinite(value):
            raise ValueError(f"{table_name}: non-finite source value")
        grouped.setdefault(key, []).append(value)
    order = tuple(spec["order"])
    if set(grouped) != set(order):
        raise ValueError(f"{table_name}: unexpected grouping cardinality")
    output: list[dict[str, str]] = []
    for sort_order, key in enumerate(order, start=1):
        values = grouped[key]
        if not values:
            raise ValueError(f"{table_name}: empty group")
        mean = statistics.mean(values)
        if not math.isfinite(mean):
            raise ValueError(f"{table_name}: non-finite aggregate")
        output.append({
            "logical_figure": str(spec["figure"]), "target_panel": str(spec["panel"]),
            "estimand_id": str(spec["estimand"]), "status": "READY", "group_1": key[0],
            "group_2": key[1] if len(key) > 1 else "", "value": _fmt_derived(mean), "units": str(spec["units"]),
            "cohort_id": COHORT, "n": str(len(values)), "source_path": source_path, "source_sha256": source_sha,
            "method": str(spec["method"]), "uncertainty_method": "not applicable; locked summary rows",
            "cluster_unit": str(spec["cluster"]), "seed": SEED, "sort_order": str(sort_order),
            "claim": str(spec["claim"]), "nonclaim": str(spec["nonclaim"]),
        })
    return output


def build_renderer_aggregates(root: Path) -> dict[str, list[dict[str, str]]]:
    return {name: _derived_rows(root, name, spec) for name, spec in DERIVED_TABLES.items()}


def _validate_registry_bindings(root: Path) -> None:
    registry_path = root / "manuscript" / "registries" / "figure_panels.csv"
    with registry_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for table_name, spec in DERIVED_TABLES.items():
        matches = [row for row in rows if row["figure_id"] == spec["figure"] and row["panel_id"] == spec["panel"]]
        expected_output = f"figures/plotdata/{table_name}.csv"
        if len(matches) != 1 or matches[0]["estimand_id"] != spec["estimand"] or matches[0]["output"] != expected_output:
            raise ValueError(f"{table_name}: registry binding mismatch")


def _row(root: Path, figure: str, panel: str, estimand: str, source_name: str, value: float | None, *,
         units: str, method: str, uncertainty: str, cluster: str, claim: str, nonclaim: str,
         interval: tuple[float | None, float | None] = (None, None), multiplicity: str = "not applicable",
         status: str = "READY", blocked_reason: str = "", downstream_action: str = "", downstream_target: str = "") -> dict[str, str]:
    if status not in STATUSES:
        raise ValueError(f"unsupported status: {status}")
    path, digest = _source(root, source_name)
    if status != "READY":
        if not blocked_reason:
            raise ValueError("blocked row requires a reason")
        downstream_action, downstream_target = (downstream_action, downstream_target) if downstream_action else BLOCKED_ACTIONS.get((figure, panel), ("", ""))
        if downstream_action not in {"replace_with", "merge_into", "retire"} or not downstream_target:
            raise ValueError("blocked row requires a structured downstream action")
        value, interval = None, (None, None)
    return {
        "logical_figure": figure, "target_panel": panel, "estimand_id": estimand, "status": status,
        "value": _fmt(value), "interval_lo": _fmt(interval[0]), "interval_hi": _fmt(interval[1]),
        "units": units, "cohort_id": COHORT, "input_id": f"A3-{figure}-{panel}", "source_path": path,
        "source_sha256": digest, "method": method, "uncertainty_method": uncertainty,
        "cluster_unit": cluster, "seed": SEED, "tie_handling": "average ranks; stable input order",
        "multiplicity": multiplicity, "claim": claim, "nonclaim": nonclaim, "blocked_reason": blocked_reason,
        "downstream_action": downstream_action, "downstream_target": downstream_target,
    }


def build_payloads(root: Path) -> dict[str, list[dict[str, str]]]:
    paired = _rows(root, "F1_scale_paired_650M_3B.csv")
    paired_mean, paired_interval = f1_paired_statistics(paired)
    scoring = _rows(root, "F1_full217_scoring.csv")
    scoring_by_name = {row["scoring"]: row for row in scoring if row["model"] == "ESM-2 650M"}
    scoring_delta = _float(scoring_by_name["WT-marginal"]["protein_mean_auroc"]) - _float(scoring_by_name["masked-marginal"]["protein_mean_auroc"])
    clinvar = _rows(root, "F3_clinvar.csv")
    gene_size_rho = _spearman([_float(r["wt_len"]) for r in clinvar], [_float(r["auroc_clinvar"]) for r in clinvar])
    p53 = next(r for r in _rows(root, "F7_headtohead.csv") if r["gene"] == "P53_HUMAN")
    p53_drop = _float(p53["auroc_esm2"]) - _float(p53["auroc_esm1v"])
    coverage = _rows(root, "F7_coverage.csv")
    coverage_range = max(_float(r["accuracy"]) for r in coverage) - min(_float(r["accuracy"]) for r in coverage)
    stability_disc = _rows(root, "F10_discrimination.csv")
    stability_disc_drop = statistics.mean(_float(r["spearman"]) for r in stability_disc if r["split"] == "homology-held-out") - statistics.mean(_float(r["spearman"]) for r in stability_disc if r["split"] == "train")
    stability_cal = _rows(root, "F10_calibration.csv")
    stability_cal_drop = statistics.mean(_float(r["calib_rmse"]) for r in stability_cal if r["split"] == "homology-held-out") - statistics.mean(_float(r["calib_rmse"]) for r in stability_cal if r["split"] == "train")
    binding = {r["key"]: r for r in _rows(root, "F4_binding.csv")}
    general_minus_ab = _float(binding["general_nonAB"]["spearman"]) - _float(binding["antibody_ABAG"]["spearman"])
    antibody_width = _float(binding["antibody_ABAG"]["ci_hi"]) - _float(binding["antibody_ABAG"]["ci_lo"])
    abn4 = _rows(root, "F8_antibody_n4.csv")
    effect = _float(abn4[1]["rho"]) - _float(abn4[0]["rho"])
    atlas = _rows(root, "F11_atlas_spearman.csv")
    atlas_spread = max(_float(r["value"]) for r in atlas) - min(_float(r["value"]) for r in atlas)
    fs1 = _rows(root, "FS1_proteingym_atlas.csv")
    fs1_ci_fraction = sum(int(r["ci_excludes_chance"]) for r in fs1) / len(fs1)
    fs2 = _rows(root, "FS2_clinvar_vs_gnomad.csv")
    fs2_delta = statistics.mean(_float(r["delta_gnomad_minus_clinvar"]) for r in fs2)
    fs3 = _rows(root, "FS3_antibody_robustness.csv")
    subset_pairs = {}
    for row in fs3:
        subset_pairs.setdefault(row["subset"], {})[row["method"]] = _float(row["rho"])
    fs3_delta = statistics.mean(v["ESM-2 3B (WT)"] - v["ESM-2 650M (WT)"] for v in subset_pairs.values())
    fs4 = _rows(root, "FS4_abstention_coverage.csv")
    fs4_gap = max(abs(_float(r["conf_auroc"]) - _float(r["shuffle_auroc"])) for r in fs4)
    fs5 = _rows(root, "FS5_vep_comparison.csv")
    fs5_range = max(_float(r["median_abs_spearman"]) for r in fs5) - min(_float(r["median_abs_spearman"]) for r in fs5)
    fs5_incomplete = sum(1 for r in fs5 if _float(r["n_proteins"]) < 10)
    f2 = _rows(root, "F2_crossmodel.csv")
    f2_delta = max(_float(r["protein_mean_auroc"]) for r in f2) - min(_float(r["protein_mean_auroc"]) for r in f2)
    f9 = _rows(root, "F9_summary.csv")
    f9_max = max(_float(r["mean_rho_unc_highshift"]) for r in f9)
    f9bo = _rows(root, "F9_bo_acquisition.csv")
    f9bo_delta = statistics.mean(_float(r["ucb"]) - _float(r["greedy"]) for r in f9bo)
    n9 = _rows(root, "F8_antibody_n9.csv")
    n9_gap = sum(int(r["antiberty_bh_sig"]) - int(r["esm2_bh_sig"]) for r in n9)
    predictor = _rows(root, "F10_reliability_predictor.csv")
    predictor_complete = [r for r in predictor if r["length"] not in {"", "NA"} and r["auroc"] not in {"", "NA"}]
    if len(predictor_complete) < 3:
        raise ValueError("F10 reliability predictor has too few complete length/AUROC records")
    pred_rho = _spearman([_float(r["length"]) for r in predictor_complete], [_float(r["auroc"]) for r in predictor_complete])
    pred_proxy = sum(r["length_is_proxy"].lower() == "true" for r in predictor) / len(predictor)

    return {
        "F1": [_row(root, "F1", "E", "A3-EQ-F1-PAIRED-01", "F1_scale_paired_650M_3B.csv", paired_mean, units="AUROC difference (3B−650M)", method="mean paired protein difference over 69 locked proteins", uncertainty="2,000-replicate fixed-seed protein-block percentile bootstrap", cluster="protein", interval=paired_interval, claim="The locked paired protein table has a near-zero mean 3B−650M difference.", nonclaim="This is descriptive, restricted to the 69 locked proteins, and does not establish a scale law."), _row(root, "F1", "F", "A3-EQ-F1-SCORE-01", "F1_full217_scoring.csv", scoring_delta, units="protein-mean AUROC difference", method="WT-marginal minus masked-marginal across full-217 scoring summary", uncertainty="not applicable; two locked aggregate rows", cluster="assay summary", claim="Score convention changes the reported aggregate discrimination summary.", nonclaim="This sensitivity does not establish a per-protein interval.")],
        "F3": [_row(root, "F3", "E", "A3-EQ-F3-SIZE-01", "F3_clinvar.csv", gene_size_rho, units="Spearman rho", method="Spearman correlation of gene length and ClinVar AUROC across 20 genes", uncertainty="not available from locked gene summary", cluster="gene", claim="Gene-size association is a descriptive sensitivity over the committed 20-gene table.", nonclaim="It is not a causal or deployable calibration claim."), _row(root, "F3", "F", "A3-EQ-F3-LOGO-01", "F3_clinvar.csv", None, units="AUROC", method="leave-one-gene-out calibration requires per-variant row scores", uncertainty="unavailable", cluster="gene", claim="No row-score leave-one-gene-out result is emitted.", nonclaim="Summary gene AUROCs cannot recover variant coordinates.", status="BLOCKED_NO_COORDINATES", blocked_reason="Locked F3 table contains gene-level summaries only, not per-variant predicted/observed coordinates.")],
        "F4": [_row(root, "F4", "E", "A3-EQ-F4-P53-01", "F7_headtohead.csv", p53_drop, units="AUROC difference (ESM-2−ESM-1v)", method="locked P53 overlap contrast", uncertainty="locked model-specific intervals are not paired", cluster="gene", claim="The P53 overlap summary provides a bounded head-to-head contrast.", nonclaim="This is not a full-20 reliability curve or a calibration refit."), _row(root, "F4", "F", "A3-EQ-F4-COVERAGE-01", "F7_coverage.csv", coverage_range, units="accuracy range", method="range across ten committed coverage thresholds", uncertainty="not resampled; locked curve summary", cluster="coverage threshold", claim="Accuracy varies across the committed selective-coverage thresholds.", nonclaim="Threshold sensitivity is descriptive and does not validate deployment.")],
        "F5": [_row(root, "F5", "E", "A3-EQ-F5-DISC-01", "F10_discrimination.csv", stability_disc_drop, units="Spearman difference (held-out−train)", method="mean homology-held-out minus train Spearman across locked WT/masked variants", uncertainty="locked two-variant summary; no cluster bootstrap", cluster="model variant", claim="Stability discrimination differs in the locked homology-held-out split.", nonclaim="This split summary does not establish a causal stability mechanism."), _row(root, "F5", "F", "A3-EQ-F5-CAL-01", "F10_calibration.csv", stability_cal_drop, units="calibration RMSE difference (held-out−train)", method="mean homology-held-out minus train RMSE; F10_meta/F10_reference are retained in provenance", uncertainty="locked two-variant summary; no residual analysis", cluster="model variant", claim="Stability calibration error differs between locked train and held-out summaries.", nonclaim="No residual or cluster-sensitivity claim is made without row scores/retained coefficients.")],
        "F7": [_row(root, "F7", "E", "A3-EQ-F7-SUBGROUP-01", "F4_binding.csv", None, units="Spearman difference", method="direct subgroup bootstrap requires binding-row coordinates", uncertainty="unavailable", cluster="binding example", claim="No subgroup bootstrap contrast is emitted from aggregate binding summaries.", nonclaim="The two 3,424 partitions are preserved; aggregate rows cannot be resampled as individual bindings.", status="BLOCKED_NO_ROW_LEVEL", blocked_reason="Committed binding table stores subgroup estimates and intervals, not binding-level coordinates."), _row(root, "F7", "F", "A3-EQ-F7-WIDTH-01", "F4_binding.csv", antibody_width, units="locked CI width", method="upper minus lower locked antibody–antigen interval", uncertainty="locked interval; no new bootstrap", cluster="binding subgroup", claim="The antibody–antigen estimate has a wider locked uncertainty interval than a precision-free point claim conveys.", nonclaim="Range width is an applicability diagnostic, not evidence of a positive binding effect.")],
        "F8": [_row(root, "F8", "E", "A3-EQ-F8-EFFECT-01", "F5_developability.csv", None, units="developability effect-size interval", method="property-level BH summaries lack dataset-level effects/independent-unit identifiers required for an interval", uncertainty="unavailable", cluster="dataset", claim="No developability effect-size interval is emitted from property summaries.", nonclaim="Seven property aggregates do not supply dataset-level interval evidence.", status="BLOCKED_NO_ROW_LEVEL", blocked_reason="F5_developability.csv has aggregate BH counts, not per-dataset effect estimates or independent-unit identifiers."), _row(root, "F8", "F", "A3-EQ-F8-MDE-01", "F5_developability.csv", None, units="power/MDE", method="alpha=0.05 and power=0.80 cannot be used without source-supported independent dataset effects and variance", uncertainty="unavailable", cluster="dataset", claim="No power or MDE calculation is emitted from aggregate property summaries.", nonclaim="This is not equivalence evidence; absence of discoveries cannot prove no effect.", status="BLOCKED_NO_ROW_LEVEL", blocked_reason="The committed developability table lacks independent dataset outcomes/variance for a valid alpha=0.05, power=0.80 calculation.")],
        "F11": [_row(root, "F11", "E", "A3-EQ-F11-SEPARATION-01", "F11_atlas_spearman.csv", atlas_spread, units="Spearman range", method="max-minus-min over six locked atlas Spearman axes", uncertainty="inherited only where locked intervals exist", cluster="task axis", claim="The atlas separates task-specific association magnitudes.", nonclaim="It is a synthesis of existing axes, not an independent experiment."), _row(root, "F11", "F", "A3-EQ-F11-PROVENANCE-01", "F11_atlas_auroc.csv", 1.0, units="provenance inventory fraction", method="all three committed AUROC atlas rows enumerated", uncertainty="not applicable", cluster="task axis", claim="The panel is an inventory of committed atlas provenance.", nonclaim="Inventory coverage does not establish an atlas law.")],
        "S1": [_row(root, "S1", "E", "A3-EQ-S1-CI-01", "FS1_proteingym_atlas.csv", fs1_ci_fraction, units="fraction of assays", method="fraction of 19 locked assay intervals excluding chance", uncertainty="locked assay intervals", cluster="assay", claim="ProteinGym uncertainty coverage is summarized across locked assays.", nonclaim="Assay summaries are not an external transport validation."), _row(root, "S1", "F", "A3-EQ-S1-EXPOSURE-01", "FS1_proteingym_pooled.csv", None, units="exposure audit", method="dataset/exposure sensitivity requires assay-level exposure metadata", uncertainty="unavailable", cluster="assay", claim="No exposure sensitivity is promoted.", nonclaim="The pooled summary has no row-level exposure labels.", status="BLOCKED_NO_ROW_LEVEL", blocked_reason="Committed pooled ProteinGym summary lacks exposure labels needed for an exposure sensitivity.")],
        "S2": [_row(root, "S2", "E", "A3-EQ-S2-GENE-01", "FS2_clinvar_vs_gnomad.csv", fs2_delta, units="mean AUROC difference", method="mean gnomAD-minus-ClinVar difference over 20 locked genes", uncertainty="not resampled from summary rows", cluster="gene", claim="The matching sensitivity is a descriptive gene-level contrast.", nonclaim="It does not resolve ascertainment causally."), _row(root, "S2", "F", "A3-EQ-S2-CLASS-01", "FS2_class_means.csv", _float(_rows(root, "FS2_class_means.csv")[1]["mean_clinvar_auroc"]) - _float(_rows(root, "FS2_class_means.csv")[0]["mean_clinvar_auroc"]), units="class-mean AUROC difference", method="other minus gain-of-function locked class means", uncertainty="not available", cluster="gene class", claim="Class means provide a limited ascertainment sensitivity.", nonclaim="Class labels do not identify a mechanism.")],
        "S3": [_row(root, "S3", "E", "A3-EQ-S3-PAIRED-01", "FS3_antibody_robustness.csv", fs3_delta, units="mean Spearman difference", method="mean 3B-WT minus 650M-WT difference over five locked subsets", uncertainty="locked subset summaries", cluster="subset", claim="Cross-size differences are summarized across the five committed subsets.", nonclaim="This is not protein-block bootstrap evidence."), _row(root, "S3", "F", "A3-EQ-S3-INFLUENCE-01", "FS3_antibody_robustness.csv", None, units="influence audit", method="leave-one-binding influence requires row-level data", uncertainty="unavailable", cluster="binding example", claim="No influence analysis is emitted.", nonclaim="Aggregate subset intervals cannot identify influential observations.", status="BLOCKED_NO_ROW_LEVEL", blocked_reason="No binding-row coordinates are committed in the robustness summary.")],
        "S4": [_row(root, "S4", "E", "A3-EQ-S4-RISK-01", "FS4_abstention_coverage.csv", fs4_gap, units="maximum AUROC gap", method="maximum absolute confidence-minus-shuffle AUROC over six locked coverages", uncertainty="locked summary curve", cluster="coverage threshold", claim="Selective-risk sensitivity is reported across committed thresholds.", nonclaim="The table does not establish a deployable abstention policy."), _row(root, "S4", "F", "A3-EQ-S4-THRESHOLD-01", "FS4_abstention_tiers.csv", None, units="threshold audit", method="tier-specific interval requires per-protein scores", uncertainty="unavailable", cluster="protein", claim="No new tier uncertainty interval is emitted.", nonclaim="Tier aggregates cannot be resampled as protein observations.", status="BLOCKED_NO_ROW_LEVEL", blocked_reason="Committed tier table stores aggregate AUROCs without protein-level scores.")],
        "S5": [_row(root, "S5", "E", "A3-EQ-S5-RANGE-01", "FS5_vep_comparison.csv", fs5_range, units="median absolute Spearman range", method="max-minus-min across seven locked predictors", uncertainty="not available", cluster="predictor", claim="Predictor summaries show a bounded spread in median association.", nonclaim="Comparisons use unequal task coverage and are not a universal ranking."), _row(root, "S5", "F", "A3-EQ-S5-MISSING-01", "FS5_vep_comparison.csv", fs5_incomplete, units="predictors with incomplete 10-protein coverage", method="count n_proteins below 10", uncertainty="not applicable", cluster="predictor", claim="The table records incomplete task coverage for some predictors.", nonclaim="Missingness is a provenance limitation, not a performance effect.")],
        "S6": [_row(root, "S6", "E", "A3-EQ-S6-PAIR-01", "F2_crossmodel.csv", f2_delta, units="protein-mean AUROC range", method="max-minus-min over three locked model summaries", uncertainty="not available", cluster="model", claim="Model summaries differ by a bounded aggregate range.", nonclaim="Aggregate model values do not give a paired protein uncertainty interval."), _row(root, "S6", "F", "A3-EQ-S6-DIRECTION-01", "F2_agreement.csv", _float(_rows(root, "F2_agreement.csv")[1]["spearman"]) - _float(_rows(root, "F2_agreement.csv")[0]["spearman"]), units="Spearman difference", method="WT minus masked agreement summary", uncertainty="not available", cluster="model comparison", claim="Score-direction provenance is explicitly separated between WT and masked summaries.", nonclaim="This does not identify a preferred scoring direction across tasks.")],
        "S7": [_row(root, "S7", "E", "A3-EQ-S7-SHIFT-01", "F9_summary.csv", f9_max, units="mean rho_unc", method="maximum locked high-shift uncertainty association across methods", uncertainty="locked split summaries", cluster="method", claim="The strongest committed high-shift uncertainty association is reported descriptively.", nonclaim="It is not a replication across independent experiments."), _row(root, "S7", "F", "A3-EQ-S7-NEGATIVE-01", "F9_bo_acquisition.csv", f9bo_delta, units="mean acquisition difference (UCB−greedy)", method="mean difference across seven locked splits", uncertainty="not resampled", cluster="split", claim="The acquisition comparison supplies a negative-control-style summary.", nonclaim="No claim is made that UCB generally improves acquisition.")],
        "S8": [_row(root, "S8", "E", "A3-EQ-S8-EFFECT-01", "F8_antibody_n4.csv", effect, units="Spearman difference (ESM-2−AntiBERTy)", method="difference of locked summary correlations", uncertainty="not recomputed", cluster="binding dataset", claim="The antibody-language-model contrast is bounded to the committed binding summary.", nonclaim="It does not validate developability prediction."), _row(root, "S8", "F", "A3-EQ-S8-FDR-01", "F8_antibody_n9.csv", float(n9_gap), units="BH-significant count difference", method="AntiBERTy minus ESM-2 counts over two committed property families", uncertainty="BH family inherited; no re-correction", cluster="property family", claim="The committed BH counts are retained as a family audit.", nonclaim="Counts from two property families are not a new multiplicity analysis.", multiplicity="inherited BH family from committed input")],
        "S9": [_row(root, "S9", "E", "A3-EQ-S9-PARTIAL-01", "F10_reliability_predictor.csv", pred_rho, units="Spearman rho", method=f"complete-case length versus AUROC Spearman association over {len(predictor_complete)} locked records", uncertainty="not resampled", cluster="protein-axis record", claim="Length association is reported as a descriptive proxy sensitivity.", nonclaim="Missing length records are excluded explicitly; length is a proxy and does not establish reliability causally."), _row(root, "S9", "F", "A3-EQ-S9-MISSING-01", "F10_reliability_predictor.csv", pred_proxy, units="fraction proxy-marked", method="fraction of records with length_is_proxy=True", uncertainty="not applicable", cluster="protein-axis record", claim="The provenance table quantifies the proxy annotation.", nonclaim="Proxy annotation does not resolve missingness or causal confounding.")],
    }


def _atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", delete=False, dir=path.parent) as handle:
        handle.write(content)
        temp = Path(handle.name)
    os.replace(temp, path)


def _csv_text(rows: list[dict[str, str]], columns: tuple[str, ...] = COLUMNS) -> str:
    from io import StringIO
    output = StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=columns, lineterminator="\n", extrasaction="raise")
    writer.writeheader(); writer.writerows(rows)
    return output.getvalue()


def export(root: Path, output_dir: Path) -> dict[str, str]:
    manifest = source_manifest(root)
    payloads = build_payloads(root)
    aggregates = build_renderer_aggregates(root)
    _validate_registry_bindings(root)
    if set(payloads) != {"F1", "F3", "F4", "F5", "F7", "F8", "F11", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9"}:
        raise ValueError("incomplete A3 figure payload set")
    exporter_hash = _sha256(Path(__file__))
    outputs: dict[str, str] = {}
    schema = {"schema_version": 1, "columns": {column: ("number|null" if column in {"value", "interval_lo", "interval_hi"} else "string") for column in COLUMNS}, "allowed_statuses": list(STATUSES), "required_null_behavior": "CSV empty value/interval fields represent null; BLOCKED rows require all three numeric fields empty."}
    for figure in sorted(payloads, key=lambda x: (x.startswith("S"), int(x[1:]))):
        rows = payloads[figure]
        if [row["target_panel"] for row in rows] != ["E", "F"]:
            raise ValueError(f"{figure} must contain exactly E/F rows")
        csv_path = output_dir / f"{figure}_EF.csv"
        content = _csv_text(rows)
        _atomic_write(csv_path, content)
        output_hash = _sha256(csv_path)
        source_paths = {row["source_path"] for row in rows}
        if figure == "F5":
            source_paths |= {f"{PLOTDATA}/F10_meta.csv", f"{PLOTDATA}/F10_reference.csv"}
        blocked_actions = [{"panel": row["target_panel"], "action": row["downstream_action"], "target": row["downstream_target"], "reason": row["blocked_reason"]} for row in rows if row["status"] != "READY"]
        provenance = {"schema_version": 1, "figure": figure, "exporter_sha256": exporter_hash, "source_contract_sha256": _sha256(SOURCE_CONTRACT), "command": "export_a3_expansion_tables.py --root <capsule> --output-dir <dir>", "cpu_only": True, "seed": SEED, "source_date_policy": "deterministic; no wall-clock timestamp", "sources": [item for item in manifest if item["path"] in source_paths], "output_sha256": output_hash, "scientific_payload_sha256": output_hash, "blocked_reasons": [row["blocked_reason"] for row in rows if row["blocked_reason"]], "blocked_downstream_actions": blocked_actions, "merge_or_retire_recommendation": "Do not treat E/F payloads as canonical panels until Task 6 renderer and Task 7 manuscript gates pass."}
        _atomic_write(output_dir / f"{figure}_EF.schema.json", json.dumps(schema, sort_keys=True, indent=2) + "\n")
        _atomic_write(output_dir / f"{figure}_EF.provenance.json", json.dumps(provenance, sort_keys=True, indent=2) + "\n")
        outputs[figure] = output_hash
    derived_schema = {
        "schema_version": 1,
        "columns": list(DERIVED_COLUMNS),
        "column_types": {column: ("finite-number" if column == "value" else "positive-integer" if column in {"n", "sort_order"} else "string") for column in DERIVED_COLUMNS},
        "status": "READY means only a deterministic aggregate of contracted locked summaries; it is not panel promotion.",
        "sort_order": "Rows follow the explicit table-specific order recorded in provenance.",
    }
    for table_name, rows in aggregates.items():
        spec = DERIVED_TABLES[table_name]
        csv_path = output_dir / f"{table_name}.csv"
        columns = F8C_COLUMNS if table_name == "F8C_liability_pooled_counts" else DERIVED_COLUMNS
        _atomic_write(csv_path, _csv_text(rows, columns))
        output_hash = _sha256(csv_path)
        registry_binding = {"figure_id": spec["figure"], "panel_id": spec["panel"], "estimand_id": spec["estimand"], "output": f"figures/plotdata/{table_name}.csv"}
        provenance = {
            "schema_version": 1, "table": table_name, "exporter_sha256": exporter_hash,
            "source_contract_sha256": _sha256(SOURCE_CONTRACT), "cpu_only": True, "seed": SEED,
            "definition": spec["method"], "units": spec["units"], "cohort_id": COHORT,
            "grouping": list(spec["group_keys"]), "sort_order": [list(item) for item in spec["order"]],
            "registry_binding": registry_binding,
            "sources": [item for item in manifest if item["path"] == f"{PLOTDATA}/{spec['source']}"],
            "output_sha256": output_hash, "scientific_payload_sha256": output_hash,
            "source_date_policy": "deterministic; no wall-clock timestamp",
            "claim_ceiling": "descriptive; locked derived table only, not a promoted panel",
        }
        table_schema = derived_schema if table_name != "F8C_liability_pooled_counts" else {"schema_version": 1, "columns": list(F8C_COLUMNS), "column_types": {"n_sig": "nonnegative-integer", "n_datasets": "positive-integer", "frac": "finite-number", "sort_order": "positive-integer"}, "definition": "group is prespecified liability membership; n_sig and n_datasets are explicit locked counts; frac=n_sig/n_datasets", "status": "READY means deterministic locked-summary pooling only; it is not a promoted panel."}
        _atomic_write(output_dir / f"{table_name}.schema.json", json.dumps(table_schema, sort_keys=True, indent=2) + "\n")
        _atomic_write(output_dir / f"{table_name}.provenance.json", json.dumps(provenance, sort_keys=True, indent=2) + "\n")
        outputs[table_name] = output_hash
    _atomic_write(
        output_dir / "A3_EF_input_manifest.json",
        json.dumps(
            {"schema_version": 2, "cpu_only": True, "seed": SEED, "sources": manifest,
             "derived_outputs": sorted(aggregates), "output_count": 1 + len(payloads) * 3 + len(aggregates) * 3},
            sort_keys=True, indent=2,
        ) + "\n",
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2], help="capsule root containing locked local inputs")
    parser.add_argument("--output-dir", required=True, type=Path, help="explicit destination; inputs are read-only")
    args = parser.parse_args()
    export(args.root.resolve(), args.output_dir.resolve())


if __name__ == "__main__":
    main()
