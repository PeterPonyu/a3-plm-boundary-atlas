"""Contract tests for CPU-only A3 E/F expansion payloads."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import subprocess
import shutil
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
EXPORTER = ROOT / "manuscript" / "scripts" / "export_a3_expansion_tables.py"
SPEC = importlib.util.spec_from_file_location("a3_expansion", EXPORTER)
if SPEC and SPEC.loader:
    expansion = importlib.util.module_from_spec(SPEC)
else:  # pragma: no cover - gives a useful RED failure before implementation.
    expansion = None

FIGURES = ("F1", "F3", "F4", "F5", "F7", "F8", "F11", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9")
REQUIRED_COLUMNS = {
    "logical_figure", "target_panel", "estimand_id", "status", "value", "interval_lo",
    "interval_hi", "units", "cohort_id", "input_id", "source_path", "source_sha256",
    "method", "uncertainty_method", "cluster_unit", "seed", "tie_handling",
    "multiplicity", "claim", "nonclaim", "blocked_reason", "downstream_action", "downstream_target",
}
BLOCKED = {"BLOCKED_NO_COORDINATES", "BLOCKED_NO_ROW_LEVEL", "BLOCKED_NO_RETAINED_COEFFICIENTS"}
DERIVED_TABLES = {
    "S5B_supervised_group_mean_median_abs_spearman": {
        "panel": ("S5", "B", "EQ-S5"),
        "source": "FS5_vep_comparison.csv",
        "groups": [("no", "", 4, 0.4205), ("yes", "", 3, 0.4136666666666667)],
    },
    "S6D_family_mean_protein_auroc": {
        "panel": ("S6", "D", "EQ-S6"),
        "source": "F2_crossmodel.csv",
        "groups": [("BERT/UniRef100+BFD", "", 1, 0.711), ("ESM/UR50", "", 2, 0.7265)],
    },
    "S7A_method_shift_mean_rho_unc": {
        "panel": ("S7", "A", "EQ-S7"),
        "source": "F9_rho_unc.csv",
        "groups": [
            ("BayesianRidge", "extreme", 1, 0.004), ("BayesianRidge", "moderate-shift", 2, -0.0005), ("BayesianRidge", "in-distribution", 4, 0.03575),
            ("Deep Ensemble", "extreme", 1, 0.003), ("Deep Ensemble", "moderate-shift", 2, 0.1025), ("Deep Ensemble", "in-distribution", 4, 0.13075),
            ("GP-RBF", "extreme", 1, 0.063), ("GP-RBF", "moderate-shift", 2, -0.0095), ("GP-RBF", "in-distribution", 4, -0.04525),
            ("MC-Dropout", "extreme", 1, -0.021), ("MC-Dropout", "moderate-shift", 2, 0.1935), ("MC-Dropout", "in-distribution", 4, 0.13),
        ],
    },
    "F8C_liability_pooled_counts": {
        "panel": ("F8", "C", "A3-EQ-F8-POOL-01"),
        "source": "F5_developability.csv",
        "groups": [("liability", 0, 21, 0.0), ("non-liability", 23, 42, 23 / 42)],
    },
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@unittest.skipIf(expansion is None, "exporter has not been created")
class A3ExpansionEstimandTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        assert SPEC and SPEC.loader
        SPEC.loader.exec_module(expansion)

    def export(self, output: Path, root: Path = ROOT) -> None:
        subprocess.run([sys.executable, str(EXPORTER), "--root", str(root), "--output-dir", str(output)], cwd=ROOT, check=True)

    def test_source_manifest_is_complete_and_hash_bound(self) -> None:
        manifest = expansion.source_manifest(ROOT)
        self.assertEqual(len(manifest), len(expansion.SOURCE_NAMES))
        self.assertIn("F1_scale.csv", expansion.SOURCE_NAMES)
        self.assertTrue(expansion.SOURCE_CONTRACT.is_file())
        for item in manifest:
            source = ROOT / item["path"]
            self.assertTrue(source.is_file(), item)
            self.assertEqual(item["sha256"], sha256(source), item)
            self.assertIn("rows", item)
            self.assertIn("columns", item)

    def test_payloads_are_complete_typed_and_honest_about_blocks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            self.export(output)
            input_manifest = json.loads((output / "A3_EF_input_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(len(input_manifest["sources"]), len(expansion.SOURCE_NAMES))
            self.assertEqual(input_manifest["seed"], expansion.SEED)
            for figure in FIGURES:
                csv_path = output / f"{figure}_EF.csv"
                schema_path = output / f"{figure}_EF.schema.json"
                provenance_path = output / f"{figure}_EF.provenance.json"
                self.assertTrue(csv_path.is_file(), figure)
                self.assertTrue(schema_path.is_file(), figure)
                self.assertTrue(provenance_path.is_file(), figure)
                rows = list(csv.DictReader(csv_path.open(encoding="utf-8", newline="")))
                self.assertEqual([row["target_panel"] for row in rows], ["E", "F"], figure)
                self.assertEqual(set(rows[0]), REQUIRED_COLUMNS, figure)
                for row in rows:
                    self.assertEqual(row["logical_figure"], figure)
                    self.assertIn(row["status"], {"READY", *BLOCKED})
                    self.assertRegex(row["source_sha256"], r"^[0-9a-f]{64}$")
                    self.assertTrue(row["claim"] and row["nonclaim"])
                    if row["status"] != "READY":
                        self.assertTrue(row["blocked_reason"])
                        self.assertEqual(row["value"], "")
                        self.assertEqual(row["interval_lo"], "")
                        self.assertEqual(row["interval_hi"], "")
                        self.assertIn(row["downstream_action"], {"replace_with", "merge_into", "retire"})
                        self.assertTrue(row["downstream_target"])
                schema = json.loads(schema_path.read_text(encoding="utf-8"))
                self.assertEqual(set(schema["columns"]), REQUIRED_COLUMNS)
                self.assertEqual(set(schema["allowed_statuses"]), {"READY", *BLOCKED})
                provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
                self.assertEqual(provenance["output_sha256"], sha256(csv_path))
                self.assertEqual(provenance["seed"], expansion.SEED)
                self.assertTrue(provenance["scientific_payload_sha256"])
                self.assertEqual(provenance["source_contract_sha256"], sha256(expansion.SOURCE_CONTRACT))

    def test_renderer_aggregates_are_locked_source_bound_and_registry_bound(self) -> None:
        """Task 6 must consume these exact tables instead of re-aggregating raw inputs."""
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            self.export(output)
            manifest = json.loads((output / "A3_EF_input_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["output_count"], 61)
            self.assertEqual(set(manifest["derived_outputs"]), set(DERIVED_TABLES))
            registry = list(csv.DictReader((ROOT / "manuscript" / "registries" / "figure_panels.csv").open(encoding="utf-8", newline="")))
            for table_name, expected in DERIVED_TABLES.items():
                csv_path = output / f"{table_name}.csv"
                schema_path = output / f"{table_name}.schema.json"
                provenance_path = output / f"{table_name}.provenance.json"
                self.assertTrue(csv_path.is_file(), table_name)
                self.assertTrue(schema_path.is_file(), table_name)
                self.assertTrue(provenance_path.is_file(), table_name)
                rows = list(csv.DictReader(csv_path.open(encoding="utf-8", newline="")))
                observed = ([(row["group"], int(row["n_sig"]), int(row["n_datasets"]), float(row["frac"])) for row in rows]
                            if table_name == "F8C_liability_pooled_counts" else
                            [(row["group_1"], row["group_2"], int(row["n"]), float(row["value"])) for row in rows])
                self.assertEqual(len(observed), len(expected["groups"]))
                for actual, wanted in zip(observed, expected["groups"]):
                    self.assertEqual(actual[:3], wanted[:3])
                    self.assertAlmostEqual(actual[3], wanted[3], places=12)
                self.assertTrue(all(row["source_path"].endswith(expected["source"]) for row in rows))
                self.assertTrue(all(row["source_sha256"] == sha256(ROOT / row["source_path"]) for row in rows))
                self.assertEqual([int(row["sort_order"]) for row in rows], list(range(1, len(rows) + 1)))
                self.assertTrue(all(row["cohort_id"] == "A3-frozen" for row in rows))
                self.assertTrue(all(row["status"] == "READY" for row in rows))
                schema = json.loads(schema_path.read_text(encoding="utf-8"))
                self.assertEqual(schema["columns"], list(expansion.F8C_COLUMNS if table_name == "F8C_liability_pooled_counts" else expansion.DERIVED_COLUMNS))
                provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
                self.assertEqual(provenance["output_sha256"], sha256(csv_path))
                self.assertEqual(provenance["source_contract_sha256"], sha256(expansion.SOURCE_CONTRACT))
                self.assertEqual(provenance["registry_binding"], {"figure_id": expected["panel"][0], "panel_id": expected["panel"][1], "estimand_id": expected["panel"][2], "output": f"figures/plotdata/{table_name}.csv"})
                matches = [row for row in registry if (row["figure_id"], row["panel_id"]) == expected["panel"][:2]]
                self.assertEqual(len(matches), 1)
                self.assertEqual(matches[0]["estimand_id"], expected["panel"][2])
                self.assertEqual(matches[0]["output"], f"figures/plotdata/{table_name}.csv")
                if table_name == "F8C_liability_pooled_counts":
                    self.assertEqual(int(matches[0]["n"]), sum(int(row["n_datasets"]) for row in rows))

    def test_registry_bound_derived_outputs_are_materialized_and_match_clean_export(self) -> None:
        """Registry paths must resolve to the exact exporter-derived payload triples."""
        with tempfile.TemporaryDirectory() as temporary:
            clean = Path(temporary)
            self.export(clean)
            for table_name in DERIVED_TABLES:
                for suffix in (".csv", ".schema.json", ".provenance.json"):
                    canonical = ROOT / "manuscript" / "figures" / "plotdata" / f"{table_name}{suffix}"
                    generated = clean / f"{table_name}{suffix}"
                    self.assertTrue(canonical.is_file(), canonical)
                    self.assertEqual(canonical.read_bytes(), generated.read_bytes(), canonical)
                csv_path = ROOT / "manuscript" / "figures" / "plotdata" / f"{table_name}.csv"
                provenance = json.loads((ROOT / "manuscript" / "figures" / "plotdata" / f"{table_name}.provenance.json").read_text(encoding="utf-8"))
                self.assertEqual(provenance["output_sha256"], sha256(csv_path))

    def test_representative_hand_calculated_estimands_and_hard_locks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            self.export(Path(temporary))
            f1 = list(csv.DictReader((Path(temporary) / "F1_EF.csv").open()))
            self.assertAlmostEqual(float(f1[0]["value"]), -0.000153993095, places=9)
            self.assertEqual(f1[0]["cluster_unit"], "protein")
            self.assertLess(float(f1[0]["interval_lo"]), float(f1[0]["value"]))
            self.assertGreater(float(f1[0]["interval_hi"]), float(f1[0]["value"]))
            self.assertIn("bootstrap", f1[0]["uncertainty_method"])
            f3 = list(csv.DictReader((Path(temporary) / "F3_EF.csv").open()))
            self.assertEqual(f3[1]["status"], "BLOCKED_NO_COORDINATES")
            f4 = list(csv.DictReader((Path(temporary) / "F4_EF.csv").open()))
            self.assertNotIn("per-bin", " ".join(f4[0].values()).lower())
            self.assertIn("P53", f4[0]["estimand_id"])
            f7 = list(csv.DictReader((Path(temporary) / "F7_EF.csv").open()))
            self.assertEqual(f7[0]["status"], "BLOCKED_NO_ROW_LEVEL")
            self.assertIn("3,424", f7[0]["nonclaim"])
            f8 = list(csv.DictReader((Path(temporary) / "F8_EF.csv").open()))
            self.assertTrue(all(row["source_path"].endswith("F5_developability.csv") for row in f8))
            self.assertIn("alpha=0.05", f8[1]["method"])
            self.assertIn("power=0.80", f8[1]["method"])
            self.assertIn("not equivalence", f8[1]["nonclaim"])

    def test_f5_uses_stability_sources_and_f8_excludes_antiberty(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            self.export(Path(temporary))
            f5 = list(csv.DictReader((Path(temporary) / "F5_EF.csv").open()))
            self.assertEqual(f5[0]["source_path"], "manuscript/figures/plotdata/F10_discrimination.csv")
            self.assertEqual(f5[1]["source_path"], "manuscript/figures/plotdata/F10_calibration.csv")
            self.assertNotIn("AntiBERTy", " ".join(value for row in f5 for value in row.values()))

    def test_input_contract_rejects_tamper_header_value_count_and_keeps_output_empty(self) -> None:
        mutations = {
            "hash": lambda p: p.write_text(p.read_text().replace("0.725", "0.726", 1), encoding="utf-8"),
            "header": lambda p: p.write_text(p.read_text().replace("protein_mean_auroc", "wrong_header", 1), encoding="utf-8"),
            "value": lambda p: p.write_text(p.read_text().replace("0.725", "nan", 1), encoding="utf-8"),
            "count": lambda p: p.write_text(p.read_text() + "3B,5,0.725,69,63\n", encoding="utf-8"),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary) / "root"
                shutil.copytree(ROOT / "manuscript" / "figures" / "plotdata", root / "manuscript" / "figures" / "plotdata")
                target = root / "manuscript" / "figures" / "plotdata" / "F1_scale.csv"
                mutate(target)
                output = Path(temporary) / "output"
                result = subprocess.run([sys.executable, str(EXPORTER), "--root", str(root), "--output-dir", str(output)], cwd=ROOT)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(output.exists() and any(output.iterdir()))

    def test_f1_rejects_duplicate_or_inconsistent_delta_and_is_reorder_invariant(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "root"
            shutil.copytree(ROOT / "manuscript" / "figures" / "plotdata", root / "manuscript" / "figures" / "plotdata")
            paired = root / "manuscript" / "figures" / "plotdata" / "F1_scale_paired_650M_3B.csv"
            lines = paired.read_text().splitlines()
            paired.write_text("\n".join([lines[0], lines[1], *lines[1:]]) + "\n")
            with self.assertRaises(ValueError): expansion.source_manifest(root)
        original = expansion._rows(ROOT, "F1_scale_paired_650M_3B.csv")
        self.assertEqual(expansion.f1_paired_statistics(original), expansion.f1_paired_statistics(list(reversed(original))))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "root"
            shutil.copytree(ROOT / "manuscript" / "figures" / "plotdata", root / "manuscript" / "figures" / "plotdata")
            paired = root / "manuscript" / "figures" / "plotdata" / "F1_scale_paired_650M_3B.csv"
            lines = paired.read_text().splitlines()
            fields = lines[1].split(","); fields[-1] = "0.9"; lines[1] = ",".join(fields)
            paired.write_text("\n".join(lines) + "\n")
            with self.assertRaises(ValueError): expansion.source_manifest(root)

    def test_reruns_are_byte_deterministic_and_inputs_are_immutable(self) -> None:
        before = {item["path"]: sha256(ROOT / item["path"]) for item in expansion.source_manifest(ROOT)}
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            self.export(Path(first))
            self.export(Path(second))
            for figure in FIGURES:
                self.assertEqual(
                    (Path(first) / f"{figure}_EF.csv").read_bytes(),
                    (Path(second) / f"{figure}_EF.csv").read_bytes(),
                    figure,
                )
                first_provenance = json.loads((Path(first) / f"{figure}_EF.provenance.json").read_text())
                second_provenance = json.loads((Path(second) / f"{figure}_EF.provenance.json").read_text())
                self.assertEqual(first_provenance["scientific_payload_sha256"], second_provenance["scientific_payload_sha256"])
            first_names = sorted(path.name for path in Path(first).iterdir())
            second_names = sorted(path.name for path in Path(second).iterdir())
            self.assertEqual(first_names, second_names)
            self.assertEqual(len(first_names), 61)
            for name in first_names:
                self.assertEqual((Path(first) / name).read_bytes(), (Path(second) / name).read_bytes(), name)
        after = {item["path"]: sha256(ROOT / item["path"]) for item in expansion.source_manifest(ROOT)}
        self.assertEqual(before, after)

    def test_missing_or_malformed_source_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "manuscript" / "figures" / "plotdata").mkdir(parents=True)
            with self.assertRaises((FileNotFoundError, ValueError)):
                expansion.source_manifest(root)


if __name__ == "__main__":
    unittest.main()
