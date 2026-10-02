import csv
import hashlib
import re
import tempfile
import unittest
from pathlib import Path


PANEL_FIELDS = {
    "figure_id", "panel_id", "role", "estimand_id", "input_id", "renderer",
    "output", "cohort_id", "n", "claim_level", "nonclaim", "status",
}
ESTIMAND_FIELDS = {
    "equation_id", "estimand", "variables", "unit", "independent_unit", "filter",
    "estimator", "uncertainty", "multiplicity", "code", "input", "output", "claim_ceiling", "input_sha256", "output_sha256",
}
INPUT_FIELDS = {"input_id", "sha256", "path", "status"}
PANEL_IDS = tuple("ABCDEFGH")
CLAIM_LEVELS = {"descriptive", "associational", "robustness", "transport", "mechanistic/causal", "deployment"}
CLAIM_RANK = {level: index for index, level in enumerate(("descriptive", "associational", "robustness", "transport", "mechanistic/causal", "deployment"))}
DECLARED_COHORTS = {"A3-frozen"}


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def validate_registry(panels, estimands, inputs):
    errors = []
    if not panels:
        return ["figure_panels registry is empty"]
    if set(panels[0]) != PANEL_FIELDS:
        errors.append("figure_panels schema mismatch")
    if not estimands or set(estimands[0]) != ESTIMAND_FIELDS:
        errors.append("figure_estimands schema mismatch")
    if not inputs or set(inputs[0]) != INPUT_FIELDS:
        errors.append("figure_inputs schema mismatch")
    equation_ids = {row.get("equation_id") for row in estimands}
    input_ids = {row.get("input_id") for row in inputs}
    equations = {row.get("equation_id"): row for row in estimands}
    grouped = {}
    for row in panels:
        grouped.setdefault(row.get("figure_id"), []).append(row)
        for field in PANEL_FIELDS:
            if not row.get(field):
                errors.append(f"{row.get('figure_id')}:{row.get('panel_id')} missing {field}")
        if row.get("estimand_id") not in equation_ids:
            errors.append(f"{row.get('figure_id')}:{row.get('panel_id')} missing equation")
        if row.get("input_id") not in input_ids:
            errors.append(f"{row.get('figure_id')}:{row.get('panel_id')} missing input")
        if row.get("claim_level") not in CLAIM_LEVELS:
            errors.append(f"{row.get('figure_id')}:{row.get('panel_id')} invalid claim level")
        elif row.get("estimand_id") in equations and equations[row["estimand_id"]]["claim_ceiling"] in CLAIM_RANK and CLAIM_RANK[row["claim_level"]] > CLAIM_RANK[equations[row["estimand_id"]]["claim_ceiling"]]:
            errors.append(f"{row.get('figure_id')}:{row.get('panel_id')} exceeds claim ceiling")
        if row.get("cohort_id") not in DECLARED_COHORTS:
            errors.append(f"{row.get('figure_id')}:{row.get('panel_id')} unregistered cohort")
        if not re.fullmatch(r"[1-9][0-9]*|NA", row.get("n", "")):
            errors.append(f"{row.get('figure_id')}:{row.get('panel_id')} invalid n")
    for figure_id, rows in grouped.items():
        ids = [row["panel_id"] for row in rows]
        expected = list(PANEL_IDS[:len(ids)])
        if len(ids) not in {2, 3, 4, 5, 6, 7, 8}:
            errors.append(f"{figure_id} panel count must be 2-8 consecutive reader panels")
        if ids != expected:
            errors.append(f"{figure_id} panel IDs must be consecutive A-H without duplicates")
        if figure_id == "F8" and any(row.get("cohort_id", "").startswith("B1-") for row in rows):
            if {row["cohort_id"] for row in rows if row["panel_id"] in "ABC"} != {"B1-182-development"}:
                errors.append("F8 A-C must use only B1-182-development")
            d_rows = [row for row in rows if row["panel_id"] == "D"]
            if len(d_rows) != 1 or d_rows[0]["cohort_id"] != "B1-33-external":
                errors.append("F8 D must use only B1-33-external")
            if any("182+33" in row["input_id"] or row["cohort_id"] == "B1-pooled" for row in rows):
                errors.append("F8 must not pool 182 and 33 systems")
    for row in estimands:
        for field in ESTIMAND_FIELDS:
            if not row.get(field):
                errors.append(f"{row.get('equation_id')} missing {field}")
    return errors


def canonical_errors(root, panels, estimands, inputs):
    errors = []
    input_map = {row["input_id"]: row for row in inputs}
    for row in panels:
        if row["status"] != "PROMOTED":
            errors.append(f"{row['figure_id']}:{row['panel_id']} not promoted")
            continue
        source = input_map.get(row["input_id"], {})
        if source.get("status") != "PROMOTED":
            errors.append(f"{row['figure_id']}:{row['panel_id']} input not promoted")
        input_path = root / source.get("path", "")
        if not input_path.is_file():
            errors.append(f"{row['figure_id']}:{row['panel_id']} input path missing")
        if not re.fullmatch(r"[0-9a-f]{64}", source.get("sha256", "")):
            errors.append(f"{row['figure_id']}:{row['panel_id']} immutable input hash required")
        elif input_path.is_file() and hashlib.sha256(input_path.read_bytes()).hexdigest() != source["sha256"]:
            errors.append(f"{row['figure_id']}:{row['panel_id']} input hash mismatch")
        if not (root / row["output"]).is_file():
            errors.append(f"{row['figure_id']}:{row['panel_id']} canonical output missing")
    return errors


class MasterFigureRegistryContract(unittest.TestCase):
    def setUp(self):
        self.root = Path(__file__).resolve().parents[2]
        registry = self.root / "registries"
        self.panels = read_csv(registry / "figure_panels.csv")
        self.estimands = read_csv(registry / "figure_estimands.csv")
        self.inputs = read_csv(registry / "figure_inputs.sha256")

    def test_planning_registry_meets_schema_and_panel_contract(self):
        self.assertEqual(validate_registry(self.panels, self.estimands, self.inputs), [])

    def test_current_canonical_audit_is_explicit_expected_red(self):
        self.assertTrue(canonical_errors(self.root, self.panels, self.estimands, self.inputs))

    def test_all_blocked_figure_fails_canonical_acceptance(self):
        self.assertTrue(canonical_errors(self.root, self.panels[:6], self.estimands, self.inputs))

    def test_promoted_existing_files_reject_wrong_content_hash(self):
        panels, inputs = [dict(self.panels[0])], [dict(self.inputs[0])]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "input.csv").write_text("real input\n", encoding="utf-8")
            (root / "output.csv").write_text("real output\n", encoding="utf-8")
            panels[0].update(status="PROMOTED", output="output.csv")
            inputs[0].update(status="PROMOTED", path="input.csv", sha256="0" * 64)
            errors = canonical_errors(root, panels, self.estimands, inputs)
        self.assertIn("F1:A input hash mismatch", errors)
        self.assertFalse(any("path missing" in error or "immutable input hash" in error for error in errors))

    def test_contract_rejects_overclaim_fake_ready_output_and_unregistered_cohort(self):
        panels = [dict(row) for row in self.panels[:6]]
        panels[0].update(claim_level="deployment", status="PROMOTED", cohort_id="UNREGISTERED", output="missing.csv")
        inputs = [dict(row) for row in self.inputs]
        inputs[0]["sha256"] = "not-a-hash"
        errors = validate_registry(panels, self.estimands, inputs) + canonical_errors(self.root, panels, self.estimands, inputs)
        for expected in ("exceeds claim ceiling", "unregistered cohort", "input not promoted", "immutable input hash required", "canonical output missing"):
            self.assertTrue(any(expected in error for error in errors), expected)

    def test_contract_rejects_missing_equation_input_renderer_cohort_and_bad_panels(self):
        panels = [{field: "x" for field in PANEL_FIELDS} for _ in range(4)]
        for index, row in enumerate(panels):
            row.update(figure_id="F1", panel_id=("A", "B", "D", "D")[index], estimand_id="EQ-1", input_id="IN-1", n="1", claim_level="descriptive")
        panels[0]["renderer"] = ""
        panels[1]["cohort_id"] = ""
        panels[2]["input_id"] = "MISSING"
        panels[3]["estimand_id"] = "MISSING"
        estimands = [{field: "x" for field in ESTIMAND_FIELDS}]
        estimands[0]["equation_id"] = "EQ-1"
        inputs = [{"input_id": "IN-1", "sha256": "PENDING", "path": "planned.csv", "status": "BLOCKED_NO_INPUT"}]
        errors = validate_registry(panels, estimands, inputs)
        self.assertTrue(any("missing renderer" in error for error in errors))
        self.assertTrue(any("missing cohort_id" in error for error in errors))
        self.assertTrue(any("missing input" in error for error in errors))
        self.assertTrue(any("missing equation" in error for error in errors))
        self.assertTrue(any("consecutive" in error for error in errors))

    def test_b1_f8_fixture_rejects_mixed_182_and_33_system_inputs(self):
        panels = [row for row in self.panels if row["figure_id"] == "F8"]
        if not panels or not any(row["cohort_id"] == "B1-182-development" for row in panels):
            self.skipTest("B1-only contract")
        self.assertEqual(validate_registry(panels, self.estimands, self.inputs), [])
        panels[3]["cohort_id"] = "B1-pooled"
        panels[3]["input_id"] = "B1-182+33"
        errors = validate_registry(panels, self.estimands, self.inputs)
        self.assertTrue(any("F8" in error and ("D must" in error or "pool" in error) for error in errors))


if __name__ == "__main__":
    unittest.main()
