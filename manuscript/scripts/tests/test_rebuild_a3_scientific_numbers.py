import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "rebuild_a3_scientific_numbers.py"
SPEC = importlib.util.spec_from_file_location("rebuild_a3_scientific_numbers", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_auroc_uses_half_credit_for_ties():
    assert MODULE.auroc([0.0, 0.5, 0.5, 1.0], [0, 0, 1, 1]) == 0.875


def test_protein_key_matches_canonical_first_two_tokens():
    assert MODULE.protein_key("P53_HUMAN_Kotler_2018") == "P53_HUMAN"


def test_paired_bootstrap_is_deterministic_and_reports_zero_span():
    first = MODULE.paired_bootstrap([-0.1, 0.1], replicates=200, seed=7)
    second = MODULE.paired_bootstrap([-0.1, 0.1], replicates=200, seed=7)
    assert first == second
    assert first["estimate"] == 0.0
    assert first["spans_zero"] is True
