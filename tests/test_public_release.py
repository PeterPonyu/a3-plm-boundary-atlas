"""Portable saved-evidence contracts; no model inference or clinical fitting."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location("public_verifier", ROOT / "tools/verify_public_release.py")
assert spec and spec.loader
verify = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verify)


def test_complete_public_evidence_and_checksums():
    report = verify.verify(ROOT)
    assert report["all_pass"] is True
    assert report["rows"] == report["oof_rows"] == 4677
    assert report["genes"] == 20
    assert report["refit_metric_records"] == 5000
    assert report["n6c_gates"] == 24
    assert report["model_inference_performed"] is False
    assert report["calibrator_fits_performed"] is False
    assert report["scientific_ready_attested"] is False


def test_auc_half_credit_for_ties_and_order_invariance():
    assert verify.auc([0, 0, 1, 1], [0.0, 0.5, 0.5, 1.0]) == 0.875
    assert verify.auc([1, 0, 1, 0], [1.0, 0.5, 0.5, 0.0]) == 0.875


def test_single_class_auc_is_not_estimable():
    with pytest.raises(ValueError, match="both classes"):
        verify.auc([1, 1], [0.0, 1.0])


def test_saved_probability_metrics_do_not_fit_a_map():
    result = verify.metrics([0.25] * 4, [0, 0, 0, 1])
    assert result["ece"] == 0.0
    assert result["brier"] == 0.1875


def test_saved_refit_percentile_interpolation():
    assert verify.percentile([0.0, 1.0, 2.0, 3.0], 0.5) == 1.5
    assert verify.percentile([1.0, 1.0, 1.0], 0.025) == 1.0


def test_invalid_saved_probability_is_a_failure():
    for p in (-0.01, 1.01, float("nan")):
        with pytest.raises(ValueError, match="probability"):
            verify.metrics([p], [0])


def test_public_builder_never_calls_private_packaging_or_fitting():
    build = (ROOT / "tools/build_public.sh").read_text()
    for forbidden in ("validate_clinical.py", "analyze_clinical.py", "rescore_clinvar_candidate.py", "build_submission_pdf.sh", "--lock-inputs"):
        assert forbidden not in build
    assert build.index("verify_public_release.py") < build.index("build_clinical_plotdata.py")


def test_exact_scorer_adapter_is_explicit_input_and_not_default():
    source = (ROOT / "tools/score_clinvar.py").read_text()
    assert '--clinical-csvs' in source and '--n6c-summary' in source
    assert 'source_rows()' in source and 'scorer.main()' in source
    assert 'load_scorer' not in (ROOT / "tools/verify_public_release.py").read_text()
