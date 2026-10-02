"""Current author-review contract, not the untouched historical submission."""
import csv
import hashlib
import json
from pathlib import Path
import re

import fitz

CAP=Path(__file__).resolve().parents[2]
MAN=CAP/"manuscript"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_current_source_and_two_pdfs_bound_without_ready_attestation():
    current=json.loads((MAN/"release-current.json").read_text())
    assert current==json.loads((MAN/"release-spec.json").read_text())
    assert current["scientific_ready_attested"] is False
    assert {r["path"] for r in current["documents"]}=={"manuscript.pdf","SI.pdf"}
    for r in current["build_inputs"]+current["documents"]:
        assert sha(MAN/r["path"])==r["sha256"],r["path"]
    for r in current["documents"]:
        doc=fitz.open(MAN/r["path"])
        assert len(doc)==r["pages"]
        assert "not submission-ready" in doc[0].get_text()


def test_current_figures_have_no_old_calibration_input():
    renderer=(MAN/"scripts/R/make_four_panel_all.R").read_text()
    clinical=renderer.split("# ===================== Fig 4 / F4_clinical_calibration")[1].split("# ===================== Fig 5")[0]
    assert "clinical_reliability_bins_20260923.csv" in clinical
    assert all(name not in clinical for name in ("F7_reliability.csv","F7_coverage.csv","F7_headtohead.csv"))
    assert "clinical_macro_20260923.csv" in renderer
    assert "F3_clinvar_20260923.csv" in renderer


def test_current_rebuild_checks_before_writes_and_never_packages_submission():
    build=(MAN/"scripts/rebuild_current.sh").read_text()
    assert build.index("verify_method_contract.py")<build.index("build_clinical_plotdata.py")<build.index("Rscript")
    assert build.index("latexmk -g -pdf -interaction=nonstopmode -halt-on-error SI.tex")<build.index("bind_author_review.py")
    assert "--lock-inputs" not in build
    assert "submission" not in build.lower()


def test_clinical_scoring_analysis_and_figure_inputs_same_version():
    root=CAP/"analysis/revision_20260923"
    summary=json.loads((root/"clinical-unified/summary.json").read_text())
    prov=json.loads((root/"row-scores/provenance.json").read_text())
    current=json.loads((MAN/"release-current.json").read_text())
    assert sha(root/"row-scores/scores.csv")==summary["input_sha256"]==prov["scores_sha256"]==current["clinical_input_sha256"]
    assert sha(root/"analyze_clinical.py")==summary["analysis_sha256"]
    a=list(csv.DictReader((root/"clinical-unified/per_gene_discrimination.csv").open()))
    b=list(csv.DictReader((MAN/"figures/plotdata/F3_clinvar_20260923.csv").open()))
    assert {r["gene"]:float(r["auroc"]) for r in a}=={r["gene_id"]:float(r["auroc_clinvar"]) for r in b}
    assert all(float(r["ci_lo"])<=float(r["auroc"])<=float(r["ci_hi"]) for r in a)


def test_no_latex_warning_or_missing_references():
    for name in ("manuscript","SI"):
        assert not re.search(r"Overfull|undefined|LaTeX Warning",(MAN/f"{name}.log").read_text())
