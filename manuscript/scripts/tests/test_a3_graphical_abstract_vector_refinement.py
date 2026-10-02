"""GA-only regression tests.

These tests deliberately do not compile the manuscript, touch method binding,
or assert anything about another figure's visual design.  Their purpose is to
keep the graphical abstract's evidence identity and disjoint file ownership
explicit while the main capsule remains concurrently edited.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import fitz
import pytest

MANUSCRIPT = Path(__file__).resolve().parents[2]
CAPSULE = MANUSCRIPT.parent
FIGURES = MANUSCRIPT / "figures"
HELPER = MANUSCRIPT / "scripts" / "make_a3_graphical_abstract_vector.py"
PROVENANCE = FIGURES / "A3_graphical_abstract.provenance.json"

spec = importlib.util.spec_from_file_location("a3_ga_vector_refinement", HELPER)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload() -> dict:
    return json.loads(PROVENANCE.read_text(encoding="utf-8"))


def test_current_ga_hashes_are_bound_to_sources_and_outputs() -> None:
    record = payload()
    for relative, expected in record["sources_sha256"].items():
        assert digest(CAPSULE / relative) == expected, relative
    expected_outputs = {
        "A3_graphical_abstract.pdf",
        "A3_graphical_abstract.png",
        "A3_graphical_abstract.svg",
        "A3_graphical_abstract.display.svg",
    }
    assert set(record["output_sha256"]) == expected_outputs
    for filename, expected in record["output_sha256"].items():
        assert digest(FIGURES / filename) == expected, filename


def test_same_score_clinical_chain_and_native_endpoint_units_are_retained() -> None:
    evidence = module.load_inputs()
    record = payload()["numerical_bindings"]
    clinical = record["clinical"]
    assert clinical["n_variants"] == 4677
    assert clinical["n_genes"] == 20
    assert clinical["macro_auroc"] == pytest.approx(0.7840631994738414, abs=1e-14)
    assert clinical["brier"] == {
        name: evidence["clinical"]["maps"][name]["point"]["brier"]
        for name in ("platt", "isotonic", "gene_prevalence")
    }
    assert record["paired_dms"]["paired_n_proteins"] == 69
    assert record["paired_dms"]["paired_bootstrap"]["estimate"] == pytest.approx(
        -0.00015399309505886248, abs=1e-14
    )
    assert float(record["archived_binding"]["general"]["rho"]) == pytest.approx(
        0.2222255919909328, abs=1e-14
    )
    assert float(record["archived_binding"]["antibody"]["rho"]) == pytest.approx(
        0.07884932747536114, abs=1e-14
    )


def test_changed_bound_input_is_rejected_before_export(tmp_path: Path) -> None:
    for relative in module.BOUND_INPUTS:
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(CAPSULE / relative, target)
    target = tmp_path / "manuscript/figures/plotdata/clinical_macro_20260923.csv"
    target.write_text(target.read_text(encoding="utf-8").replace("0.784063", "0.984063"), encoding="utf-8")
    with pytest.raises(ValueError, match="bound input changed"):
        module.load_inputs(tmp_path)


def test_only_panel_labels_are_bold_and_all_pdf_text_is_arial() -> None:
    document = fitz.open(FIGURES / "A3_graphical_abstract.pdf")
    assert len(document) == 1
    page = document[0]
    assert page.rect.width == pytest.approx(468.0, abs=0.01)
    assert page.rect.height == pytest.approx(309.6, abs=0.01)
    bold: list[str] = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line.get("spans", []):
                text = span["text"].strip()
                if not text:
                    continue
                assert "Arial" in span["font"]
                assert span["size"] >= 7.5
                if "Bold" in span["font"]:
                    bold.append(text)
    assert bold == list("ABCDE")
    module.validate_pdf(FIGURES / "A3_graphical_abstract.pdf")


def test_text_bounds_have_no_collisions_or_out_of_panel_records() -> None:
    record = payload()
    boxes = []
    for row in record["typography"]:
        box = row["bbox_pt"]
        assert 2 <= box[0] < box[2] <= 466
        assert 2 <= box[1] < box[3] <= 307.6
        if row["region"] is not None:
            x0, y0, x1, y1 = row["region"]
            assert x0 <= box[0] and y0 <= box[1]
            assert box[2] <= x1 and box[3] <= y1
        for previous in boxes:
            assert (
                min(box[2], previous[2]) - max(box[0], previous[0]) <= 0.25
                or min(box[3], previous[3]) - max(box[1], previous[1]) <= 0.25
            )
        boxes.append(box)
    checks = record["automatic_geometry_checks"]
    assert checks["text_collision_pairs"] == 0
    assert checks["out_of_panel_text"] == 0
    assert checks["non_label_bold"] == 0
    assert min(row["size_pt"] for row in record["typography"]) == 7.5
    assert checks["minimum_text_pt"] >= 7.5


def test_editable_and_outline_svg_are_native_and_disjoint() -> None:
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    editable = ET.parse(FIGURES / "A3_graphical_abstract.svg").getroot()
    outline = ET.parse(FIGURES / "A3_graphical_abstract.display.svg").getroot()
    assert editable.get("viewBox") == outline.get("viewBox") == "0 0 468 309.6"
    assert len(editable.findall(".//svg:text", namespace)) == len(payload()["typography"])
    assert not outline.findall(".//svg:text", namespace)
    assert len(outline.findall(".//svg:path", namespace)) >= 100
    assert not editable.findall(".//svg:image", namespace)
    assert not outline.findall(".//svg:image", namespace)


def svg_style(text: ET.Element) -> dict[str, str]:
    return {
        key.strip(): value.strip()
        for item in text.get("style", "").split(";") if item.strip()
        for key, value in [item.split(":", 1)]
    }


def test_editable_svg_actual_styles_have_all_five_panel_letters_bold_only() -> None:
    path = FIGURES / "A3_graphical_abstract.svg"
    texts = ET.parse(path).getroot().findall(".//{http://www.w3.org/2000/svg}text")
    bold = []
    for text in texts:
        content = "".join(text.itertext()).strip()
        style = svg_style(text)
        assert "Arial" in style["font-family"]
        assert float(style["font-size"].removesuffix("px")) >= 7.5
        weight = style.get("font-weight", "normal")
        if content in set("ABCDE"):
            assert weight == "700", (content, style)
            bold.append(content)
        else:
            assert weight in {"normal", "400"}, (content, style)
    assert bold == list("ABCDE")
    checks = module.svg_profile(path)
    assert checks["bold_texts"] == list("ABCDE")
    assert checks["non_label_bold"] == 0
    assert checks["minimum_text_pt"] >= 7.5


@pytest.mark.parametrize("defect, message", (
    ("missing_panel_weight", "panel label A.*font-weight: 700"),
    ("body_bold", "Only panel labels may be bold"),
    ("small_font", "SVG text below floor"),
))
def test_svg_export_rejects_actual_style_regressions(
    tmp_path: Path, defect: str, message: str,
) -> None:
    tree = ET.parse(FIGURES / "A3_graphical_abstract.svg")
    texts = tree.getroot().findall(".//{http://www.w3.org/2000/svg}text")
    if defect == "missing_panel_weight":
        target = next(text for text in texts if text.text == "A")
        style = svg_style(target)
        style.pop("font-weight", None)
    else:
        target = next(text for text in texts if text.text not in set("ABCDE"))
        style = svg_style(target)
        style["font-weight" if defect == "body_bold" else "font-size"] = (
            "700" if defect == "body_bold" else "7.49px"
        )
    target.set("style", "; ".join(f"{key}: {value}" for key, value in style.items()))
    path = tmp_path / "regressed.svg"
    tree.write(path, encoding="utf-8", xml_declaration=True)
    with pytest.raises(ValueError, match=message):
        module.svg_profile(path)


def test_archived_binding_uses_complementary_partition_a_without_new_claims() -> None:
    evidence = module.load_inputs()
    archive = payload()["numerical_bindings"]["archived_binding"]
    assert archive["source"] == "manuscript/figures/plotdata/FS3_antibody_robustness.csv"
    for name in ("general", "antibody"):
        assert archive[name] == evidence["binding"][name]
    assert archive["general"]["n"] == "2797"
    assert archive["antibody"]["n"] == "627"
    assert archive["alternate_archive"] == evidence["alternate_archive"]
    assert archive["subset_semantics"] == {
        "cohort": "same SKEMPI set",
        "partition": "Table 2 partition A",
        "relationship": "complementary_non_overlapping",
        "general": "General (non-Ab)",
        "antibody": "Ab-antigen",
        "n_total": 3424,
        "independent_experiments": False,
    }
    assert "archived points only" in archive["interpretation"]
    assert "row/scorer/resampling identity unresolved" in archive["interpretation"]
    assert "partitions overlap" not in archive["interpretation"]
    svg = ET.parse(FIGURES / "A3_graphical_abstract.svg").getroot()
    svg_text = "\n".join(
        "".join(text.itertext())
        for text in svg.findall(".//{http://www.w3.org/2000/svg}text")
    )
    page = fitz.open(FIGURES / "A3_graphical_abstract.pdf")[0]
    pdf_text = page.get_text()
    for text in (svg_text, pdf_text):
        for required in (
            "General (non-Ab)", "Ab-antigen", "0.222", "0.079",
            "Same SKEMPI set; archived points", "only; no model claim.",
            "Row/scorer/resampling unresolved.",
        ):
            assert required in text, required
        assert "Subsets overlap" not in text
    required_e = {
        "General (non-Ab)", "Row/scorer/resampling unresolved.",
        "Same SKEMPI set; archived points", "only; no model claim.",
    }
    spans = [
        span for block in page.get_text("dict")["blocks"]
        for line in block.get("lines", []) for span in line.get("spans", [])
        if span["text"] in required_e
    ]
    assert len(spans) == len(required_e)
    for span in spans:
        x0, y0, x1, y1 = span["bbox"]
        assert 322 <= x0 < x1 <= 449 and 176 <= y0 < y1 <= 288, span
        assert span["size"] >= 7.5 and "Bold" not in span["font"]


def test_boundary_statements_and_distinct_archived_intervals_remain_visible() -> None:
    text = fitz.open(FIGURES / "A3_graphical_abstract.pdf")[0].get_text()
    for required in (
        "LLR is a score",
        "within-panel",
        "Not an equivalence test",
        "Spearman, not AUROC",
        "unresolved",
        "Same SKEMPI set; archived points",
        "only; no model claim.",
        "germline control N/A",
        "no common reliability scale",
        "external clinical validation",
    ):
        assert required in text, required
    record = payload()
    assert record["scientific_ready_attested"] is False
    archive = record["numerical_bindings"]["archived_binding"]
    assert float(archive["antibody"]["ci_lo"]) != float(archive["alternate_archive"]["ci_lo"])
    assert any("not pooled" in item for item in record["scope"])


def test_ga_export_does_not_modify_non_ga_formal_figures(tmp_path: Path) -> None:
    formal = sorted(
        path for path in FIGURES.iterdir()
        if path.suffix in {".pdf", ".png"} and path.name != "A3_graphical_abstract.pdf"
        and path.name != "A3_graphical_abstract.png"
    )
    before = {path.relative_to(FIGURES): digest(path) for path in formal}
    output = tmp_path / "ga-only-output"
    output.mkdir()
    subprocess.run(
        ["python3", "-B", str(HELPER), "--output", str(output)],
        check=True,
        capture_output=True,
        text=True,
    )
    after = {path.relative_to(FIGURES): digest(path) for path in formal}
    assert after == before
    assert set(path.name for path in output.iterdir()) == {
        "A3_graphical_abstract.pdf",
        "A3_graphical_abstract.png",
        "A3_graphical_abstract.svg",
        "A3_graphical_abstract.display.svg",
        "A3_graphical_abstract.provenance.json",
    }


def test_repeated_ga_only_exports_are_byte_identical(tmp_path: Path) -> None:
    first, second = tmp_path / "one", tmp_path / "two"
    first.mkdir()
    second.mkdir()
    for output in (first, second):
        subprocess.run(
            ["python3", "-B", str(HELPER), "--output", str(output)],
            check=True,
            capture_output=True,
            text=True,
        )
    names = {
        "A3_graphical_abstract.pdf",
        "A3_graphical_abstract.png",
        "A3_graphical_abstract.svg",
        "A3_graphical_abstract.display.svg",
        "A3_graphical_abstract.provenance.json",
    }
    for name in names:
        assert digest(first / name) == digest(second / name) == digest(FIGURES / name)
