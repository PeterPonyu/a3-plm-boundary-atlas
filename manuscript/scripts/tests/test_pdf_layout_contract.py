"""Fail-closed page-layout contract for staged PDF sidecar metadata."""
import json
import hashlib
from pathlib import Path
import fitz

def audit_layout(pdf, sidecar, expected=None):
    if sidecar is None or expected is None:
        document = fitz.open(pdf)
        return {"status": "UNSUPPORTED", "errors": [f"no layout sidecar; PDF has {document.page_count} pages"]}
    if isinstance(sidecar, (str, Path)):
        sidecar = json.loads(Path(sidecar).read_text(encoding="utf-8"))
    observed = [(row.get("page"), row.get("figure")) for row in sidecar.get("pages", [])]
    required = {(row["page"], row["figure"]) for row in expected}
    if (sidecar.get("schema_version") != 1 or sidecar.get("pdf_sha256") != hashlib.sha256(Path(pdf).read_bytes()).hexdigest() or sidecar.get("page_count") != fitz.open(pdf).page_count or set(observed) != required or len(observed) != len(set(observed)) or any(not row.get("content_present") for row in sidecar.get("pages", []))):
        return {"status": "EXPECTED_RED", "errors": ["invalid sidecar PDF binding or externally required coverage"]}
    errors = []
    for page in sidecar.get("pages", []):
        waiver = page.get("journal_class_forced") and page.get("waiver_reason")
        if page.get("orphan_heading") or page.get("orphan_caption"): errors.append("orphan heading/caption")
        if page.get("equation_overflow"): errors.append("equation overflow")
        if page.get("blank_page"): errors.append("blank page")
        if page.get("legend_occlusion") or page.get("cropped_text"): errors.append("legend occlusion/cropped text")
        if page.get("low_resolution_raster"): errors.append("low-resolution raster")
        height, whitespace = page.get("printable_height_pt", 0), page.get("content_free_height_pt", 0)
        if height <= 0: errors.append("missing printable height")
        elif whitespace / height > .25 and not waiver: errors.append("section-transition whitespace >25%")
    return {"status": "FAIL" if errors else "PASS", "errors": errors}

def authority_pdf(root):
    paths = {"a3-plm-boundary-atlas": "manuscript/submission/flat_Proteins/manuscript.pdf", "b1-cofolding-confidence": "manuscript/submission/flat_JCIM/manuscript.pdf", "spgd-deconv": "manuscript/submission/portal_CBC/main.pdf"}
    return root / paths[root.name]

def test_layout_contract_enforces_blank_orphan_overflow_and_whitespace(tmp_path):
    pdf = tmp_path / "staged.pdf"; document = fitz.open(); document.new_page(); document.save(pdf)
    bound = {"schema_version": 1, "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(), "page_count": 1, "pages": [{"page": 1, "figure": "F1", "content_present": True}]}
    good = bound | {"pages": [bound["pages"][0] | {"printable_height_pt": 600, "content_free_height_pt": 150}]}
    bad = bound | {"pages": [bound["pages"][0] | {"printable_height_pt": 600, "content_free_height_pt": 151, "orphan_caption": True, "equation_overflow": True, "legend_occlusion": True, "low_resolution_raster": True}]}
    waiver = bound | {"pages": [bound["pages"][0] | {"printable_height_pt": 600, "content_free_height_pt": 300, "journal_class_forced": True, "waiver_reason": "journal class fresh reference page"}]}
    staged_sidecar = tmp_path / "layout.sidecar.json"; staged_sidecar.write_text(json.dumps(good), encoding="utf-8")
    assert audit_layout(pdf, staged_sidecar, [{"page": 1, "figure": "F1"}])["status"] == "PASS"
    errors = audit_layout(pdf, bad, [{"page": 1, "figure": "F1"}])["errors"]
    assert any("whitespace" in e for e in errors) and any("orphan" in e for e in errors) and any("overflow" in e for e in errors)
    assert audit_layout(pdf, waiver, [{"page": 1, "figure": "F1"}])["status"] == "PASS"

def test_authority_layout_audit_is_unsupported_without_staged_sidecar():
    pdf = authority_pdf(Path(__file__).resolve().parents[3]); assert pdf.is_file()
    assert audit_layout(pdf, None, [{"page": 15, "figure": "F5"}])["status"] == "UNSUPPORTED"
