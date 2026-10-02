"""Font-role and embedding contract; role evidence is staged sidecar metadata."""
import json
import hashlib
import re
import subprocess
from pathlib import Path
import fitz

def font_report_errors(report, returncode=0):
    errors = []
    if returncode: errors.append("pdffonts failed")
    if re.search(r"Type\s*3", report, re.I): errors.append("Type 3 font")
    if re.search(r"\bno\s+no\b", report): errors.append("unembedded font")
    return errors

def pdffonts_errors(pdf):
    result = subprocess.run(["pdffonts", str(pdf)], text=True, capture_output=True, check=False)
    return font_report_errors(result.stdout + result.stderr, result.returncode)

def audit_font_roles(pdf, sidecar, expected=None):
    errors = pdffonts_errors(pdf)
    if sidecar is None or expected is None:
        return {"status": "UNSUPPORTED", "errors": errors + ["no font-role sidecar or external expected coverage"]}
    if isinstance(sidecar, (str, Path)):
        sidecar = json.loads(Path(sidecar).read_text(encoding="utf-8"))
    observed = [(row.get("page"), row.get("figure")) for row in sidecar.get("pages", [])]
    required = {(row["page"], row["figure"]) for row in expected}
    if (sidecar.get("schema_version") != 1 or sidecar.get("pdf_sha256") != hashlib.sha256(Path(pdf).read_bytes()).hexdigest() or sidecar.get("page_count") != fitz.open(pdf).page_count or set(observed) != required or len(observed) != len(set(observed)) or any(not row.get("content_present") for row in sidecar.get("pages", []))):
        return {"status": "EXPECTED_RED", "errors": errors + ["invalid sidecar PDF binding or externally required coverage"]}
    by_role = {}
    for item in sidecar.get("text", []):
        if not isinstance(item.get("size_pt"), (int, float)) or not item.get("role"):
            errors.append("missing text role or point size"); continue
        by_role.setdefault(item["role"], []).append(item["size_pt"])
        if item.get("bold") and item["role"] != "panel_label": errors.append(f"bold {item['role']}")
    if not by_role: errors.append("empty font-role metadata")
    for role, sizes in by_role.items():
        if max(sizes) - min(sizes) > .5: errors.append(f"{role} size drift >0.5 pt")
    return {"status": "FAIL" if errors else "PASS", "errors": errors}

def authority_pdf(root):
    paths = {"a3-plm-boundary-atlas": "manuscript/submission/flat_Proteins/manuscript.pdf", "b1-cofolding-confidence": "manuscript/submission/flat_JCIM/manuscript.pdf", "spgd-deconv": "manuscript/submission/portal_CBC/main.pdf"}
    return root / paths[root.name]

def test_font_roles_reject_size_drift_and_non_label_bold(tmp_path):
    # The real PDF's embedding is separately audited; these synthetic sidecars exercise role rules.
    pdf = authority_pdf(Path(__file__).resolve().parents[3])
    good = {"text": [{"role": "panel_label", "size_pt": 10, "bold": True}, {"role": "panel_label", "size_pt": 10.5, "bold": True}, {"role": "tick_label", "size_pt": 8, "bold": False}]}
    drift = {"text": [{"role": "tick_label", "size_pt": 7.5, "bold": False}, {"role": "tick_label", "size_pt": 8.1, "bold": False}]}
    bold_axis = {"text": [{"role": "y_axis_title", "size_pt": 8.5, "bold": True}]}
    staged_sidecar = tmp_path / "figure.sidecar.json"; staged_sidecar.write_text(json.dumps(good), encoding="utf-8")
    bound = {"schema_version": 1, "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(), "page_count": fitz.open(pdf).page_count, "pages": [{"page": 1, "figure": "F1", "content_present": True}]}
    staged_sidecar.write_text(json.dumps(bound | good), encoding="utf-8")
    assert audit_font_roles(pdf, staged_sidecar, [{"page": 1, "figure": "F1"}])["status"] == "PASS"
    assert any("size drift" in e for e in audit_font_roles(pdf, bound | drift, [{"page": 1, "figure": "F1"}])["errors"])
    assert any("bold y_axis_title" in e for e in audit_font_roles(pdf, bound | bold_axis, [{"page": 1, "figure": "F1"}])["errors"])

def test_pdffonts_contract_rejects_type3_and_unembedded_output():
    assert font_report_errors("Demo Type 3 Custom yes yes yes") == ["Type 3 font"]
    assert font_report_errors("Demo Type 1 Custom no no yes") == ["unembedded font"]

def test_authority_font_audit_is_expected_red_without_role_metadata():
    pdf = authority_pdf(Path(__file__).resolve().parents[3]); assert pdf.is_file()
    result = audit_font_roles(pdf, None, [{"page": 15, "figure": "F5"}])
    assert result["status"] == "UNSUPPORTED"
    assert any("no font-role" in error for error in result["errors"])
