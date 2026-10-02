"""Fail-closed panel-label geometry contract for staged PDF sidecar metadata."""
import json
import hashlib
from pathlib import Path
import fitz

POINTS_PER_MM = 72 / 25.4

def rect_gap_points(left, right):
    dx = max(left[0] - right[2], right[0] - left[2], 0)
    dy = max(left[1] - right[3], right[1] - left[3], 0)
    return (dx * dx + dy * dy) ** .5

def audit_clearance(pdf, sidecar, expected):
    if sidecar is None or expected is None:
        return {"status": "UNSUPPORTED", "errors": ["no glyph-bbox sidecar"]}
    if isinstance(sidecar, (str, Path)):
        sidecar = json.loads(Path(sidecar).read_text(encoding="utf-8"))
    required = {(row["page"], row["figure"]) for row in expected}
    observed = [(row.get("page"), row.get("figure")) for row in sidecar.get("pages", [])]
    binding_errors = []
    if sidecar.get("schema_version") != 1: binding_errors.append("unsupported schema")
    if sidecar.get("pdf_sha256") != hashlib.sha256(Path(pdf).read_bytes()).hexdigest(): binding_errors.append("stale PDF sha256")
    if sidecar.get("page_count") != fitz.open(pdf).page_count: binding_errors.append("wrong PDF page count")
    if len(observed) != len(set(observed)) or set(observed) != required: binding_errors.append("missing/extra/duplicate expected pages")
    if any(not row.get("content_present") for row in sidecar.get("pages", [])): binding_errors.append("empty sidecar page")
    if binding_errors: return {"status": "EXPECTED_RED", "errors": binding_errors}
    text = sidecar.get("text", [])
    labels = [row for row in text if row.get("role") == "panel_label"]
    # Every visible semantic text role is a physical blocker. Axis titles and
    # ticks need the 3-mm scientific-label margin; all other visible text
    # needs the 2-mm general margin.
    blockers = [row for row in text if row.get("role") != "panel_label"]
    mandatory = {"panel_label", "y_axis_title", "tick_label"}
    if not mandatory <= {row.get("role") for row in text} or not labels or not blockers or any(len(row.get("bbox", [])) != 4 for row in labels + blockers):
        return {"status": "EXPECTED_RED", "errors": ["incomplete panel/y-axis/tick bbox metadata"]}
    errors = []
    for label in labels:
        for other in blockers:
            gap_mm = rect_gap_points(label["bbox"], other["bbox"]) / POINTS_PER_MM
            threshold = 3.0 if other.get("role") in {"y_axis_title", "tick_label"} else 2.0
            if gap_mm < threshold:
                errors.append(f"{label.get('text', '?')}->{other['role']}: {gap_mm:.3f} mm")
    return {"status": "FAIL" if errors else "PASS", "errors": errors}

def regression_fixtures(root):
    if root.name == "a3-plm-boundary-atlas": return [(14, "F5 panel B/y-axis overlap"), (16, "F7 rotated-label crowding")]
    if root.name == "spgd-deconv": return [(10, "F2 gutter/header"), (16, "F6 gutter/header"), (18, "F8 gutter/header")]
    return []

def test_three_mm_and_overlap_are_failures(tmp_path):
    pdf = tmp_path / "staged.pdf"; doc = fitz.open(); doc.new_page(); doc.save(pdf)
    bound = {"schema_version": 1, "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(), "page_count": 1, "pages": [{"page": 1, "figure": "F1", "content_present": True}]}
    good = bound | {"text": [{"role": "panel_label", "text": "A", "bbox": [0, 0, 8, 8]}, {"role": "y_axis_title", "bbox": [17, 0, 25, 8]}, {"role": "tick_label", "bbox": [30, 0, 38, 8]}]}
    close = bound | {"text": good["text"][:2] + [{"role": "tick_label", "bbox": [16, 0, 24, 8]}]}
    overlap = bound | {"text": [{"role": "panel_label", "text": "A", "bbox": [0, 0, 8, 8]}, {"role": "y_axis_title", "bbox": [6, 0, 14, 8]}, {"role": "tick_label", "bbox": [30, 0, 38, 8]}]}
    staged_sidecar = tmp_path / "figure.sidecar.json"; staged_sidecar.write_text(json.dumps(good), encoding="utf-8")
    assert audit_clearance(pdf, staged_sidecar, [{"page": 1, "figure": "F1"}])["status"] == "PASS"
    assert audit_clearance(pdf, close, [{"page": 1, "figure": "F1"}])["status"] == "FAIL"  # 2.822 mm
    assert audit_clearance(pdf, overlap, [{"page": 1, "figure": "F1"}])["status"] == "FAIL"

def test_current_authority_is_not_falsely_closed_without_glyph_bboxes():
    root = Path(__file__).resolve().parents[3]
    pdf = root / "manuscript/submission/flat_Proteins/manuscript.pdf"
    for page, figure in ((14, "F5"), (16, "F7")):
        assert audit_clearance(pdf, None, [{"page": page, "figure": figure}])["status"] == "UNSUPPORTED"
    assert regression_fixtures(root) == [(14, "F5 panel B/y-axis overlap"), (16, "F7 rotated-label crowding")]

def test_red_sidecar_cannot_self_declare_coverage_or_pdf_identity(tmp_path):
    pdf = tmp_path / "staged.pdf"; doc = fitz.open(); doc.new_page(); doc.save(pdf)
    sidecar = {"schema_version": 1, "text": [{"role": "panel_label", "bbox": [0, 0, 8, 8]}, {"role": "tick_label", "bbox": [17, 0, 25, 8]}]}
    path = tmp_path / "sidecar.json"; path.write_text(json.dumps(sidecar), encoding="utf-8")
    # RED: current interface accepts an unbound sidecar and lets it choose coverage.
    assert audit_clearance(pdf, path, [{"page": 1, "figure": "F1"}])["status"] != "PASS"

def test_direct_missing_mandatory_roles_is_expected_red(tmp_path):
    pdf = tmp_path / "staged.pdf"; doc = fitz.open(); doc.new_page(); doc.save(pdf)
    bound = {"schema_version": 1, "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(), "page_count": 1, "pages": [{"page": 1, "figure": "F1", "content_present": True}]}
    panel_only = bound | {"text": [{"role": "panel_label", "bbox": [0, 0, 8, 8]}]}
    no_tick = bound | {"text": [{"role": "panel_label", "bbox": [0, 0, 8, 8]}, {"role": "y_axis_title", "bbox": [17, 0, 25, 8]}]}
    assert audit_clearance(pdf, panel_only, [{"page": 1, "figure": "F1"}])["status"] == "EXPECTED_RED"
    assert audit_clearance(pdf, no_tick, [{"page": 1, "figure": "F1"}])["status"] == "EXPECTED_RED"

def test_authority_regression_sidecars_are_sha_bound_and_fail():
    root = Path(__file__).resolve().parents[3]; pdf = root / "manuscript/submission/flat_Proteins/manuscript.pdf"
    base = {"schema_version": 1, "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(), "page_count": fitz.open(pdf).page_count}
    document = fitz.open(pdf)
    # Locate the current reader-facing Fig. 4 page from its caption rather than
    # freezing a page number from the superseded 39-page/old-label build.
    page_number = next(
        number for number, page in enumerate(document, 1) if "Fig 4." in page.get_text()
    )
    words = document[page_number - 1].get_text("words")
    panel_candidates = [word for word in words if word[4] == "B" and word[1] < document[page_number - 1].rect.height * 0.4]
    assert panel_candidates, "current Fig. 4 page has no top-level B panel tag"
    panel_b = list(min(panel_candidates, key=lambda word: (word[1], word[0]))[:4])
    # This is a deliberately synthetic expected-red fixture: bind it to the
    # current PDF/page, then place a y-axis title exactly on the panel tag.  The
    # test therefore exercises the fail-closed geometry rule without relying on
    # stale token names or coordinates from an earlier renderer.
    axis = panel_b[:]
    tick = [panel_b[2] + 40, panel_b[1], panel_b[2] + 48, panel_b[3]]
    sidecar = base | {
        "pages": [{"page": page_number, "figure": "F5", "content_present": True}],
        "text": [
            {"role": "panel_label", "text": "B", "bbox": panel_b},
            {"role": "y_axis_title", "text": "synthetic", "bbox": axis},
            {"role": "tick_label", "text": "0", "bbox": tick},
        ],
    }
    assert audit_clearance(pdf, sidecar, [{"page": page_number, "figure": "F5"}])["status"] == "FAIL"
    # No glyph sidecar is shipped for the authority PDF, so an unbound audit
    # remains explicitly unsupported rather than being presented as PASS.
    assert audit_clearance(pdf, None, [{"page": 16, "figure": "F7"}])["status"] == "UNSUPPORTED"
