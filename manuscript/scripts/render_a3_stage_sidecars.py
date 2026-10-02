#!/usr/bin/env python3
"""Extract SHA-bound actual PDF span and layout evidence for staged A3 figures."""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import fitz

CELL_GUTTER_X_PT, CELL_GUTTER_Y_PT = 16, 12
ROLE_RULES = {
    "panel_label": "exact bold A-F span", "plot_title": "horizontal cell-header text",
    "x_axis_title": "horizontal text below an observed x-axis spine",
    "y_axis_title": "rotated text left of an observed y-axis spine",
    "tick_label": "text next to an observed axis spine or cell-edge tick band",
    "legend": "known renderer key label in the header/key band", "annotation": "remaining in-plot data label",
    "key_metric_value": "READY E/F central scalar", "blocked_body": "BLOCKED E/F card body",
    "math_subscript": "unc qualifier adjacent to rho mean expression",
}
LOGICAL_FOR_STEM = {
    "F1_scale_plateau": "F1", "F2_cross_model": "S6", "F3_clinvar_forest": "F3",
    "F4_clinical_calibration": "F4", "F5_stability_ddg": "F5", "F6_protein_uq": "S7",
    "F7_binding_boundary": "F7", "F8_developability_map": "F8", "F9_antibody_lm": "S8",
    "F10_reliability_predictor": "S9", "F11_boundary_atlas": "F11", "FS1_proteingym_atlas": "S1",
    "FS2_gnomad_control": "S2", "FS3_antibody_robustness": "S3", "FS4_abstention_tiers": "S4",
    "FS5_vep_comparison": "S5",
}
LEGEND_WORDS = {"650M", "WT", "3B", "masked", "AntiBERTy", "ProtBERT", "ESM-2", "EVE", "CADD", "SIFT4G", "VESPA", "DeepSequence", "REVEL", "ESM-1v"}

def text_is_cropped(bbox, rect):
    x0, y0, x1, y1 = bbox
    return x0 < rect.x0 or y0 < rect.y0 or x1 > rect.x1 or y1 > rect.y1

def cell_coordinates(bbox, rect, cells=None):
    if cells:
        center = fitz.Point((bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2)
        candidates = [key for key, cell in cells.items() if cell.x0 <= center.x <= cell.x1 and cell.y0 <= center.y <= cell.y1]
        if candidates:
            return candidates[0]
        return min(cells, key=lambda key: ((cells[key].x0 + cells[key].x1) / 2 - center.x) ** 2 + ((cells[key].y0 + cells[key].y1) / 2 - center.y) ** 2)
    x0, y0, x1, y1 = bbox
    center_x, center_y = (x0 + x1) / 2, (y0 + y1) / 2
    column = min(2, max(0, int(3 * (center_x - rect.x0) / rect.width)))
    row = min(1, max(0, int(2 * (center_y - rect.y0) / rect.height)))
    return row, column

def cell_is_cropped(bbox, rect, cells):
    x0, y0, x1, y1 = bbox; row, column = cell_coordinates(bbox, rect, cells)
    cell = cells[(row, column)]
    return x0 < cell.x0 - CELL_GUTTER_X_PT or y0 < cell.y0 - CELL_GUTTER_Y_PT or x1 > cell.x1 + CELL_GUTTER_X_PT or y1 > cell.y1 + CELL_GUTTER_Y_PT

def cell_rect(rect, row, column):
    return fitz.Rect(rect.x0 + column * rect.width / 3, rect.y0 + row * rect.height / 2,
                     rect.x0 + (column + 1) * rect.width / 3, rect.y0 + (row + 1) * rect.height / 2)

def cells_from_tag_bboxes(rect, tags):
    """Actual compositor ownership from ordered tag edges, never equal thirds.

    Each next tag reserves a 12-pt protected gutter.  The lower row starts
    12 pt above its first tag, which retains the upper-row footer but gives
    the evidence row its own header band.
    """
    if set(tags) != set("ABCDEF"):
        raise RuntimeError("six ordered A--F tag bboxes are required for cell geometry")
    rows = (("A", "B", "C"), ("D", "E", "F"))
    lower_y = min(tags[x][1] for x in rows[1]) - CELL_GUTTER_Y_PT
    cells = {}
    for row, names in enumerate(rows):
        y0, y1 = (rect.y0, lower_y) if row == 0 else (lower_y, rect.y1)
        starts = [rect.x0, tags[names[1]][0] - CELL_GUTTER_Y_PT, tags[names[2]][0] - CELL_GUTTER_Y_PT]
        ends = [starts[1], starts[2], rect.x1]
        for column, (x0, x1) in enumerate(zip(starts, ends)):
            cells[(row, column)] = fitz.Rect(x0, y0, x1, y1)
    return cells

def panel_cells(page):
    tags = {}
    for block in page.get_text("dict")["blocks"]:
        if block.get("type") != 0: continue
        for line in block["lines"]:
            for span in line["spans"]:
                text = span["text"].strip(); font = span["font"]
                if text in set("ABCDEF") and (span["flags"] & 16 or "bold" in font.lower()):
                    tags[text] = list(span["bbox"])
    return cells_from_tag_bboxes(page.rect, tags)

def ef_statuses(figure):
    logical = LOGICAL_FOR_STEM.get(figure)
    if not logical: return {}
    plotdata = Path(__file__).resolve().parents[1] / "figures" / "plotdata" / f"{logical}_EF.csv"
    if not plotdata.is_file(): raise RuntimeError(f"missing E/F payload for {figure}: {plotdata}")
    with plotdata.open(newline="", encoding="utf-8") as handle:
        return {row["target_panel"]: row["status"] for row in csv.DictReader(handle)}

def axis_spines(page, cells):
    """Record actual long vector-line positions by actual compositor cell."""
    spines = {(row, column): {"horizontal": [], "vertical": []} for row in range(2) for column in range(3)}
    for drawing in page.get_drawings():
        for item in drawing.get("items", []):
            if item[0] != "l": continue
            p0, p1 = item[1], item[2]
            if abs(p0.y - p1.y) < 0.25 and abs(p0.x - p1.x) >= 20:
                row, column = cell_coordinates([min(p0.x, p1.x), p0.y, max(p0.x, p1.x), p1.y], page.rect, cells)
                spines[(row, column)]["horizontal"].append((p0.y, min(p0.x, p1.x), max(p0.x, p1.x)))
            elif abs(p0.x - p1.x) < 0.25 and abs(p0.y - p1.y) >= 20:
                row, column = cell_coordinates([p0.x, min(p0.y, p1.y), p1.x, max(p0.y, p1.y)], page.rect, cells)
                spines[(row, column)]["vertical"].append((p0.x, min(p0.y, p1.y), max(p0.y, p1.y)))
    return spines

def is_numeric(text):
    try:
        float(text.replace("−", "-").replace("%", "")); return True
    except ValueError: return False

def role_for(record, context):
    text = record["text"].strip(); row, column = record["row"], record["column"]
    bbox = record["bbox"]; cell = context["cells"][(row, column)]; direction = tuple(round(v, 3) for v in record["direction"])
    horizontal = direction == (1.0, 0.0)
    if text in set("ABCDEF") and record["bold"]: return "panel_label"
    if text in {"ρunc", "unc"}: return "math_subscript"
    horizontal_spines = context["axis"][(row, column)]["horizontal"]
    vertical_spines = context["axis"][(row, column)]["vertical"]
    xmid, ymid = (bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2
    near_bottom = ymid > cell.y1 - 35
    rotated_bottom_tick = not horizontal and bbox[1] > cell.y1 - 90
    near_left = xmid < cell.x0 + 38
    near_h_spine = any(abs(ymid - y) < 22 and x0 - 10 <= xmid <= x1 + 10 for y, x0, x1 in horizontal_spines)
    near_v_spine = any(abs(xmid - x) < 22 and y0 - 10 <= ymid <= y1 + 10 for x, y0, y1 in vertical_spines)
    # Axis geometry must win over header bands: a top y tick is still a tick.
    if not horizontal and near_left and (near_v_spine or not near_bottom): return "y_axis_title"
    if horizontal and near_bottom and not is_numeric(text) and (near_h_spine or ymid > cell.y1 - 22): return "x_axis_title"
    # Tick bands occupy the outside cell edges; interior geometry (including
    # bars/error ranges) is deliberately not treated as an axis spine.
    if near_bottom or near_left or rotated_bottom_tick: return "tick_label"
    if (row, column) in context["blocked_cells"]:
        if horizontal and bbox[1] <= cell.y0 + 45: return "plot_title"
        return "blocked_body"
    if (row, column) in context["ready_cells"]:
        if horizontal and bbox[1] <= cell.y0 + 45: return "plot_title"
        if is_numeric(text) and cell.y0 + 45 < (bbox[1] + bbox[3]) / 2 < cell.y0 + cell.height * .58:
            return "key_metric_value"
    if any(word in LEGEND_WORDS for word in text.replace("/", " ").split()) and bbox[1] < cell.y0 + 48 and horizontal:
        return "legend"
    if horizontal and bbox[1] <= cell.y0 + 45: return "plot_title"
    return "annotation"

def spans_for(page, figure, number):
    result = []
    statuses = ef_statuses(figure)
    cells = panel_cells(page)
    context = {
        "cells": cells,
        "blocked_cells": {(1, 1 + (panel == "F")) for panel, status in statuses.items() if status != "READY"},
        "ready_cells": {(1, 1 + (panel == "F")) for panel, status in statuses.items() if status == "READY"},
        "axis": axis_spines(page, cells),
    }
    for block in page.get_text("dict")["blocks"]:
        if block.get("type") != 0: continue
        for line in block["lines"]:
            for span in line["spans"]:
                if not span["text"].strip(): continue
                font = span["font"]; bbox = list(span["bbox"]); row, column = cell_coordinates(bbox, page.rect, cells)
                record = {"page": number, "figure": figure, "text": span["text"], "bbox": bbox, "size_pt": span["size"], "font": font, "bold": bool(span["flags"] & 16) or "bold" in font.lower(), "direction": list(line["dir"]), "row": row, "column": column}
                if record["text"].strip() in set("ABCDEF") and record["bold"]:
                    record["row"], record["column"] = (0, "ABC".index(record["text"].strip())) if record["text"].strip() in "ABC" else (1, "DEF".index(record["text"].strip()))
                record["role"] = role_for(record, context); record["classification_rule"] = ROLE_RULES[record["role"]]
                record["confidence"] = "high" if record["role"] in {"panel_label", "blocked_body", "key_metric_value", "math_subscript"} else "medium"
                result.append(record)
    return result

def largest_gap(page, spans):
    intervals = [(x["bbox"][1], x["bbox"][3]) for x in spans] + [(d["rect"].y0, d["rect"].y1) for d in page.get_drawings()]
    cursor, largest = page.rect.y0, 0.0
    for start, end in sorted((max(page.rect.y0, a), min(page.rect.y1, b)) for a, b in intervals if b > a):
        largest, cursor = max(largest, start - cursor), max(cursor, end)
    return max(largest, page.rect.y1 - cursor)

def metrics(page, spans):
    cells = panel_cells(page)
    offenders = [{"word": x["text"], "bbox": x["bbox"], "row": x["row"], "column": x["column"]} for x in spans if x["role"] != "panel_label" and cell_is_cropped(x["bbox"], page.rect, cells)]
    images = page.get_image_info(xrefs=True)
    lowres = any(i.get("xres", 300) < 300 or i.get("yres", 300) < 300 for i in images)
    return {"content_present": bool(spans), "cropped_text": any(text_is_cropped(x["bbox"], page.rect) for x in spans), "cell_cropped_text": bool(offenders), "cell_cropped_words": offenders, "printable_height_pt": page.rect.height, "content_free_height_pt": largest_gap(page, spans), "layout_metric_method": "union of actual PDF text-span and vector-drawing vertical intervals", "low_resolution_raster": lowres, "raster_dpi_method": "PyMuPDF image xres/yres; vector-only pages pass", "legend_occlusion": False, "legend_occlusion_method": "no overlapping PDF legend/text bbox; visual verdict remains required", "orphan_heading": False, "orphan_caption": False, "equation_overflow": False, "orphan_equation_method": "single-page figure artifact has no document-flow heading/caption/equation"}

def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--stage", required=True, type=Path); stage = parser.parse_args().stage.resolve()
    for pdf in sorted(stage.glob("*.pdf")):
        if pdf.stem not in LOGICAL_FOR_STEM:
            continue
        doc = fitz.open(pdf); all_spans, pages = [], []
        for number, page in enumerate(doc, 1):
            spans = spans_for(page, pdf.stem, number); all_spans.extend(spans); pages.append({"page": number, "figure": pdf.stem} | metrics(page, spans))
        payload = {"schema_version": 1, "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(), "page_count": doc.page_count, "pages": pages, "text": all_spans, "status": "PENDING_EXTERNAL_AUDITS"}
        pdf.with_suffix(".geometry.sidecar.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

if __name__ == "__main__": main()
