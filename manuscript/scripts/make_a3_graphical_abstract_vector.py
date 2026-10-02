#!/usr/bin/env python3
"""Render the current A3 graphical abstract as native vector artwork.

The figure is an evidence map, not a new analysis.  It reads only the saved
current score, calibration, paired-DMS and archived-binding records.  AUROC,
Brier and Spearman retain their native units and are deliberately not merged
into one visual scale.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="a3-ga-mpl-"))

import fitz
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, ft2font
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


MANUSCRIPT = Path(__file__).resolve().parents[1]
CAPSULE = MANUSCRIPT.parent
STEM = "A3_graphical_abstract"
WIDTH_PT, HEIGHT_PT = 468.0, 309.6  # 6.5 x 4.3 in
MIN_TEXT_PT = 7.5
INK, MUTED, RULE = "#17232E", "#495866", "#CCD5DD"
BLUE, GREEN, AMBER = "#0072B2", "#007E63", "#9B4A10"
PALE = {"score": "#F5F8FA", "audit": "#F5F8FA",
        "clinical": "#F1F7FA", "dms": "#F0F8F5", "binding": "#FCF6F0"}

# These are the current source identities at handoff.  A changed input is a
# hard error: the GA must never silently drift away from the bound evidence.
BOUND_INPUTS = {
    "analysis/revision_20260923/row-scores/provenance.json":
        "f23da8217a42424851dd19cb457b0283db4340647b32ec818571aa7e887bb687",
    "analysis/revision_20260923/row-scores/scores.csv":
        "0fac1a086b5b4797b336484f8f0a6ee6eb00334b83607c353b4ea4cd5788eccd",
    "analysis/revision_20260923/clinical-unified/summary.json":
        "4d3029e9cf0a98f3dcd1672d7a4d0ce26d354c3e97d915542ebcb3403acd38ae",
    "analysis/revision_20260923/clinical-unified/oof_predictions.csv":
        "e658cff982102527a8d84c799390ccfdc4774d32c7ca668bbc2985997027a800",
    "manuscript/figures/plotdata/clinical_macro_20260923.csv":
        "88b40c659a5be6beb6f7531e19ac37db2af9f2aa11f74ac1bc1fd147b7f63854",
    "manuscript/figures/plotdata/clinical_metrics_20260923.csv":
        "40e9a0aaefd18fe87c68285f880fcb40f46c864c5ba38b5ace8c1dbe88671fbd",
    "manuscript/figures/plotdata/F1_scale_paired_650M_3B.json":
        "1eef587cd6a1ecbe370c8426ed5d15fc1cd4f25de8b9d40cf88cf0456f248df3",
    "manuscript/figures/plotdata/F1_scale_paired_650M_3B.csv":
        "6496354cf21a4fc87954ace959d754bcc0f470f3dd6c7dbcde3177ab1a2c3d8c",
    "manuscript/figures/plotdata/FS3_antibody_robustness.csv":
        "96205842095ca818a7bc8744c19acfe6a17ae448f4ab135c5b99d1ba740971fa",
    "manuscript/figures/plotdata/F8_antibody_n4.csv": None,
    "manuscript/figures/plotdata/F5_developability.csv":
        "17c2f64b133c446ca316169d727f2e114308013fc515f5a32743c57b1b8ab878",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def close(a: float, b: float, tolerance: float = 1e-13) -> None:
    if not math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tolerance):
        raise ValueError(f"Unreconciled bound values: {a} != {b}")


def load_inputs(capsule: Path = CAPSULE) -> dict:
    source_hashes: dict[str, str] = {}
    for relative, expected in BOUND_INPUTS.items():
        path = capsule / relative
        if not path.is_file():
            raise ValueError(f"Missing bound input: {relative}")
        actual = sha256(path)
        if expected is not None and actual != expected:
            raise ValueError(f"GA bound input changed: {relative}")
        source_hashes[relative] = actual

    summary = json.loads(
        (capsule / "analysis/revision_20260923/clinical-unified/summary.json").read_text()
    )
    provenance = json.loads(
        (capsule / "analysis/revision_20260923/row-scores/provenance.json").read_text()
    )
    scores = read_csv(capsule / "analysis/revision_20260923/row-scores/scores.csv")
    oof = read_csv(capsule / "analysis/revision_20260923/clinical-unified/oof_predictions.csv")
    score_ids = {row["example_id"] for row in scores}
    oof_ids = {row["example_id"] for row in oof}
    if not (
        len(scores) == len(score_ids) == len(oof) == 4677
        and score_ids == oof_ids
        and len({row["gene"] for row in scores}) == 20
        and summary["n_rows"] == provenance["n_rows"] == 4677
        and summary["n_genes"] == 20
        and summary["folds"] == 5
        and provenance["scorer"] == "N6c fp32 WT-marginal"
        and provenance["scores_sha256"] == source_hashes["analysis/revision_20260923/row-scores/scores.csv"]
    ):
        raise ValueError("Clinical same-score / identity / panel contract failed")

    macro = read_csv(capsule / "manuscript/figures/plotdata/clinical_macro_20260923.csv")[0]
    close(summary["macro_auroc"], float(macro["value"]))
    close(summary["gene_mean_ci"][0], float(macro["lo"]))
    close(summary["gene_mean_ci"][1], float(macro["hi"]))

    metrics = {
        row["map"]: row
        for row in read_csv(capsule / "manuscript/figures/plotdata/clinical_metrics_20260923.csv")
    }
    for name in ("platt", "isotonic", "gene_prevalence"):
        close(summary["maps"][name]["point"]["brier"], float(metrics[name]["brier"]))

    dms = json.loads(
        (capsule / "manuscript/figures/plotdata/F1_scale_paired_650M_3B.json").read_text()
    )
    pairs = read_csv(capsule / "manuscript/figures/plotdata/F1_scale_paired_650M_3B.csv")
    contrast = math.fsum(float(row["delta_3b_minus_650m"]) for row in pairs) / len(pairs)
    if len(pairs) != dms["paired_n_proteins"] or len(pairs) != 69:
        raise ValueError("DMS paired protein count failed")
    close(contrast, dms["paired_bootstrap"]["estimate"], 2e-12)
    if not dms["paired_bootstrap"]["ci_low"] < 0 < dms["paired_bootstrap"]["ci_high"]:
        raise ValueError("Saved DMS contrast no longer spans zero")

    binding_rows = read_csv(capsule / "manuscript/figures/plotdata/FS3_antibody_robustness.csv")
    binding = {
        row["subset"]: row
        for row in binding_rows
        if row["method"] == "ESM-2 650M (WT)"
    }
    alternate = next(
        row
        for row in read_csv(capsule / "manuscript/figures/plotdata/F8_antibody_n4.csv")
        if row["model"] == "ESM-2 650M"
    )
    close(float(binding["antibody"]["rho"]), float(alternate["rho"]), 0.00005)
    if float(binding["antibody"]["ci_lo"]) == float(alternate["ci_lo"]):
        raise ValueError("Distinct archived intervals unexpectedly collapsed")

    return {
        "sources": source_hashes,
        "clinical": summary,
        "score_provenance": provenance,
        "dms": dms,
        "binding": binding,
        "alternate_archive": alternate,
    }


def required_fonts() -> dict[str, Path]:
    found: dict[str, Path] = {}
    for key, pattern, is_bold in (
        ("regular", "Arial:style=Regular", False),
        ("bold", "Arial:style=Bold", True),
    ):
        match = subprocess.run(
            ["fc-match", "-f", "%{file}", pattern],
            text=True,
            capture_output=True,
            check=True,
        ).stdout.strip()
        path = Path(match)
        face = ft2font.FT2Font(str(path))
        if face.family_name != "Arial" or ("Bold" in face.style_name) != is_bold:
            raise ValueError(f"Arial {key} is unavailable; font fallback refused")
        found[key] = path
    return found


class Drawing:
    def __init__(self, fonts: dict[str, Path]):
        self.fonts = fonts
        self.figure = plt.figure(
            figsize=(WIDTH_PT / 72.0, HEIGHT_PT / 72.0), dpi=300, facecolor="white"
        )
        self.axis = self.figure.add_axes((0, 0, 1, 1))
        self.axis.set(xlim=(0, WIDTH_PT), ylim=(HEIGHT_PT, 0))
        self.axis.axis("off")
        self.texts: list[dict] = []

    def text(
        self,
        x: float,
        y: float,
        content: str,
        size: float = 8.0,
        color: str = INK,
        role: str = "body",
        align: str = "left",
        region: tuple[float, float, float, float] | None = None,
    ):
        if role == "panel_label" and content not in set("ABCDE"):
            raise ValueError("Only panel labels may use bold")
        if size < MIN_TEXT_PT:
            raise ValueError(f"GA authored text below floor: {size}")
        key = "bold" if role == "panel_label" else "regular"
        face = ft2font.FT2Font(str(self.fonts[key]))
        missing = sorted(
            {
                ord(char)
                for char in content
                if char not in "\n\t" and ord(char) not in face.get_charmap()
            }
        )
        if missing:
            raise ValueError(f"Arial missing glyphs: {missing}")
        artist = self.axis.text(
            x,
            y,
            content,
            fontsize=size,
            color=color,
            # A font file selects the PDF/outline face, but editable SVG CSS
            # follows FontProperties.weight and otherwise silently stays normal.
            fontproperties=font_manager.FontProperties(
                fname=str(self.fonts[key]), weight="bold" if key == "bold" else "normal"
            ),
            ha=align,
            va="top",
            linespacing=1.16,
            clip_on=False,
        )
        self.texts.append(
            {
                "artist": artist,
                "text": content,
                "size_pt": size,
                "role": role,
                "font_face": key,
                "region": region,
            }
        )
        return artist

    def line(self, x0, y0, x1, y1, color=RULE, width=0.7, dashed=False):
        self.axis.plot(
            (x0, x1),
            (y0, y1),
            color=color,
            lw=width,
            linestyle=(0, (2, 2)) if dashed else "-",
            solid_capstyle="round",
            clip_on=False,
        )

    def circle(self, x, y, radius, color, fill=None, width=0.8):
        self.axis.add_patch(
            Circle((x, y), radius, facecolor=fill or "white", edgecolor=color, linewidth=width)
        )

    def rectangle(self, x, y, width, height, fill="white", color=RULE, lw=0.7):
        self.axis.add_patch(
            Rectangle((x, y), width, height, facecolor=fill, edgecolor=color, linewidth=lw)
        )

    def box(self, x, y, width, height, fill):
        self.axis.add_patch(
            FancyBboxPatch(
                (x, y),
                width,
                height,
                boxstyle="round,pad=0,rounding_size=2",
                facecolor=fill,
                edgecolor=RULE,
                linewidth=0.65,
            )
        )

    def arrow(self, x0, y0, x1, y1, color=MUTED):
        self.axis.add_patch(
            FancyArrowPatch(
                (x0, y0),
                (x1, y1),
                arrowstyle="-|>",
                mutation_scale=5.8,
                linewidth=0.8,
                color=color,
                shrinkA=0,
                shrinkB=0,
            )
        )

    def title(self, x, y, letter, title):
        self.text(x, y, letter, 10.8, role="panel_label")
        self.text(x + 17, y + 0.6, title, 8.8, role="title")

    def text_bounds(self) -> list[dict]:
        self.figure.canvas.draw()
        renderer = self.figure.canvas.get_renderer()
        records: list[dict] = []
        bounds: list[tuple[str, tuple[float, float, float, float]]] = []
        transform = self.axis.transData.inverted()
        for row in self.texts:
            display = row["artist"].get_window_extent(renderer)
            corners = transform.transform(
                [(display.x0, display.y0), (display.x0, display.y1),
                 (display.x1, display.y0), (display.x1, display.y1)]
            )
            bbox = (
                float(corners[:, 0].min()),
                float(corners[:, 1].min()),
                float(corners[:, 0].max()),
                float(corners[:, 1].max()),
            )
            if not (2 <= bbox[0] < bbox[2] <= WIDTH_PT - 2 and
                    2 <= bbox[1] < bbox[3] <= HEIGHT_PT - 2):
                raise ValueError(f"Text exceeds page safe area: {row['text']}: {bbox}")
            region = row["region"]
            if region is not None:
                x0, y0, x1, y1 = region
                if not (x0 <= bbox[0] and y0 <= bbox[1] and bbox[2] <= x1 and bbox[3] <= y1):
                    raise ValueError(f"Text exceeds its panel: {row['text']}: {bbox}")
            for previous, other in bounds:
                intersect_w = min(bbox[2], other[2]) - max(bbox[0], other[0])
                intersect_h = min(bbox[3], other[3]) - max(bbox[1], other[1])
                if intersect_w > 0.25 and intersect_h > 0.25:
                    raise ValueError(f"Text collision: {row['text']} / {previous}")
            bounds.append((row["text"], bbox))
            records.append({key: value for key, value in row.items() if key != "artist"} |
                           {"bbox_pt": bbox})
        return records


def draw(evidence: dict, fonts: dict[str, Path]) -> Drawing:
    scene = Drawing(fonts)
    scene.text(12, 9, "Task-specific reliability of zero-shot protein scores", 10.5, role="title")
    scene.text(12, 24, "Sequence likelihood must be interpreted through endpoint-specific evidence.", 7.5, color=MUTED)

    # A/B: score object and the audit that determines which downstream claim
    # is actually supported.
    scene.title(12, 42, "A", "Sequence likelihood score")
    scene.title(236, 42, "B", "Evidence audit")
    scene.box(12, 57, 202, 78, PALE["score"])
    scene.box(236, 57, 220, 78, PALE["audit"])

    for y, label, fill in ((79, "WT", BLUE), (102, "mut", GREEN)):
        scene.text(23, y - 5, label, 7.5, color=MUTED)
        points = [(50 + i * 6.5, y + (i % 2) * 2.5) for i in range(5)]
        for first, second in zip(points, points[1:]):
            scene.line(*first, *second, color=BLUE, width=1.05)
        for index, (x, yy) in enumerate(points):
            scene.circle(x, yy, 2.25, fill if index == 2 else BLUE,
                         fill if index == 2 else "white")
    scene.text(95, 73, "L = log p(mut | context)", 7.8)
    scene.text(95, 87, "    - log p(wt | context)", 7.8)
    scene.text(95, 104, "Context-specific scoring", 7.5, color=MUTED)
    scene.text(23, 123, "LLR is a score, not a probability.", 7.5, color=MUTED)
    scene.arrow(219, 96, 230, 96)

    scene.rectangle(248, 72, 16, 43, fill="white", color=MUTED)
    for y in (80, 93, 106):
        scene.line(252, y, 260, y, color=MUTED)
        scene.circle(272, y, 1.5, MUTED, "white", width=0.65)
    scene.text(279, 75, "Endpoint & native metric", 7.8)
    scene.text(279, 91, "Score version & row identity", 7.8)
    scene.text(279, 107, "Control applicability", 7.8)
    scene.text(248, 122, "Single-parent: germline control N/A", 7.5, color=MUTED)

    # Common audit spine; the three endpoint panels below are deliberately not
    # plotted on a shared numerical axis.
    scene.line(346, 135, 346, 145, color=MUTED)
    scene.line(82, 145, 385, 145, color=MUTED)
    for x in (82, 234, 385):
        scene.arrow(x, 145, x, 152)

    panels = {
        "C": (12, 143, 140, "Clinical", "clinical"),
        "D": (163, 142, 140, "Paired DMS", "dms"),
        "E": (314, 143, 142, "Binding / antibody", "binding"),
    }
    panel_regions = {}
    for label, (x, width, _y, title, kind) in panels.items():
        scene.title(x, 155, label, title)
        scene.box(x, 171, width, 121, PALE[kind])
        panel_regions[label] = (x + 8, 176, x + width - 8, 288)

    clinical = evidence["clinical"]
    c = panel_regions["C"]
    for i in range(5):
        scene.rectangle(23 + (i % 3) * 6, 181 + (i // 3) * 7, 4.5, 5,
                        fill=BLUE if i < 4 else GREEN,
                        color=BLUE if i < 4 else GREEN, lw=0.4)
    scene.text(52, 179, f"{clinical['n_rows']:,} variants", 7.8, color=BLUE, region=c)
    scene.text(52, 192, f"{clinical['n_genes']} genes; same LLR", 7.5, region=c)
    scene.line(23, 208, 144, 208)
    scene.text(23, 214, f"Macro-AUROC  {clinical['macro_auroc']:.4f}", 7.9, color=BLUE, region=c)
    scene.text(23, 227, "5-fold within-panel calibration", 7.5, color=MUTED, region=c)
    scene.text(23, 242, "Brier", 7.7, color=MUTED, region=c)
    scene.text(144, 242, "lower is better", 7.5, color=MUTED, align="right", region=c)
    for y, name, label in ((255, "platt", "Platt"), (267, "isotonic", "Isotonic"),
                           (279, "gene_prevalence", "Gene prior")):
        scene.text(23, y, label, 7.5, region=c)
        scene.text(144, y, f"{clinical['maps'][name]['point']['brier']:.4f}",
                   7.5, color=BLUE, align="right", region=c)

    d = panel_regions["D"]
    for y in (181, 194):
        for i in range(3):
            if i:
                scene.line(175 + (i - 1) * 6, y, 175 + i * 6, y, color=GREEN)
            scene.circle(175 + i * 6, y, 1.5, GREEN, "white", width=0.7)
    scene.line(181, 184, 181, 191, color=GREEN)
    dms = evidence["dms"]
    effect = dms["paired_bootstrap"]
    scene.text(202, 179, f"{dms['paired_n_proteins']} paired proteins", 7.7, color=GREEN, region=d)
    scene.text(202, 192, "3B minus 650M", 7.5, region=d)
    scene.line(174, 208, 294, 208)
    scene.text(174, 214, f"AUROC difference  {effect['estimate']:.5f}", 7.5, color=GREEN, region=d)
    scale = lambda value: 234 + 60 * value / 0.02
    scene.line(scale(0), 228, scale(0), 241, color=MUTED, width=0.6, dashed=True)
    scene.line(scale(effect["ci_low"]), 234, scale(effect["ci_high"]), 234, color=GREEN, width=1.05)
    for endpoint in ("ci_low", "ci_high"):
        scene.line(scale(effect[endpoint]), 231, scale(effect[endpoint]), 237, color=GREEN, width=0.8)
    scene.circle(scale(effect["estimate"]), 234, 1.9, GREEN, GREEN)
    scene.text(174, 246, f"95% CI [{effect['ci_low']:.4f}, {effect['ci_high']:.4f}]", 7.5, region=d)
    scene.text(174, 263, "No gain detected", 7.8, color=GREEN, region=d)
    scene.text(174, 278, "Not an equivalence test", 7.5, color=MUTED, region=d)

    e = panel_regions["E"]
    for x0, y0, x1, y1 in ((324, 181, 332, 189), (340, 181, 332, 189),
                            (332, 189, 332, 199), (323, 185, 329, 191),
                            (341, 185, 335, 191)):
        scene.line(x0, y0, x1, y1, color=AMBER, width=1.2)
    scene.text(352, 179, "Archived binding", 7.7, color=AMBER, region=e)
    scene.text(352, 192, "Spearman, not AUROC", 7.5, region=e)
    scene.line(324, 208, 445, 208)
    scene.text(324, 214, "ESM-2 650M WT-marginal", 7.5, color=MUTED, region=e)
    for y, name, label in ((228, "general", "General (non-Ab)"), (241, "antibody", "Ab-antigen")):
        scene.text(324, y, label, 7.5, region=e)
        scene.text(445, y, f"{float(evidence['binding'][name]['rho']):.3f}",
                   7.5, color=AMBER, align="right", region=e)
    scene.line(324, 255, 445, 255)
    scene.text(324, 261, "Row/scorer/resampling unresolved.", 7.5, color=AMBER, region=e)
    scene.text(324, 270, "Same SKEMPI set; archived points", 7.5, color=MUTED, region=e)
    scene.text(324, 280, "only; no model claim.", 7.5, color=MUTED, region=e)

    scene.text(12, 299, "Native AUROC, Brier and Spearman units; no common reliability scale or external clinical validation.",
               7.5, color=MUTED)
    return scene


def validate_pdf(path: Path) -> dict:
    document = fitz.open(path)
    if len(document) != 1:
        raise ValueError("GA PDF must contain one page")
    page = document[0]
    if abs(page.rect.width - WIDTH_PT) > 0.02 or abs(page.rect.height - HEIGHT_PT) > 0.02:
        raise ValueError(f"GA page geometry changed: {page.rect}")
    bold = []
    sizes: list[float] = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line.get("spans", []):
                text = span["text"].strip()
                if not text:
                    continue
                if "Arial" not in span["font"]:
                    raise ValueError(f"Non-Arial PDF face: {span}")
                if span["size"] < MIN_TEXT_PT:
                    raise ValueError(f"GA glyph below floor: {span}")
                sizes.append(float(span["size"]))
                if "Bold" in span["font"]:
                    bold.append(text)
    if bold != list("ABCDE"):
        raise ValueError(f"Only panel labels may be bold, found: {bold}")
    report = subprocess.run(["pdffonts", str(path)], text=True, capture_output=True, check=True).stdout
    rows = [line for line in report.splitlines()[2:] if line.strip() and not line.startswith("-")]
    if not rows:
        raise ValueError("pdffonts returned no fonts")
    for row in rows:
        fields = row.split()
        if ("Type 3" in row or "Arial" not in row or len(fields) < 9 or
                fields[4:7] != ["yes", "yes", "yes"]):
            raise ValueError(f"PDF font embedding check failed: {row}")
    return {
        "page_count": 1,
        "page_size_pt": [WIDTH_PT, HEIGHT_PT],
        "embedded_fonts": True,
        "non_label_bold": 0,
        "minimum_text_pt": min(sizes),
        "font_report": report,
    }


def svg_profile(path: Path) -> dict:
    root = ET.fromstring(path.read_bytes())
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    text_nodes = root.findall(".//svg:text", namespace)
    path_nodes = root.findall(".//svg:path", namespace)
    image_nodes = root.findall(".//svg:image", namespace)
    if image_nodes:
        raise ValueError("GA SVG unexpectedly contains bitmap content")
    bold: list[str] = []
    sizes: list[float] = []
    for node in text_nodes:
        content = "".join(node.itertext()).strip()
        style = {
            key.strip(): value.strip()
            for item in node.get("style", "").split(";") if item.strip()
            for key, value in [item.split(":", 1)]
        }
        if "Arial" not in style.get("font-family", ""):
            raise ValueError(f"Non-Arial SVG face: {content}")
        size = float(style["font-size"].removesuffix("px"))
        if size < MIN_TEXT_PT:
            raise ValueError(f"SVG text below floor: {content}: {size}")
        sizes.append(size)
        weight = style.get("font-weight", "normal")
        if content in set("ABCDE"):
            if weight != "700":
                raise ValueError(f"SVG panel label {content} requires font-weight: 700")
            bold.append(content)
        elif weight not in {"normal", "400"}:
            raise ValueError(f"Only panel labels may be bold in SVG: {content}")
    if text_nodes and bold != list("ABCDE"):
        raise ValueError(f"SVG must contain exactly five bold panel labels: {bold}")
    return {
        "text_elements": len(text_nodes),
        "path_elements": len(path_nodes),
        "image_elements": len(image_nodes),
        "viewBox": root.get("viewBox"),
        "width": root.get("width"),
        "height": root.get("height"),
        "bold_texts": bold,
        "non_label_bold": 0,
        "minimum_text_pt": min(sizes) if sizes else None,
    }


def export(output: Path) -> dict:
    output = output.resolve()
    if not output.is_dir():
        raise ValueError("Explicit, existing GA output directory required")
    evidence = load_inputs()
    fonts = required_fonts()
    matplotlib.rcParams.update({
        "font.family": "Arial",
        "font.weight": "normal",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "pdf.use14corefonts": False,
        "pdf.compression": 9,
        "svg.hashsalt": "A3-GA-20260930-v1",
        "text.usetex": False,
    })
    scene = draw(evidence, fonts)
    typography = scene.text_bounds()
    fixed = dt.datetime(1970, 1, 1, tzinfo=dt.timezone.utc)
    suffixes = ("pdf", "png", "svg", "display.svg")
    with tempfile.TemporaryDirectory(prefix="a3-ga-export-") as directory:
        stage = Path(directory)
        scene.figure.savefig(
            stage / f"{STEM}.pdf",
            metadata={"Title": "A3 graphical abstract: endpoint-specific evidence audit",
                      "Author": "", "CreationDate": fixed, "ModDate": fixed,
                      "Subject": "Current author-review figure; not scientific clearance"},
        )
        scene.figure.savefig(
            stage / f"{STEM}.png", dpi=300,
            metadata={"Software": "A3 GA native vector renderer"},
        )
        for suffix, fonttype in (("svg", "none"), ("display.svg", "path")):
            with matplotlib.rc_context({"svg.fonttype": fonttype}):
                scene.figure.savefig(
                    stage / f"{STEM}.{suffix}",
                    metadata={"Title": "A3 graphical abstract: endpoint-specific evidence audit",
                              "Date": "1970-01-01",
                              "Description": "Conceptual icons and saved observations; native endpoint units."},
                )
        normalizer = MANUSCRIPT / "scripts/normalize_a3_pdf_metadata.py"
        subprocess.run(["python3", str(normalizer), str(stage / f"{STEM}.pdf")], check=True)
        pdf_checks = validate_pdf(stage / f"{STEM}.pdf")
        profiles = {
            suffix: svg_profile(stage / f"{STEM}.{suffix}")
            for suffix in ("svg", "display.svg")
        }
        if profiles["svg"]["text_elements"] != len(typography):
            raise ValueError("Editable SVG text count does not match authored typography")
        if profiles["display.svg"]["text_elements"] or profiles["display.svg"]["path_elements"] < 100:
            raise ValueError("Outline SVG text/path contract failed")
        for suffix in suffixes:
            # The temporary directory and the capsule may be different
            # filesystems in Codex. Copying preserves the bytes while avoiding
            # a cross-device rename failure.
            shutil.copyfile(stage / f"{STEM}.{suffix}", output / f"{STEM}.{suffix}")
    plt.close(scene.figure)

    sources = dict(evidence["sources"])
    sources.update({
        "manuscript/text/methods.tex": sha256(MANUSCRIPT / "text/methods.tex"),
        "manuscript/scripts/R/make_a3_graphical_abstract.R": sha256(
            MANUSCRIPT / "scripts/R/make_a3_graphical_abstract.R"
        ),
        "manuscript/scripts/make_a3_graphical_abstract_vector.py": sha256(Path(__file__)),
    })
    payload = {
        "schema": "a3-graphical-abstract-native-vector-v1",
        "status": "CURRENT_AUTHOR_REVIEW_GA_NOT_SCIENTIFIC_CLEARANCE",
        "scientific_ready_attested": False,
        "source_date_policy": "deterministic metadata; no wall-clock timestamp",
        "dimensions": {"width_in": 6.5, "height_in": 4.3,
                       "page_size_pt": [WIDTH_PT, HEIGHT_PT], "png_dpi": 300},
        "command": "A3_FIGURE_OUTPUT=<existing-directory> Rscript scripts/R/make_a3_graphical_abstract.R",
        "sources_sha256": sources,
        "output_sha256": {
            f"{STEM}.{suffix}": sha256(output / f"{STEM}.{suffix}") for suffix in suffixes
        },
        "fonts": {
            key: {"family": "Arial", "face": key, "local_path": str(path),
                  "sha256": sha256(path), "redistributed": False}
            for key, path in fonts.items()
        },
        "typography": typography,
        "svg_profiles": profiles,
        "automatic_geometry_checks": pdf_checks | {
            "text_collision_pairs": 0,
            "out_of_panel_text": 0,
            "non_label_bold": 0,
            "minimum_authored_text_pt": min(row["size_pt"] for row in typography),
            "manual_visual_review": "Required separately; automatic bounds are not a visual verdict.",
        },
        "numerical_bindings": {
            "clinical": {
                "protocol": evidence["clinical"]["protocol"],
                "scorer": evidence["score_provenance"]["scorer"],
                "n_variants": 4677,
                "n_genes": 20,
                "macro_auroc": evidence["clinical"]["macro_auroc"],
                "macro_gene_bootstrap_ci": evidence["clinical"]["gene_mean_ci"],
                "brier": {
                    name: evidence["clinical"]["maps"][name]["point"]["brier"]
                    for name in ("platt", "isotonic", "gene_prevalence")
                },
                "interpretation": "Same recovered score; fitting-fold-only five-fold calibration; fixed panel, not external validation.",
            },
            "paired_dms": evidence["dms"],
            "archived_binding": {
                "source": "manuscript/figures/plotdata/FS3_antibody_robustness.csv",
                "general": evidence["binding"]["general"],
                "antibody": evidence["binding"]["antibody"],
                "alternate_archive": evidence["alternate_archive"],
                "subset_semantics": {
                    "cohort": "same SKEMPI set",
                    "partition": "Table 2 partition A",
                    "relationship": "complementary_non_overlapping",
                    "general": "General (non-Ab)",
                    "antibody": "Ab-antigen",
                    "n_total": int(evidence["binding"]["general"]["n"]) + int(
                        evidence["binding"]["antibody"]["n"]
                    ),
                    "independent_experiments": False,
                },
                "interpretation": "Same SKEMPI set; Table 2 partition A contains complementary general (non-Ab) and antibody-antigen subsets; archived points only; row/scorer/resampling identity unresolved. No independent complex-level CI or model-comparison inference.",
            },
        },
        "scope": [
            "Icons are schematic, not an experimental observation or predicted molecular structure.",
            "AUROC, Brier and Spearman are not merged, null-shifted or transformed into a common reliability score.",
            "DMS absence of detected gain is not equivalence or a universal scaling plateau.",
            "Single-parent germline control is inapplicable, not a negative discovery.",
            "Archived 650M binding and antibody scorer identities remain unreconciled; alternate intervals are not pooled.",
            "No model inference, resampling, manuscript compilation, other-figure edit or method-contract rebinding.",
        ],
    }
    (output / f"{STEM}.provenance.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = export(args.output)
    print(json.dumps({key: result[key] for key in ("status", "dimensions", "output_sha256", "svg_profiles")}, indent=2))
