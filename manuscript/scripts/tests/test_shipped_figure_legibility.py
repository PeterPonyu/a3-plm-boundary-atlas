"""Legibility geometry measured on the SHIPPED figures, not on a staged copy.

Why this file exists (2026-09-09)
---------------------------------
A3 already carries font-role, panel-label-clearance and layout contract tests,
and all of them were green while Fig 2 panel A had twenty gene names struck
through each other. Two reasons, both real:

  1. `test_figure_font_roles.py`, `test_panel_label_clearance.py` and
     `test_pdf_layout_contract.py` build their sidecars and fixture PDFs
     in-test. They audit the checker's RULES, never a shipped figure's glyphs.
     `test_authority_font_audit_is_expected_red_without_role_metadata` even
     asserts the real authority PDF comes back UNSUPPORTED. Tightening a
     threshold in those files cannot catch a shipped-figure defect.

  2. `test_a3_six_panel_compositor.py` does read real rendered PDFs -- but from
     `.superpowers/sdd/.../task-6-stage/`, which was last written 2026-08-24 and
     whose F1/F2/F3/F6 bytes all differ from the current figures. Those tests
     were green against artifacts sixteen days older than what ships.

So the gap was not a threshold. It was that nothing measured the artifact the
reader actually receives. This file measures `figures/*.pdf` directly.

What it guarantees: no axis tick label in any shipped figure vertically
overlaps its neighbour, no glyph prints below 7.5 pt after the 6.5-inch / 0.70
text-height constraint, embedded Arial with only panel letters bold, aligned
tag rows/columns, net tag-to-text/axis clearance and in-page text bounds.
It would have failed on the pre-fix Fig 2 panel A (4.07 pt pitch under a 6.86 pt
glyph height, 19 of 19 adjacent pairs overlapping).
"""
from __future__ import annotations

import re
import math
import unittest
from collections import defaultdict
from pathlib import Path

import fitz

ROOT = Path(__file__).resolve().parents[2]
FIGURES = ROOT / "figures"

# Current main/SI figures, including the real-structure binding face. The candidate
# F7_20gene_ece pair is deliberately unpromoted and is not audited here.
SHIPPED = (
    "F1_scale_plateau", "F2_cross_model", "F3_clinvar_forest",
    "F4_clinical_calibration", "F5_stability_ddg", "F6_protein_uq",
    "F7_binding_boundary", "F8_developability_map", "F9_antibody_lm",
    "F10_reliability_predictor", "F11_boundary_atlas",
    "FS1_proteingym_atlas", "FS2_gnomad_control", "FS3_antibody_robustness",
    "FS4_abstention_tiers", "FS5_vep_comparison",
)

# Measure the final includegraphics constraint, not just authored glyph size.
# Letter paper with one-inch margins: 6.5-inch width and 9-inch text height.
PRINT_WIDTH_BP = 6.5 * 72
PRINT_MAX_HEIGHT_BP = 0.70 * 9 * 72
MIN_TICK_PT = 7.5
# Allow a hairline of bbox contact: PDF glyph bboxes include ascender/descender
# padding that adjacent all-caps labels do not actually fill.
OVERLAP_TOLERANCE_BP = 1.0


# Figures with an axis-label overlap that is NOT resolvable by the renderer
# alone, each with the measurement that makes it a blocker rather than a defer.
# Entries are enforced in both directions: excluded from the overlap scan, and
# separately asserted to STILL overlap, so an entry cannot outlive its cause.
# Empty on purpose. FS1's 19-id overlap was the last shipped blocker; it is
# resolved by a3_compose_tall_a (spanning cell, same type scale). A future
# geometric collision goes here only with a measurement, never as a defer.
KNOWN_UNRESOLVED_OVERLAP = set()

# Independent cells have different legends/aspect ratios. Their tag baselines
# must nevertheless share the declared row, and each column retains its inset.
TAG_ROWS = {
    "F1_scale_plateau": ("A", "BC"),
    "F3_clinvar_forest": ("AB", "C", "DE"),
    "F6_protein_uq": ("A", "B", "CD"),
    "F8_developability_map": ("AB",),
    "F11_boundary_atlas": ("AB",),
    "FS1_proteingym_atlas": ("AB", "C"),
    "FS2_gnomad_control": ("AB", "C"),
}


def text_spans(page):
    return [
        span for block in page.get_text("dict")["blocks"]
        for line in block.get("lines", []) for span in line.get("spans", [])
        if span["text"].strip()
    ]


def printed_scale(page):
    return min(PRINT_WIDTH_BP / page.rect.width,
               PRINT_MAX_HEIGHT_BP / page.rect.height)


def bbox_clearance(a, b):
    return math.hypot(max(0, a.x0 - b.x1, b.x0 - a.x1),
                      max(0, a.y0 - b.y1, b.y0 - a.y1))


def axis_segments(page):
    """Long horizontal/vertical strokes, excluding diagonal data/reference lines."""
    for drawing in page.get_drawings():
        if drawing["type"] != "s":
            continue
        for item in drawing["items"]:
            if item[0] != "l":
                continue
            a, b = item[1:]
            dx, dy = abs(a.x - b.x), abs(a.y - b.y)
            if max(dx, dy) < 20 or (dx > 0.01 and dy > 0.01):
                continue
            rect = fitz.Rect(min(a.x, b.x), min(a.y, b.y),
                             max(a.x, b.x), max(a.y, b.y))
            half = drawing["width"] / 2
            yield rect + (-half, -half, half, half)


def horizontal_tick_runs(page):
    """Group left-of-panel horizontal text into columns sharing an x-extent.

    A y-axis tick block is a set of labels that share a right edge (they are
    right-aligned against the axis) and are stacked vertically. Grouping on the
    rounded right edge isolates one axis's tick block from another panel's.
    """
    runs = defaultdict(list)
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            # Horizontal text only; rotated tick labels stack along x instead.
            if tuple(round(v) for v in line.get("dir", (1, 0))) != (1, 0):
                continue
            for span in line.get("spans", []):
                text = span["text"].strip()
                if not text:
                    continue
                rect = fitz.Rect(span["bbox"])
                runs[round(rect.x1)].append((rect, text, span["size"]))
    return runs


def _overlaps(stem):
    """True if any two axis tick labels in this figure collide."""
    page = fitz.open(FIGURES / f"{stem}.pdf")[0]
    for items in horizontal_tick_runs(page).values():
        if len(items) < 3:
            continue
        items.sort(key=lambda item: item[0].y0)
        for (a_rect, _, _), (b_rect, _, _) in zip(items, items[1:]):
            if min(a_rect.x1, b_rect.x1) - max(a_rect.x0, b_rect.x0) <= 0:
                continue
            if a_rect.y1 - b_rect.y0 > OVERLAP_TOLERANCE_BP:
                return True
    return False


class ShippedFigureLegibility(unittest.TestCase):
    def test_every_shipped_figure_pdf_is_present(self) -> None:
        missing = [s for s in SHIPPED if not (FIGURES / f"{s}.pdf").is_file()]
        self.assertEqual(missing, [], f"missing shipped figures: {missing}")

    def test_known_unresolved_overlaps_still_overlap(self) -> None:
        """Strict guard: a known blocker may not be silently fixed or forgotten.

        If one of these figures stops overlapping, this test fails and forces
        the entry to be removed, so KNOWN_UNRESOLVED_OVERLAP cannot rot into a
        permanent excuse the way release-spec's documents.main.sha256 did.
        """
        for stem in KNOWN_UNRESOLVED_OVERLAP:
            self.assertTrue(
                _overlaps(stem),
                f"{stem} no longer overlaps -- remove it from "
                f"KNOWN_UNRESOLVED_OVERLAP and from the renderer comment",
            )

    def test_no_axis_tick_label_overlaps_its_neighbour(self) -> None:
        offenders = []
        for stem in SHIPPED:
            if stem in KNOWN_UNRESOLVED_OVERLAP:
                continue
            page = fitz.open(FIGURES / f"{stem}.pdf")[0]
            for right_edge, items in horizontal_tick_runs(page).items():
                if len(items) < 3:
                    # Two-item groups are usually a title plus a stray label,
                    # not a tick block; requiring three avoids false positives.
                    continue
                items.sort(key=lambda item: item[0].y0)
                for (a_rect, a_text, _), (b_rect, b_text, _) in zip(items, items[1:]):
                    # Two labels collide only if their boxes intersect in BOTH
                    # axes. Checking y alone misreads a horizontal row of x-tick
                    # labels that happen to share a rounded right edge.
                    if min(a_rect.x1, b_rect.x1) - max(a_rect.x0, b_rect.x0) <= 0:
                        continue
                    overlap = a_rect.y1 - b_rect.y0
                    if overlap > OVERLAP_TOLERANCE_BP:
                        offenders.append(
                            f"{stem}: {a_text!r} over {b_text!r} by {overlap:.2f} bp "
                            f"(x1~{right_edge})"
                        )
        self.assertEqual(offenders, [], "overlapping axis labels:\n" + "\n".join(offenders))

    def test_no_numeric_tick_labels_are_rendered_abutting(self) -> None:
        """Catch tick labels crammed together with no separating whitespace.

        The pairwise overlap scan cannot see this case: when adjacent numeric
        tick labels touch, the PDF text extractor merges them into ONE span, so
        there is no pair left to compare. It shows up instead as a single span
        whose text is several decimals run together, e.g. '0.50.60.70.80.91.0'
        -- which is what FS2 panel A produced when its forest was halved in
        width and it kept six default breaks.
        """
        abutting = re.compile(r"\d\.\d+\d\.\d")
        offenders = []
        for stem in SHIPPED:
            page = fitz.open(FIGURES / f"{stem}.pdf")[0]
            for block in page.get_text("dict")["blocks"]:
                for line in block.get("lines", []):
                    for span in line.get("spans", []):
                        text = span["text"].strip()
                        if abutting.search(text):
                            offenders.append(f"{stem}: {text!r}")
        self.assertEqual(offenders, [], "abutting tick labels:\n" + "\n".join(offenders))

    def test_no_shipped_glyph_prints_below_the_legibility_floor(self) -> None:
        offenders = []
        for stem in SHIPPED:
            page = fitz.open(FIGURES / f"{stem}.pdf")[0]
            for block in page.get_text("dict")["blocks"]:
                for line in block.get("lines", []):
                    for span in line.get("spans", []):
                        printed_pt = span["size"] * printed_scale(page)
                        if span["text"].strip() and printed_pt < MIN_TICK_PT:
                            offenders.append(
                                f"{stem}: {span['text'].strip()!r} at {printed_pt:.2f} printed pt"
                            )
        self.assertEqual(offenders, [], "sub-floor type:\n" + "\n".join(offenders))

    def test_fonts_are_embedded_arial_and_only_panel_letters_are_bold(self) -> None:
        offenders = []
        for stem in SHIPPED:
            with fitz.open(FIGURES / f"{stem}.pdf") as document:
                for font in document[0].get_fonts(full=True):
                    if "Arial" not in font[3] or not document.extract_font(font[0])[3]:
                        offenders.append(f"{stem}: unembedded/non-Arial font {font[3]}")
                for span in text_spans(document[0]):
                    if "Arial" not in span["font"]:
                        offenders.append(f"{stem}: non-Arial text {span['text']!r}")
                    if "Bold" in span["font"] and span["text"].strip() not in tuple("ABCDEF"):
                        offenders.append(f"{stem}: non-panel bold text {span['text']!r}")
        self.assertEqual(offenders, [], "font-role defects:\n" + "\n".join(offenders))

    def test_panel_letter_rows_and_column_insets_are_aligned(self) -> None:
        for stem in SHIPPED:
            page = fitz.open(FIGURES / f"{stem}.pdf")[0]
            tags = {s["text"].strip(): s for s in text_spans(page)
                    if "Bold" in s["font"] and s["text"].strip() in tuple("ABCDEF")}
            rows = TAG_ROWS.get(stem, ("AB", "CD"))
            self.assertEqual(set(tags), set("".join(rows)), stem)
            for row in rows:
                ys = [tags[t]["origin"][1] for t in row]
                self.assertLessEqual((max(ys) - min(ys)) * printed_scale(page), 0.25, stem)
            xs = [s["origin"][0] for s in tags.values()]
            for x in xs:
                edge = min(xs) if x < page.rect.width / 2 else max(xs)
                self.assertLessEqual(abs(x - edge) * printed_scale(page), 0.25, stem)

    def test_panel_letters_have_net_text_and_axis_clearance(self) -> None:
        for stem in SHIPPED:
            page = fitz.open(FIGURES / f"{stem}.pdf")[0]
            spans = text_spans(page)
            tags = [s for s in spans if "Bold" in s["font"]
                    and s["text"].strip() in tuple("ABCDEF")]
            axes = list(axis_segments(page))
            self.assertTrue(axes, f"{stem}: no measurable axis strokes")
            for tag in tags:
                rect = fitz.Rect(tag["bbox"])
                text_gap = min(bbox_clearance(rect, fitz.Rect(s["bbox"]))
                               for s in spans if s is not tag) * printed_scale(page)
                axis_gap = min(bbox_clearance(rect, s) for s in axes) * printed_scale(page)
                self.assertGreaterEqual(text_gap, 6.0, f"{stem} {tag['text']}: text clearance")
                self.assertGreaterEqual(axis_gap, 6.0, f"{stem} {tag['text']}: axis clearance")

    def test_no_text_is_outside_the_current_pdf_canvas(self) -> None:
        for stem in SHIPPED:
            page = fitz.open(FIGURES / f"{stem}.pdf")[0]
            bounds = page.rect + (-0.25, -0.25, 0.25, 0.25)
            for span in text_spans(page):
                self.assertTrue(bounds.contains(fitz.Rect(span["bbox"])),
                                f"{stem}: out-of-canvas text {span['text']!r}")

    def test_fig2_panel_a_names_all_twenty_genes(self) -> None:
        """Regression lock for the twice-deferred Fig 2 panel A defect.

        Panel A is a twenty-gene ClinVar forest. It was resolved by splitting
        the genes across two sub-columns of ten; the failure mode to guard is a
        future 'fix' that restores legibility by dropping gene labels.
        """
        import csv

        genes = [
            row["gene"]
            for row in csv.DictReader(
                (FIGURES / "plotdata" / "F3_clinvar.csv").open(encoding="utf-8")
            )
        ]
        self.assertEqual(len(genes), 20)
        text = fitz.open(FIGURES / "F3_clinvar_forest.pdf")[0].get_text()
        missing = [gene for gene in genes if gene not in text]
        self.assertEqual(missing, [], f"gene labels dropped from Fig 2 panel A: {missing}")

    def test_fs1_panel_a_names_all_nineteen_assays(self) -> None:
        """Tall-A forest must keep every ProteinGym identifier, not drop labels."""
        import csv
        import re

        assays = []
        for row in csv.DictReader((FIGURES / "plotdata" / "FS1_proteingym_atlas.csv").open(encoding="utf-8")):
            short = re.sub(r"_[0-9]{4}$", "", re.sub(r"_[A-Za-z]+$", "", row["assay"]))[:18]
            assays.append(short)
        self.assertEqual(len(assays), 19)
        text = fitz.open(FIGURES / "FS1_proteingym_atlas.pdf")[0].get_text()
        missing = [assay for assay in assays if assay not in text]
        self.assertEqual(missing, [], f"assay labels dropped from FS1 panel A: {missing}")


if __name__ == "__main__":
    unittest.main()
