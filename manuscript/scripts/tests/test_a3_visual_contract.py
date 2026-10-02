"""Shipped-figure visual contract for the A3 compositor rewrite.

Measures manuscript/figures/*.pdf — the artifacts the reader receives —
not the August task-6 stage and not synthetic sidecars.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import fitz

import importlib.util

ROOT = Path(__file__).resolve().parents[2]
FIGURES = ROOT / "figures"
F1 = FIGURES / "F1_scale_plateau.pdf"

spec = importlib.util.spec_from_file_location(
    "a3_figure_release", ROOT / "scripts" / "verify_figure_release.py"
)
contract = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(contract)
PANEL_COUNTS = contract.PANEL_COUNTS

SHIPPED = (
    "F1_scale_plateau", "F2_cross_model", "F3_clinvar_forest",
    "F4_clinical_calibration", "F5_stability_ddg", "F6_protein_uq",
    "F7_binding_boundary", "F8_developability_map", "F9_antibody_lm",
    "F10_reliability_predictor", "F11_boundary_atlas",
    "FS1_proteingym_atlas", "FS2_gnomad_control", "FS3_antibody_robustness",
    "FS4_abstention_tiers", "FS5_vep_comparison",
)

# Authored scale at the printed 6.5-inch width (theme_a3.R).
MIN_PRINT_PT = 7.0
AXIS_PT = 9.0
TICK_PT = 8.0
LABEL_PT = 11.0


def _page(stem: str):
    return fitz.open(FIGURES / f"{stem}.pdf")[0]


def _tags(page):
    tags = [(w[4], fitz.Rect(w[:4])) for w in page.get_text("words") if w[4] in set("ABCDEF")]
    return tags


def test_every_canonical_figure_has_row_major_letters_outside_ggplot_tags():
    theme = (ROOT / "scripts" / "R" / "theme_a3.R").read_text(encoding="utf-8")
    compose = (ROOT / "scripts" / "R" / "compose_a3.R").read_text(encoding="utf-8")
    assert "plot.tag = element_blank()" in theme
    assert "a3_labelled_cell" in compose
    assert 'plot.tag.position = "topleft"' not in theme
    for stem in SHIPPED:
        page = _page(stem)
        tags = _tags(page)
        letters = [t for t, _ in tags]
        n = PANEL_COUNTS.get(stem, 6)
        assert letters.count("A") == 1 and len(tags) == n, (stem, letters)
        expected = list("ABCDEF"[:n])
        ordered = contract.panel_labels(FIGURES / f"{stem}.pdf")
        assert ordered == expected, (stem, ordered)


def test_fs1_tall_first_letters_keep_reserved_bands():
    """Current FS1: A spans the left column, B/C stack at the right."""
    page = _page("FS1_proteingym_atlas")
    by = {letter: rect for letter, rect in _tags(page)}
    assert set(by) == set("ABC")
    mid_x = page.rect.width * 0.5
    assert by["A"].x0 < mid_x and by["A"].y0 < 16
    assert by["B"].x0 > mid_x and by["C"].x0 > mid_x
    assert abs(by["B"].y0 - by["A"].y0) < 1
    assert by["C"].y0 > by["B"].y0 + page.rect.height * .4
    for letter, rect in by.items():
        assert rect.y1 - rect.y0 < 14, (letter, rect)


def test_panel_letters_sit_in_the_reserved_band_not_on_spines():
    """A--F must sit above the panel body, not on the y-spine or in the data box."""
    page = _page("F1_scale_plateau")
    tags = _tags(page)
    assert tags
    # Letters occupy the top ~14 pt of each 157 pt cell.
    for letter, rect in tags:
        cell_top = 0 if rect.y0 < 80 else page.rect.height / 2
        assert rect.y0 - cell_top < 16, (letter, rect.y0, cell_top)
        assert rect.y1 - cell_top < 28, (letter, rect)


def test_f1_type_scale_is_the_contracted_small_face():
    page = _page("F1_scale_plateau")
    sizes = defaultdict(list)
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line.get("spans", []):
                text = span["text"].strip()
                if text:
                    sizes[round(span["size"], 1)].append(text)
    assert any(abs(s - LABEL_PT) < 0.3 for s in sizes), sizes.keys()
    assert any(abs(s - AXIS_PT) < 0.3 for s in sizes), sizes.keys()
    assert not any(s < MIN_PRINT_PT - 0.05 for s in sizes), sizes.keys()
    titles = [s for s, texts in sizes.items() if "protein-mean AUROC" in texts]
    assert titles and max(titles) <= AXIS_PT + 0.3, titles
    ticks = [s for s, texts in sizes.items() if "35M" in texts]
    assert ticks and max(ticks) <= TICK_PT + 0.3, ticks


def test_f1_ylab_sits_next_to_its_own_axis_not_a_sibling_gutter():
    """Independent grobs: each ylab's right edge is close to its own ticks.

    The patchwork over-alignment defect put every column's ylab on one x
    and left a blank band between a short title and its axis. After the
    compositor rewrite the two left-column titles (A and D) must each sit
    next to their own ticks, so their title-to-tick gaps stay small and
    need not be identical.
    """
    page = _page("F1_scale_plateau")
    titles = []
    ticks = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            direction = tuple(round(v) for v in line.get("dir", (1, 0)))
            for span in line.get("spans", []):
                text = span["text"].strip()
                rect = fitz.Rect(span["bbox"])
                if text in {"protein-mean AUROC"} and direction != (1, 0):
                    titles.append((text, rect, direction))
                if text in {"0.66", "0.70", "0.75", "0.00", "0.50", "1.00", "0.0", "0.5", "1.0"}:
                    if rect.x0 < 80:
                        ticks.append((text, rect))
    assert len(titles) >= 1, titles
    # Rotated y-titles have dir (0, -1) or similar; either way the gap
    # from title to the nearest left-column tick must stay under 18 pt.
    gaps = []
    for text, trect, _ in titles:
        near = [krect for _, krect in ticks if abs(krect.y0 - trect.y0) < 80 or abs(krect.y1 - trect.y1) < 80]
        if not near:
            continue
        gap = min(krect.x0 - trect.x1 for krect in near)
        gaps.append((text, gap))
        assert 0 <= gap < 18, (text, gap, trect, near[:3])
    assert len(gaps) >= 1, gaps
