"""Source-level contract for the staged A3 six-panel compositor."""
from pathlib import Path
import csv
import re
import fitz


ROOT = Path(__file__).resolve().parents[2]
RENDERER = ROOT / "scripts" / "R" / "make_four_panel_all.R"
COMPOSE = ROOT / "scripts" / "R" / "compose_a3.R"
EVIDENCE = ROOT / "scripts" / "R" / "evidence_panels_a3.R"
THEME = ROOT / "scripts" / "R" / "theme_a3.R"
SIDECARS = ROOT / "scripts" / "render_a3_stage_sidecars.py"
STAGE = ROOT.parent / ".superpowers" / "sdd" / "three-paper-results-figures-experiments-expansion-implementation-20260824" / "task-6-stage"


def test_canonical_renderer_uses_fail_closed_six_to_eight_panel_compositor():
    source = RENDERER.read_text(encoding="utf-8")
    compose = COMPOSE.read_text(encoding="utf-8")
    assert "save_panels <- function" in source
    assert "exactly 6--8 panels" in source
    assert "a3_compose_readable" in source
    assert "a3_compose_six <- function" in compose
    assert "a3_compose_tall_a <- function" in compose
    assert "save4 <- function" not in source
    assert "A3_EF.csv" in source
    assert 'evidence_panel(ef[1, , drop = FALSE], "E")' in source
    assert 'evidence_panel(ef[2, , drop = FALSE], "F")' in source
    assert 'save_panels(list(p8A, p8D)' in source
    assert 'save_panels(list(p8A, p8B, p8C_plot, p8D)' not in source


def test_second_iteration_uses_true_equal_2x3_cells_and_reader_safe_cards():
    source = RENDERER.read_text(encoding="utf-8")
    compose = COMPOSE.read_text(encoding="utf-8")
    assert "a3_compose_readable" in source
    assert 'layout = "tall_a"' in source
    assert "wrap_plots(" not in source
    assert "unit(1, \"null\"), gutter, unit(1, \"null\"), gutter, unit(1, \"null\")" in compose
    assert "reader_action" in source
    assert "strwrap" in source
    assert '"Action: "' not in source
    assert "replace_with" not in source and "merge_into" not in source and "retire" not in source


def test_theme_reserves_a_shared_header_gutter_and_regular_text_roles():
    theme = THEME.read_text(encoding="utf-8")
    compose = COMPOSE.read_text(encoding="utf-8")
    assert "A3_PANEL_LABEL_PT" in theme
    assert "A3_LABEL_BAND_PT" in theme
    assert "A3_LABEL_BAND_PT" in compose
    assert "axis.title = element_text(size = A3_AXIS_PT, face = \"plain\"" in theme
    # Letters live in compose_a3.R, not in ggplot plot.tag. The theme title
    # keeps a compact inset so the coordinate region keeps the cell.
    assert "plot.tag = element_blank()" in theme
    assert re.search(
        r"plot\.title\s*=.*?margin\s*=\s*margin\(\s*t\s*=\s*0\s*,\s*b\s*=\s*2\s*\)",
        theme,
        re.S,
    )
    assert "fontface = \"bold\"" in compose


def test_next_visual_batch_uses_one_tag_position_and_sparse_safe_axes():
    # One deliberate letter system: compose_a3.R grid textGrobs in a reserved
    # band above each independent panel grob. ggplot plot.tag (topleft / panel
    # / numeric) is forbidden because it sat on spines and crowded y-titles.
    theme = THEME.read_text(encoding="utf-8")
    source = RENDERER.read_text(encoding="utf-8")
    compose = COMPOSE.read_text(encoding="utf-8")
    panels = EVIDENCE.read_text(encoding="utf-8")
    joined = theme + source + compose + panels
    assert "A3_PANEL_TAG_X" not in joined
    assert "plot.tag = element_blank()" in theme
    assert 'plot.tag.position = "topleft"' not in joined
    assert 'plot.tag.location = "panel"' not in joined
    assert "plot.tag.position = c(" not in source + panels + compose
    assert "a3_labelled_cell" in compose
    assert "a3_y_for_labels" in source
    assert "short_estimand" not in source


def test_ready_evidence_cells_are_data_backed_and_fail_closed():
    """A READY sensitivity cell must plot the table its scalar summarises."""
    source = RENDERER.read_text(encoding="utf-8")
    panels = EVIDENCE.read_text(encoding="utf-8")
    assert "evidence_panel_data_backed(row, tag)" in source
    # Fail-closed: pinned digest is re-checked, and the drawn observations must
    # reproduce the committed value before anything is rendered.
    assert "digest::digest(file = path, algo = \"sha256\")" in panels
    assert "evidence panel source digest drift" in panels
    assert "evidence panel disagrees with committed value" in panels
    assert "READY estimand %s has no data-backed evidence panel" in panels


def test_every_ready_estimand_has_a_registered_evidence_panel():
    panels = EVIDENCE.read_text(encoding="utf-8")
    plotdata = ROOT / "figures" / "plotdata"
    ready = set()
    for path in sorted(plotdata.glob("*_EF.csv")):
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row["status"] == "READY":
                    ready.add(row["estimand_id"])
    assert ready, "no READY estimands found"
    registered = set(re.findall(r'"(A3-EQ-[A-Z0-9-]+)"\s*=\s*function', panels))
    assert ready <= registered, f"unregistered READY estimands: {sorted(ready - registered)}"
    assert registered <= ready, f"registry entries with no READY row: {sorted(registered - ready)}"


def test_round_five_typography_contract_is_one_regular_size_and_in_page_tags():
    # 2026-09-09 visual-contract: axis 6.5 / tick 6.0 / letter 8.0, with data
    # marks grown via A3_POINT / A3_LINE / A3_BAR. Letters are compositor
    # grobs, not ggplot tags.
    theme = THEME.read_text(encoding="utf-8")
    source = RENDERER.read_text(encoding="utf-8")
    assert "A3_AXIS_PT <- 8.0" in theme
    assert "A3_TICK_PT <- 7.5" in theme
    assert "A3_TEXT_PT <- 7.5" in theme and "A3_TEXT_MM <- 2.65" in theme
    assert "A3_PANEL_LABEL_PT <- 10" in theme
    assert "A3_POINT <- 3.5" in theme
    assert "A3_LINE <- 1.10" in theme
    assert "A3_BAR <- 0.74" in theme
    assert "geom_text <- function(..., size = A3_TEXT_MM)" in source
    assert "args$size <- size" in source


def test_stage_sidecars_emit_a_machine_readable_clipping_condition():
    source = SIDECARS.read_text(encoding="utf-8")
    assert "text_is_cropped" in source
    assert '"cropped_text": any(' in source
    assert '"font": font' in source and '"direction": list(line["dir"])' in source


def test_sidecar_contract_uses_tag_edge_cell_text_boundaries():
    source = SIDECARS.read_text(encoding="utf-8")
    assert "cell_is_cropped" in source
    assert "cells_from_tag_bboxes" in source
    assert "next tag" in source
    assert '"cell_cropped_text"' in source


def test_sidecar_records_actionable_cell_offenders_and_stage_must_clear_them():
    source = SIDECARS.read_text(encoding="utf-8")
    assert "cell_cropped_words" in source
    assert '"word": x["text"]' in source
    assert '"column": column' in source
    for sidecar in STAGE.glob("*.geometry.sidecar.json"):
        payload = __import__("json").loads(sidecar.read_text(encoding="utf-8"))
        assert all(not page.get("cell_cropped_words", []) for page in payload["pages"])


def test_staged_render_has_six_atomic_cells_and_one_ordered_tag_each():
    """Regression for the prior patchwork-container expansion (A,A,B,C/C,D,E,F)."""
    pdf = STAGE / "F1_scale_plateau.pdf"
    assert pdf.is_file(), "render the Task-6 stage before topology verification"
    words = [(word[4], word[0], word[1]) for word in fitz.open(pdf)[0].get_text("words")]
    tags = [(token, x, y) for token, x, y in words if token in set("ABCDEF")]
    assert len(tags) == 6
    top = [token for token, _, _ in sorted((item for item in tags if item[2] < 0.5 * max(y for _, _, y in tags)), key=lambda item: item[1])]
    bottom = [token for token, _, _ in sorted((item for item in tags if item[2] >= 0.5 * max(y for _, _, y in tags)), key=lambda item: item[1])]
    assert top + bottom == list("ABCDEF")


def test_blocked_cards_extract_complete_concise_copy_without_ellipsis():
    """Blocked panels must communicate their evidence boundary without clipping prose."""
    expected = {
            "F3_clinvar_forest.pdf": (
                "Variant-level audit", "Not estimable", "Missing: score rows",
                "Source: F3 table", "No recovery claim.",
        ),
        "F7_binding_boundary.pdf": (
            "Binding resampling", "Not estimable", "Missing: row coords",
            "Source: F4 table", "No resampling.",
        ),
        "F8_developability_map.pdf": (
            "Effect interval", "Not estimable", "Missing: dataset effects and variance",
            "Source: F5 developability table", "No effect interval.",
            "Power / MDE", "No equivalence conclusion.",
        ),
        "FS1_proteingym_atlas.pdf": (
            "Exposure audit", "Not estimable", "Missing:", "exposure labels",
            "Source: S1 table", "No exposure audit.",
        ),
        "FS3_antibody_robustness.pdf": (
            "Influence audit", "Not estimable", "Missing: row coords",
            "Source: S3 table", "No influence audit.",
        ),
        "FS4_abstention_tiers.pdf": (
            "Tier uncertainty", "Not estimable", "Missing: protein-level scores",
            "Source: S4 tier table", "No tier uncertainty interval.",
        ),
    }
    for filename, phrases in expected.items():
        extracted = " ".join(fitz.open(STAGE / filename)[0].get_text().split())
        for phrase in phrases:
            assert phrase in extracted, (filename, phrase, extracted)
        for incomplete_marker in ("...", "Summary gen", "Aggregate partitio", "Non-claim: Summary gen"):
            assert incomplete_marker not in extracted, (filename, incomplete_marker, extracted)


def test_locked_renderer_aggregates_are_read_directly_and_rendered_as_marks():
    source = RENDERER.read_text(encoding="utf-8")
    for locked_csv in (
        "S6D_family_mean_protein_auroc.csv",
        "S7A_method_shift_mean_rho_unc.csv",
        "S5B_supervised_group_mean_median_abs_spearman.csv",
    ):
        assert f'locked_summary("{locked_csv}"' in source
    assert "f2 %>% group_by(family)" not in source
    assert "ru %>% group_by(method, shift_class)" not in source
    assert "vep %>% group_by(supervised)" not in source
    # F6's locked input is demonstrated through its plotted marks; numerical
    # value labels are intentionally absent to preserve uniform 8.5-pt text.
    assert len(fitz.open(STAGE / "F6_protein_uq.pdf")[0].get_drawings()) > 10


def test_f8c_uses_the_locked_pool_table_not_renderer_aggregation():
    source = RENDERER.read_text(encoding="utf-8")
    assert "locked_f8c()" in source
    f8c = source[source.index('p8C <-'):source.index('p8C_plot <-')]
    assert "group_by" not in f8c and "summarise" not in f8c and "sum(" not in f8c


def test_renderer_uses_shared_ten_point_panel_tags_without_f5_override():
    source = RENDERER.read_text(encoding="utf-8")
    assert "tag_theme <-" not in source
    assert "plot.tag = element_text(size = 12" not in source


def test_actual_span_roles_distinguish_data_annotations_cards_and_rotated_ticks():
    """Real staged PDFs define the semantic classifier contract, not point size."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("stage_sidecars", SIDECARS)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    observed = {}
    for filename, wanted in {
        "F2_cross_model.pdf": {"0.728", "Protein AUROC range"},
        "F9_antibody_lm.pdf": {"3/14", "0/7"},
        "FS3_antibody_robustness.pdf": set(),
        "FS5_vep_comparison.pdf": {"EVE"},
        "F3_clinvar_forest.pdf": {"Missing: score rows", "No recovery claim."},
        "F6_protein_uq.pdf": {"unc"},
    }.items():
        spans = module.spans_for(fitz.open(STAGE / filename)[0], Path(filename).stem, 1)
        for span in spans:
            if span["text"].strip() in wanted:
                observed.setdefault((filename, span["text"].strip()), []).append(span["role"])
    assert "annotation" in observed[("F2_cross_model.pdf", "0.728")]
    assert "annotation" in observed[("F2_cross_model.pdf", "Protein AUROC range")]
    assert "annotation" in observed[("F9_antibody_lm.pdf", "3/14")]
    assert "annotation" in observed[("F9_antibody_lm.pdf", "0/7")]
    assert ("FS3_antibody_robustness.pdf", "650M WT") not in observed
    assert "tick_label" in observed[("FS5_vep_comparison.pdf", "EVE")]
    assert observed[("F3_clinvar_forest.pdf", "Missing: score rows")] == ["blocked_body"]
    assert observed[("F3_clinvar_forest.pdf", "No recovery claim.")] == ["blocked_body"]
    assert observed[("F6_protein_uq.pdf", "unc")] == ["math_subscript", "math_subscript"]


def test_span_role_is_geometry_semantic_not_point_size():
    import importlib.util
    spec = importlib.util.spec_from_file_location("stage_sidecars", SIDECARS)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    page = fitz.open(STAGE / "F2_cross_model.pdf")[0]
    record = next(span for span in module.spans_for(page, "F2_cross_model", 1) if span["text"].strip() == "0.728")
    cells = module.panel_cells(page)
    context = {"cells": cells, "blocked_cells": set(), "ready_cells": set(), "axis": module.axis_spines(page, cells)}
    assert module.role_for(record, context) == module.role_for(record | {"size_pt": 72}, context)
