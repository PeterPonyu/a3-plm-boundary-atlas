#!/usr/bin/env python3
"""Fail-closed validation for staged A3 release figures."""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
import fitz

MAIN_READER_SEQUENCE = (
    (1, "F1_scale_plateau"),
    (2, "F3_clinvar_forest"),
    (3, "F4_clinical_calibration"),
    (4, "F5_stability_ddg"),
    (5, "F7_binding_boundary"),
    (6, "F8_developability_map"),
    (7, "F11_boundary_atlas"),
)

CANDIDATE_STEMS = ("F7_20gene_ece",)
PERBIN_STOP = "VERDICT: STOP. Full-20 per-bin coordinates are not recoverable from capsule files."

STEMS = (
    "F1_scale_plateau", "F2_cross_model", "F3_clinvar_forest",
    "F4_clinical_calibration", "F5_stability_ddg", "F6_protein_uq",
    "F7_binding_boundary", "F8_developability_map",
    "F9_antibody_lm", "F10_reliability_predictor", "F11_boundary_atlas",
    "FS1_proteingym_atlas", "FS2_gnomad_control", "FS3_antibody_robustness",
    "FS4_abstention_tiers", "FS5_vep_comparison", "A3_graphical_abstract",
)
STAGED_CANONICAL_STEMS = tuple(stem for stem in STEMS if stem != "A3_graphical_abstract")


def perbin_candidate_state(status_file: Path) -> str:
    """Return the only permitted state for the unrecoverable full-20 bin data."""
    if PERBIN_STOP in status_file.read_text(encoding="utf-8"):
        return "candidate_blocked"
    raise ValueError("full-20 per-bin status does not carry the required STOP verdict")


def is_canonical_stem(stem: str) -> bool:
    return stem in STEMS


def fail(message: str) -> None:
    raise SystemExit(f"A3 figure-release contract failed: {message}")


def fonts(pdf: Path) -> str:
    completed = subprocess.run(["pdffonts", str(pdf)], text=True, capture_output=True)
    if completed.returncode:
        fail(f"pdffonts failed for {pdf.name}: {completed.stderr.strip()}")
    return completed.stdout


def panel_labels(pdf: Path) -> list[str]:
    """Read unique panel tags in geometric row-major order."""
    document = fitz.open(pdf)
    if document.page_count != 1:
        return []
    tags = [(w[4], w[0], w[1]) for w in document[0].get_text("words") if w[4] in set("ABCDEF")]
    if len({v[0] for v in tags}) != len(tags):
        return []
    return [v[0] for v in sorted(tags, key=lambda v:(round(v[2]/8),v[1]))]

PANEL_COUNTS = {stem:4 for stem in STAGED_CANONICAL_STEMS}
PANEL_COUNTS.update({"F1_scale_plateau":3,"F3_clinvar_forest":5,
    "F4_clinical_calibration":4,"F8_developability_map":2,"F11_boundary_atlas":2,"FS1_proteingym_atlas":3,"FS2_gnomad_control":3})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--staged-canonical-only", action="store_true",
                        help="validate the sixteen rendered figures; the graphical abstract has its own contract")
    args = parser.parse_args()
    output = args.output.resolve()
    if not output.is_dir():
        fail(f"missing staged output directory {output}")
    status = Path(__file__).resolve().parents[1] / "figures" / "plotdata" / "F7_20gene_perbin_STATUS.txt"
    try:
        candidate_state = perbin_candidate_state(status)
    except (OSError, ValueError) as error:
        fail(str(error))
    if candidate_state != "candidate_blocked" or set(STEMS) & set(CANDIDATE_STEMS):
        fail("full-20 candidate must remain blocked and excluded from canonical stems")
    stems = STAGED_CANONICAL_STEMS if args.staged_canonical_only else STEMS
    for stem in stems:
        png, pdf = output / f"{stem}.png", output / f"{stem}.pdf"
        if not png.is_file() or png.stat().st_size < 20_000:
            fail(f"missing or implausibly small PNG for {stem}")
        if not pdf.is_file() or pdf.stat().st_size < 5_000:
            fail(f"missing or implausibly small PDF for {stem}")
        report = fonts(pdf)
        if "Type 3" in report or " no  no " in report:
            fail(f"unembedded or Type 3 font in {stem}")
        if "Arial" not in report:
            fail(f"required Arial font absent from {stem}")
        if stem != "A3_graphical_abstract" and panel_labels(pdf) != list("ABCDEF"[:PANEL_COUNTS.get(stem,6)]):
            fail(f"required row-major panel tags absent, duplicated, or misordered in {stem}")
    print(f"A3 staged figure-release contract passed ({len(stems)} assets, Arial embedded, no Type 3, required panel tags visible).")


if __name__ == "__main__":
    main()
