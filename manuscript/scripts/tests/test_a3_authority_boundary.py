"""Authority boundary tests for the A3 canonical figure spine."""

from __future__ import annotations

import importlib.util
import re
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANUSCRIPT = ROOT / "manuscript.tex"
LEDGER = ROOT / "figure_ledger.md"
PROVENANCE = ROOT / "FIGURE-PROVENANCE.md"
STATUS = ROOT / "figures" / "plotdata" / "F7_20gene_perbin_STATUS.txt"
BUILD_ASSETS = ROOT / "scripts" / "build_release_assets.sh"
BUILD_SUBMISSION = ROOT / "scripts" / "build_submission_pdf.sh"
VERIFY = ROOT / "scripts" / "verify_figure_release.py"

spec = importlib.util.spec_from_file_location("a3_figure_release", VERIFY)
contract = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(contract)


class A3AuthorityBoundaryTests(unittest.TestCase):
    def test_reader_facing_main_sequence_is_explicit_and_separate_from_logical_ids(self) -> None:
        self.assertEqual(
            contract.MAIN_READER_SEQUENCE,
            (
                (1, "F1_scale_plateau"),
                (2, "F3_clinvar_forest"),
                (3, "F4_clinical_calibration"),
                (4, "F5_stability_ddg"),
                (5, "F7_binding_boundary"),
                (6, "F8_developability_map"),
                (7, "F11_boundary_atlas"),
            ),
        )
        ledger = LEDGER.read_text(encoding="utf-8")
        self.assertIn("Reader-facing number", ledger)
        for reader_number, logical_id in contract.MAIN_READER_SEQUENCE:
            self.assertIn(f"| {reader_number} | `{logical_id}` |", ledger)
        canonical_source_stems = re.findall(
            r"\\includegraphics\[[^]]*\]\{figures/([^}]+)\.pdf\}",
            MANUSCRIPT.read_text(encoding="utf-8"),
        )
        self.assertEqual(
            canonical_source_stems[: len(contract.MAIN_READER_SEQUENCE)],
            [logical_id for _, logical_id in contract.MAIN_READER_SEQUENCE],
        )

    def test_candidate_is_not_canonical_source_or_build_output(self) -> None:
        candidate = "F7_20gene_ece"
        self.assertNotIn(candidate, contract.STEMS)
        self.assertEqual(contract.CANDIDATE_STEMS, (candidate,))
        self.assertNotIn(candidate, MANUSCRIPT.read_text(encoding="utf-8"))
        self.assertNotIn("make_f7_20gene_ece.R", BUILD_ASSETS.read_text(encoding="utf-8"))
        self.assertNotIn("make_f7_20gene_ece.R", BUILD_SUBMISSION.read_text(encoding="utf-8"))

    def test_canonical_main_source_caption_bindings_and_subjects_match_reader_numbers(self) -> None:
        expected = (
            ("F1_scale_plateau", "fig1", "Fig 1.", "Discrimination plateaus"),
            ("F3_clinvar_forest", "fig2", "Fig 2.", "ClinVar discrimination"),
            ("F4_clinical_calibration", "fig3", "Fig 3.", "calibration subset"),
            ("F5_stability_ddg", "fig5", "Fig 4.", "Stability"),
            ("F7_binding_boundary", "fig6", "Fig 5.", "binding ranker"),
            ("F8_developability_map", "fig7", "Fig 6.", "No BH-significant"),
            ("F11_boundary_atlas", "fig8", "Fig 7.", "Boundary atlas"),
        )
        bindings = re.findall(
            r"\\includegraphics\[[^]]*\]\{figures/([^}]+)\.pdf\}\s*\\input\{captions/([^}]+)\}",
            MANUSCRIPT.read_text(encoding="utf-8"),
        )
        self.assertEqual(bindings[: len(expected)], [(stem, caption) for stem, caption, _, _ in expected])
        for _, caption, title, subject in expected:
            caption_text = (ROOT / "captions" / f"{caption}.tex").read_text(encoding="utf-8")
            self.assertIn(title, caption_text)
            self.assertIn(subject, caption_text)

    def test_stop_status_blocks_candidate_from_canonical_promotion(self) -> None:
        status = STATUS.read_text(encoding="utf-8")
        self.assertIn("VERDICT: STOP. Full-20 per-bin coordinates are not recoverable from capsule files.", status)
        self.assertEqual(contract.perbin_candidate_state(STATUS), "candidate_blocked")
        self.assertFalse(contract.is_canonical_stem("F7_20gene_ece"))
        self.assertTrue(contract.is_canonical_stem("F7_binding_boundary"))

    def test_missing_stop_verdict_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory(prefix="a3-perbin-status-") as temporary:
            status = Path(temporary) / "status.txt"
            status.write_text("VERDICT: RECOVERABLE\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "required STOP verdict"):
                contract.perbin_candidate_state(status)

    def test_atlas_is_registered_only_as_a_boundary_synthesis(self) -> None:
        for document in (LEDGER, PROVENANCE):
            text = document.read_text(encoding="utf-8").lower()
            self.assertIn("atlas-not-law", text)
            self.assertIn("synthesis", text)
            self.assertIn("not an independent data axis", text)


if __name__ == "__main__":
    unittest.main()
