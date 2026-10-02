import csv
import re
import unittest
from pathlib import Path


def expected_panels(rows):
    return {figure: {row["panel_id"] for row in rows if row["figure_id"] == figure} for figure in {row["figure_id"] for row in rows}}

def parse_caption_panels(caption_text):
    found = {}
    for figure, body in re.findall(r"(?:Figure|Fig\.?)\s*(S?\d+)\s*[:.](.*?)(?=(?:Figure|Fig\.?)\s*S?\d+\s*[:.]|\Z)", caption_text, re.I | re.S):
        labels = set(re.findall(r"(?<![A-Za-z])[A-H](?![A-Za-z])", body))
        found[f"F{figure}" if not figure.startswith("S") else figure] = labels
    return found

def caption_errors(expected, observed):
    errors = []
    for figure, panels in expected.items():
        if figure not in observed: errors.append(f"{figure} caption missing")
        elif observed[figure] != panels: errors.append(f"{figure} caption panel set mismatch")
    return errors

class CaptionPanelCoverageAudit(unittest.TestCase):
    def test_current_caption_audit_is_explicit_expected_red(self):
        root = Path(__file__).resolve().parents[2]
        panels_path = root / "registries" / "figure_panels.csv"
        captions = list((root / "captions").glob("*")) if (root / "captions").is_dir() else []
        with panels_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        caption_text = "\n".join(path.read_text(encoding="utf-8", errors="ignore") for path in captions if path.is_file())
        self.assertTrue(caption_errors(expected_panels(rows), parse_caption_panels(caption_text)), "caption audit is expected RED until caption tasks complete")

    def test_caption_fixture_rejects_missing_panels_and_mismatched_figure(self):
        expected = {"F1": set("ABCDEF")}
        self.assertTrue(caption_errors(expected, parse_caption_panels("Figure 1: (A) a (B) b (C) c (D) d (E) e")))
        self.assertTrue(caption_errors(expected, parse_caption_panels("Figure 2: (A) a (B) b (C) c (D) d (E) e (F) f")))

    def test_valid_literal_caption_parses_without_errors(self):
        self.assertEqual(caption_errors({"F1": set("ABCDEF")}, parse_caption_panels("Figure 1: (A) alpha (B) beta (C) gamma (D) delta (E) epsilon (F) zeta")), [])


if __name__ == "__main__":
    unittest.main()
