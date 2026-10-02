"""Task-7 source contract for the frozen A3 six-panel manuscript integration."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CAPTIONS = ROOT / "captions"
PANELS = ROOT / "registries" / "figure_panels.csv"
ESTIMANDS = ROOT / "registries" / "figure_estimands.csv"
METHODS = ROOT / "text" / "methods.tex"
BUILD = ROOT / "scripts" / "build_submission_pdf.sh"
MANUSCRIPT = ROOT / "manuscript.tex"
VISUAL_SOURCE = ROOT / "visual-review-source.json"

FIGURE_CAPTIONS = {
    "F1": "fig1.tex", "F3": "fig2.tex", "F4": "fig3.tex", "F5": "fig5.tex",
    "F7": "fig6.tex", "F8": "fig7.tex", "F11": "fig8.tex",
    "S1": "figs1.tex", "S2": "figs2.tex", "S3": "figs3.tex", "S4": "figs4.tex",
    "S5": "figs5.tex", "S6": "figs6.tex", "S7": "figs7.tex", "S8": "figs8.tex",
    "S9": "figs9.tex",
}

PANEL_SPECIFIC_EF = {
    "F1": ("A3-EQ-F1-PAIRED-01", "A3-EQ-F1-SCORE-01"),
    "F3": ("A3-EQ-F3-SIZE-01", "A3-EQ-F3-LOGO-01"),
    "F4": ("A3-EQ-F4-P53-01", "A3-EQ-F4-COVERAGE-01"),
    "F5": ("A3-EQ-F5-DISC-01", "A3-EQ-F5-CAL-01"),
    "F7": ("A3-EQ-F7-SUBGROUP-01", "A3-EQ-F7-WIDTH-01"),
    "F8": ("A3-EQ-F8-EFFECT-01", "A3-EQ-F8-MDE-01"),
    "F11": ("A3-EQ-F11-SEPARATION-01", "A3-EQ-F11-PROVENANCE-01"),
    "S1": ("A3-EQ-S1-CI-01", "A3-EQ-S1-EXPOSURE-01"),
    "S2": ("A3-EQ-S2-GENE-01", "A3-EQ-S2-CLASS-01"),
    "S3": ("A3-EQ-S3-PAIRED-01", "A3-EQ-S3-INFLUENCE-01"),
    "S4": ("A3-EQ-S4-RISK-01", "A3-EQ-S4-THRESHOLD-01"),
    "S5": ("A3-EQ-S5-RANGE-01", "A3-EQ-S5-MISSING-01"),
    "S6": ("A3-EQ-S6-PAIR-01", "A3-EQ-S6-DIRECTION-01"),
    "S7": ("A3-EQ-S7-SHIFT-01", "A3-EQ-S7-NEGATIVE-01"),
    "S8": ("A3-EQ-S8-EFFECT-01", "A3-EQ-S8-FDR-01"),
    "S9": ("A3-EQ-S9-PARTIAL-01", "A3-EQ-S9-MISSING-01"),
}


class Task7ManuscriptContract(unittest.TestCase):
    def test_atlas_antibody_interface_row_is_cited_to_the_locked_n4_cell(self) -> None:
        """Fig 7B 0.165 is the N4 intersection cell, not a Table 2 / F4 partition."""
        f11 = list(csv.DictReader((ROOT / "figures" / "plotdata" / "F11_atlas_spearman.csv").open(encoding="utf-8")))
        f4 = list(csv.DictReader((ROOT / "figures" / "plotdata" / "F4_binding.csv").open(encoding="utf-8")))
        row = next(item for item in f11 if item["task"] == "Binding ΔΔG, antibody interface")
        self.assertEqual((row["value"], row["lo"], row["hi"]), ("0.165", "0.075", "0.249"))
        self.assertFalse(any(abs(float(item["spearman"]) - 0.165) < 1e-9 for item in f4))

        parent = ROOT.parent / "results-snapshot" / "N4_antibody_ddg" / "results.json"
        self.assertTrue(parent.is_file(), parent)
        n4 = json.loads(parent.read_text(encoding="utf-8"))["antibody_interface"]
        self.assertEqual(int(n4["n"]), 525)
        self.assertEqual(round(float(n4["rho"]), 3), 0.165)
        self.assertEqual(round(float(n4["ci_lo"]), 3), 0.075)
        self.assertEqual(round(float(n4["ci_hi"]), 3), 0.249)

        body = MANUSCRIPT.read_text(encoding="utf-8")
        caption = (CAPTIONS / "fig8.tex").read_text(encoding="utf-8")
        methods = METHODS.read_text(encoding="utf-8")
        for text in (body, caption, methods):
            self.assertIn("0.165", text)
            self.assertIn("F11", text.replace("\\_", "_"))
            self.assertIn("525", text)
            self.assertIn("not a Table 2", text)
            self.assertNotIn("F7_20gene_ece", text)

    def test_full20_headline_does_not_point_to_a_noncanonical_figure(self) -> None:
        """The full-20 token is textual; canonical Fig. 4 is stability, not calibration."""
        source = MANUSCRIPT.read_text(encoding="utf-8")
        self.assertNotRegex(
            source,
            r"full 20-gene panel[^\n]*Fig 4A--C",
        )

    def test_f4_full20_headline_has_direct_locked_equation_bindings(self) -> None:
        """Full-20 summaries are traceable without misrepresenting six-gene bins."""
        with (ROOT / "registries" / "figure_inputs.sha256").open(encoding="utf-8", newline="") as handle:
            inputs = {row["input_id"]: row for row in csv.DictReader(handle)}
        with ESTIMANDS.open(encoding="utf-8", newline="") as handle:
            estimands = {row["equation_id"]: row for row in csv.DictReader(handle)}
        required = {
            "A3-EQ-F4-FULL20-ECE-01": ("A3-F4-FULL20-ECE", "figures/plotdata/F7_20gene_ece_summary.csv"),
            "A3-EQ-F4-FULL20-COVERAGE-01": ("A3-F4-FULL20-COVERAGE", "figures/plotdata/F7_20gene_coverage.csv"),
        }
        for equation_id, (input_id, source_path) in required.items():
            with self.subTest(equation_id=equation_id):
                self.assertIn(input_id, inputs)
                self.assertEqual(inputs[input_id]["path"], source_path)
                self.assertEqual(inputs[input_id]["status"], "LOCKED")
                self.assertIn(equation_id, estimands)
                self.assertEqual(estimands[equation_id]["input"], source_path)
                self.assertEqual(estimands[equation_id]["output"], "manuscript.tex")
                self.assertRegex(estimands[equation_id]["input_sha256"], r"^[0-9a-f]{64}$")
                self.assertRegex(estimands[equation_id]["output_sha256"], r"^[0-9a-f]{64}$")

        with PANELS.open(encoding="utf-8", newline="") as handle:
            panels = list(csv.DictReader(handle))
        f4_primary = [row for row in panels if row["figure_id"] == "F4" and row["panel_id"] in "ABCD"]
        self.assertTrue(f4_primary)
        self.assertTrue(all(row["input_id"] == "A3-F4-PRIMARY" for row in f4_primary))
        self.assertTrue(all(inputs[row["input_id"]]["path"] == "figures/plotdata/F7_reliability.csv" for row in f4_primary))
        self.assertTrue(all("without recovering per-bin coordinates" in estimands[equation_id]["estimand"] for equation_id in required))

    def test_new_xfit_calibration_is_printed_beside_the_historical_stop_record(self) -> None:
        """A3-XFIT-CAL-20260918 is a new cross-fit, not a restored 0.0207 pair."""
        status = (ROOT / "figures" / "plotdata" / "F7_20gene_perbin_STATUS.txt").read_text(encoding="utf-8")
        self.assertIn("VERDICT: STOP. Full-20 per-bin coordinates are not recoverable from capsule files.", status)
        xfit = list(csv.DictReader((ROOT / "figures" / "plotdata" / "F7_20gene_xfit_reliability.csv").open(encoding="utf-8")))
        self.assertGreaterEqual(len(xfit), 20)
        self.assertTrue({"calibrator", "mean_pred", "obs_rate", "n"} <= set(xfit[0]))
        self.assertTrue(any(row["calibrator"] == "platt" for row in xfit))
        body = MANUSCRIPT.read_text(encoding="utf-8")
        methods = METHODS.read_text(encoding="utf-8")
        for text in (body, methods):
            self.assertIn("A3-XFIT-CAL-20260918", text)
            self.assertIn("0.019", text)
            self.assertNotIn("0.0207", text)
            self.assertNotIn("91%", text)
        ece = list(csv.DictReader((ROOT / "figures" / "plotdata" / "F7_20gene_xfit_ece.csv").open(encoding="utf-8")))
        platt = next(row for row in ece if row["panel"] == "full20_gof" and row["calibrator"] == "platt")
        self.assertAlmostEqual(float(platt["ece"]), 0.019244, places=5)
        genes = list(csv.DictReader((ROOT / "figures" / "plotdata" / "F7_20gene_xfit_per_gene.csv").open(encoding="utf-8")))
        worse = [
            row["group"]
            for row in genes
            if row["ece_platt"] and row["ece_uncalibrated"]
            and float(row["ece_platt"]) > float(row["ece_uncalibrated"])
        ]
        self.assertIn("P53_HUMAN", worse)
        self.assertIn("OTC_HUMAN", worse)

    def test_xfit_round2_mechanism_is_printed_and_not_a_restoration(self) -> None:
        """A3-XFIT-CAL-R2-20260918 explains 7/20 degradation; not 0.0207."""
        status = (ROOT / "figures" / "plotdata" / "F7_20gene_xfit_r2_STATUS.txt").read_text(encoding="utf-8")
        self.assertIn("A3-XFIT-CAL-R2-20260918", status)
        self.assertIn("prevalence mismatch", status)
        self.assertIn("NOT a restoration", status)
        genes = list(csv.DictReader((ROOT / "figures" / "plotdata" / "F7_20gene_xfit_r2_per_gene.csv").open(encoding="utf-8")))
        self.assertEqual(len(genes), 20)
        required = {"gene", "n", "base_rate", "ece_uncal", "ece_platt", "ece_isotonic", "delta_platt", "p_holm", "q_bh"}
        self.assertTrue(required <= set(genes[0]))
        holm_deg = [
            row["gene"]
            for row in genes
            if row["degrades_pooled_platt"] == "True" and row["reject_holm_05"] == "True"
        ]
        self.assertEqual(sorted(holm_deg), ["OTC_HUMAN", "PIK3CA_GOF"])
        p53 = next(row for row in genes if row["gene"] == "P53_HUMAN")
        self.assertEqual(p53["degrades_pooled_platt"], "True")
        self.assertEqual(p53["reject_holm_05"], "False")
        schemes = list(csv.DictReader((ROOT / "figures" / "plotdata" / "F7_20gene_xfit_r2_schemes.csv").open(encoding="utf-8")))
        pergene = next(row for row in schemes if row["scheme"] == "pergene_platt")
        pooled = next(row for row in schemes if row["scheme"] == "pooled_platt")
        self.assertAlmostEqual(float(pooled["ece"]), 0.019244, places=5)
        self.assertAlmostEqual(float(pergene["ece"]), 0.014590, places=5)
        self.assertGreater(float(pergene["delta_vs_pooled_p"]), 0.05)
        slopes = list(csv.DictReader((ROOT / "figures" / "plotdata" / "F7_20gene_xfit_r2_slope.csv").open(encoding="utf-8")))
        uncal = next(row for row in slopes if row["scheme"] == "uncalibrated")
        platt = next(row for row in slopes if row["scheme"] == "pooled_platt")
        self.assertAlmostEqual(float(uncal["slope"]), 0.850, places=2)
        self.assertAlmostEqual(float(platt["slope"]), 0.999, places=2)
        body = MANUSCRIPT.read_text(encoding="utf-8")
        methods = METHODS.read_text(encoding="utf-8")
        for text in (body, methods):
            self.assertIn("A3-XFIT-CAL-R2-20260918", text)
            self.assertNotIn("0.0207", text)
            self.assertNotIn("91%", text)
        self.assertIn("0.815", body)
        self.assertIn("Holm", body)

    def test_retained_captions_match_drawn_panels_in_main_and_si(self) -> None:
        manuscript_inputs = set(re.findall(r"\\input\{captions/([^}]+)\}", MANUSCRIPT.read_text(encoding="utf-8")))
        si_inputs = set(re.findall(r"\\input\{captions/([^}]+)\}", (ROOT / "SI.tex").read_text(encoding="utf-8")))
        main_files = {"fig1.tex", "fig2.tex", "fig3.tex", "fig5.tex", "fig6.tex", "fig7.tex", "fig8.tex"}
        si_files = {f"figs{n}.tex" for n in range(1, 10)}
        self.assertEqual(manuscript_inputs & set(FIGURE_CAPTIONS.values()), main_files)
        self.assertEqual(si_inputs & set(FIGURE_CAPTIONS.values()), si_files)
        self.assertNotIn("fig4", FIGURE_CAPTIONS.values(), "unused caption must not become authoritative")
        drawn = {
            "F1": set("ABC"), "F3": set("ABCDE"), "F4": set("ABCDEF"), "F5": set("ABCDEF"),
            "F7": set("ABCDE"), "F8": set("AB"), "F11": set("ABCDEF"),
            "S1": set("ABCDE"), "S2": set("ABCDEF"), "S3": set("ABCDE"), "S4": set("ABCDE"),
            "S5": set("ABCDEF"), "S6": set("ABCDEF"), "S7": set("ABCDEF"), "S8": set("ABCDEF"),
            "S9": set("ABCDEF"),
        }
        for figure_id, filename in FIGURE_CAPTIONS.items():
            text = (CAPTIONS / filename).read_text(encoding="utf-8")
            observed = set(re.findall(r"\(([A-F])\)", text))
            self.assertEqual(observed, drawn[figure_id], figure_id)
        f8 = (CAPTIONS / "fig7.tex").read_text(encoding="utf-8").lower()
        self.assertIn("effect interval", f8)
        self.assertIn("minimum detectable effect", f8)
        self.assertNotIn("associations are absent", f8)
        self.assertNotIn("defined in the figure registry", f8)
        self.assertNotIn("locked E-panel expansion reports", f8)
        self.assertNotIn("locked F-panel expansion reports", f8)

    def test_acknowledgments_do_not_force_references_to_a_new_page(self) -> None:
        text = MANUSCRIPT.read_text(encoding="utf-8")
        self.assertNotRegex(
            text,
            r"\\section\*\{Acknowledgments\}.*?\\clearpage\s*\\begin\{thebibliography\}",
        )

    def test_visual_source_has_only_the_seven_main_and_nine_supplement_figures(self) -> None:
        import json
        rows = json.loads(VISUAL_SOURCE.read_text(encoding="utf-8"))
        figures = [row["item"] for row in rows if row.get("document") == "main"]
        si_figures = [row["item"] for row in rows if row.get("document") == "si"]
        self.assertEqual(figures, [*(f"Fig {n}" for n in range(1, 8))])
        self.assertEqual(si_figures, [*(f"Fig S{n}" for n in range(1, 10))])
        self.assertNotIn("Fig 8", figures)

    def test_locked_derived_supplement_panels_remain_nonpromoted_audits(self) -> None:
        expected = {
            "figs5.tex": "S5B_supervised_group_mean_median_abs_spearman.csv",
            "figs6.tex": "S6D_family_mean_protein_auroc.csv",
            "figs7.tex": "S7A_method_shift_mean_rho_unc.csv",
        }
        with PANELS.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        for filename, output_name in expected.items():
            text = (CAPTIONS / filename).read_text(encoding="utf-8")
            self.assertIn(output_name.replace("_", r"\_"), text)
            self.assertIn("not an independent or promoted measurement", text)
            row = next(row for row in rows if row["output"].endswith(output_name))
            self.assertEqual((row["status"], row["role"], row["claim_level"]), ("BLOCKED_NO_INPUT", "evidence_boundary", "descriptive"))

    def test_retained_panels_bind_real_outputs_or_explicit_evidence_boundaries(self) -> None:
        with PANELS.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        with ESTIMANDS.open(newline="", encoding="utf-8") as handle:
            estimands = {row["equation_id"] for row in csv.DictReader(handle)}
        self.assertEqual({row["figure_id"] for row in rows}, set(FIGURE_CAPTIONS))
        for row in rows:
            self.assertIn(row["estimand_id"], estimands)
            self.assertNotEqual(row["role"], "planned", f"{row['figure_id']}:{row['panel_id']}")
            self.assertNotIn("planned/", row["output"])
            self.assertNotIn("planned:", row["renderer"])
            self.assertIn(row["status"], {"READY", "BLOCKED_NO_INPUT"})

    def test_panel_estimands_have_concrete_hash_bound_payloads(self) -> None:
        banned = ("locked source variables", "registered unit", "registered locked csv", "registered panel output")
        with ESTIMANDS.open(newline="", encoding="utf-8") as handle:
            equations = {row["equation_id"]: row for row in csv.DictReader(handle)}
        with PANELS.open(newline="", encoding="utf-8") as handle:
            panel_equations = {row["estimand_id"] for row in csv.DictReader(handle) if row["panel_id"] in "EF"}
        for equation_id in panel_equations:
            row = equations[equation_id]
            self.assertFalse(any(token in " ".join(row.values()).lower() for token in banned), equation_id)
            for field in ("input", "output", "input_sha256", "output_sha256"):
                self.assertTrue(row[field], f"{equation_id}: {field}")
            for path_field, hash_field in (("input", "input_sha256"), ("output", "output_sha256")):
                path = ROOT / row[path_field]
                self.assertTrue(path.is_file(), f"{equation_id}: {path_field}")
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), row[hash_field], equation_id)

    def test_methods_define_and_registry_uses_required_a3_equation_family(self) -> None:
        methods = METHODS.read_text(encoding="utf-8")
        required = {
            "A3-EQ-MEAN-01", "A3-EQ-DELTA-01", "A3-EQ-AUROC-01", "A3-EQ-SPEARMAN-01",
            "A3-EQ-RMSE-01", "A3-EQ-ECE-01", "A3-EQ-COVERAGE-01", "A3-EQ-BOOT-01",
            "A3-EQ-BH-01", "A3-EQ-POWER-01", "A3-EQ-ATLAS-01", "A3-EQ-F8-POOL-01",
        }
        for equation_id in required:
            self.assertIn(equation_id, methods)
        with ESTIMANDS.open(newline="", encoding="utf-8") as handle:
            registered = {row["equation_id"] for row in csv.DictReader(handle)}
        self.assertTrue(required <= registered)

    def test_isolated_build_requires_hash_checked_external_frozen_figure_overlay(self) -> None:
        build = BUILD.read_text(encoding="utf-8")
        self.assertIn("A3_FROZEN_FIGURE_STAGE", build)
        self.assertIn("frozen figure", build.lower())
        self.assertIn("sha256sum", build)
        self.assertIn("refusing frozen figure stage inside A3 capsule tree", build)
        self.assertNotIn("make_f7_20gene_ece.R", build)


if __name__ == "__main__":
    unittest.main()
