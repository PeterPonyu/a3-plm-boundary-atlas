from __future__ import annotations

import importlib.util
import hashlib
import os
import subprocess
import tempfile
import time
import unittest
from pathlib import Path
import re
import json


ROOT = Path(__file__).resolve().parents[2]
VERIFY = ROOT / "scripts" / "verify_figure_release.py"
THEME = ROOT / "scripts" / "R" / "theme_a3.R"
BUILD = ROOT / "scripts" / "build_release_assets.sh"
FULL_BUILD = ROOT / "scripts" / "build_submission_pdf.sh"
STAGED_VERIFY = ROOT / "scripts" / "verify_submission_package.py"
STAGE_INPUTS = ROOT / "scripts" / "stage_release_inputs.py"
RELEASE_SPEC = ROOT / "release-spec.json"

spec = importlib.util.spec_from_file_location("a3_contract", VERIFY)
contract = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(contract)

package_spec = importlib.util.spec_from_file_location("a3_package", STAGED_VERIFY)
package_contract = importlib.util.module_from_spec(package_spec)
assert package_spec and package_spec.loader
package_spec.loader.exec_module(package_contract)


class ReleaseContractTests(unittest.TestCase):
    def test_expected_release_asset_set_is_complete_and_unique(self) -> None:
        self.assertEqual(len(contract.STEMS), 17)
        self.assertEqual(len(contract.STEMS), len(set(contract.STEMS)))
        self.assertIn("F8_developability_map", contract.STEMS)
        self.assertIn("F11_boundary_atlas", contract.STEMS)
        self.assertNotIn("F7_20gene_ece", contract.STEMS)

    def test_release_pipeline_requires_explicit_staging_and_strict_tinos(self) -> None:
        theme = THEME.read_text()
        build = BUILD.read_text()
        self.assertIn('stop("A3 release rendering requires the Arial font', theme)
        self.assertIn('A3_FIGURE_OUTPUT', build)
        self.assertIn('realpath -m', build)
        self.assertIn('realpath -e', build)
        self.assertIn('CAPSULE_ROOT', build)
        self.assertIn("trap 'on_signal 130' INT", build)
        self.assertIn("trap 'on_signal 143' TERM", build)
        self.assertIn('canonical figure hashes changed during staging build', build)
        self.assertNotIn('|| true', build)

    def test_direct_release_asset_entrypoint_is_self_contained(self) -> None:
        """The direct figure-only builder must provide its PDF normalizer."""
        with tempfile.TemporaryDirectory(prefix="a3-release-assets-") as temporary:
            output = Path(temporary) / "assets"
            environment = os.environ.copy()
            environment.pop("A3_PDF_NORMALIZER", None)
            completed = subprocess.run(
                ["bash", str(BUILD), str(output)],
                text=True,
                capture_output=True,
                check=False,
                timeout=240,
                env=environment,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue((output / "A3_graphical_abstract.pdf").is_file())
            self.assertTrue((output / "A3_graphical_abstract.png").is_file())

    @staticmethod
    def _figure_digest() -> str:
        digest = hashlib.sha256()
        for path in sorted(ROOT.glob("figures/*")):
            if path.suffix not in {".pdf", ".png"}:
                continue
            digest.update(path.name.encode("utf-8"))
            digest.update(path.read_bytes())
        return digest.hexdigest()

    def test_complete_release_rejects_capsule_traversal_and_symlink_target(self) -> None:
        before = self._figure_digest()
        rejected = ROOT / "must-not-create-release"
        blocked = subprocess.run(
            ["bash", str(FULL_BUILD), str(rejected)], text=True, capture_output=True, check=False
        )
        self.assertNotEqual(blocked.returncode, 0)
        self.assertFalse(rejected.exists())
        with tempfile.TemporaryDirectory(prefix="a3-release-symlink-") as temporary:
            parent = Path(temporary)
            target = parent / "external-target"
            target.mkdir()
            link = parent / "published"
            link.symlink_to(target, target_is_directory=True)
            blocked = subprocess.run(
                ["bash", str(FULL_BUILD), str(link)], text=True, capture_output=True, check=False
            )
            self.assertNotEqual(blocked.returncode, 0)
            self.assertTrue(link.is_symlink())
        self.assertEqual(before, self._figure_digest())

    def test_complete_release_signal_leaves_no_publish_target_or_staging(self) -> None:
        before = self._figure_digest()
        with tempfile.TemporaryDirectory(prefix="a3-release-signal-") as temporary:
            parent = Path(temporary)
            output = parent / "published"
            process = subprocess.Popen(
                ["bash", str(FULL_BUILD), str(output)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
            deadline = time.monotonic() + 30
            while not list(parent.glob(".a3-submission.*")) and time.monotonic() < deadline:
                time.sleep(0.05)
            os.killpg(process.pid, 15)
            process.wait(timeout=30)
            self.assertFalse(output.exists())
            self.assertFalse(list(parent.glob(".a3-submission.*")))
        self.assertEqual(before, self._figure_digest())

    def test_complete_release_post_staging_failure_leaves_no_publish_target_or_staging(self) -> None:
        before = self._figure_digest()
        with tempfile.TemporaryDirectory(prefix="a3-release-post-staging-") as temporary:
            parent = Path(temporary)
            output = parent / "published"
            completed = subprocess.run(
                ["bash", str(FULL_BUILD), str(output)],
                text=True,
                capture_output=True,
                check=False,
                timeout=420,
                env={**os.environ, "A3_RELEASE_TEST_FAIL_AFTER_STAGING": "1"},
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("intentional post-staging validation failure", completed.stderr)
            self.assertFalse(output.exists())
            self.assertFalse(list(parent.glob(".a3-submission.*")))
        self.assertEqual(before, self._figure_digest())

    def test_complete_release_builds_a_full_checked_external_package(self) -> None:
        before = self._figure_digest()
        with tempfile.TemporaryDirectory(prefix="a3-release-full-") as temporary:
            external = Path(temporary) / "external"
            external.mkdir()
            output = external / "published"
            completed = subprocess.run(
                ["bash", str(FULL_BUILD), str(output)], text=True, capture_output=True, check=False, timeout=420
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue((output / "manuscript.pdf").is_file())
            self.assertTrue((output / "SI.pdf").is_file())
            self.assertFalse((output / "manuscript.log").exists())
            self.assertTrue((output / "attestation" / "manuscript.log").is_file())
            self.assertTrue((output / "attestation" / "SI.log").is_file())
            self.assertTrue((output / "attestation" / "log-policy.json").is_file())
            self.assertTrue((output / "attestation" / "build-attestation.json").is_file())
            self.assertTrue((output / "A3_graphical_abstract.pdf").is_file())
            self.assertTrue((output / "A3_graphical_abstract.png").is_file())
            self.assertTrue((output / "visual-registry.json").is_file())
            pages = subprocess.run(
                ["pdfinfo", str(output / "manuscript.pdf")], text=True, capture_output=True, check=True
            ).stdout
            expected_pages = json.loads(RELEASE_SPEC.read_text(encoding="utf-8"))["documents"]["main"]["pages"]
            expected_si_pages = json.loads(RELEASE_SPEC.read_text(encoding="utf-8"))["documents"]["si"]["pages"]
            self.assertRegex(pages, rf"(?m)^Pages:\s+{expected_pages}$")
            si_info = subprocess.run(
                ["pdfinfo", str(output / "SI.pdf")], text=True, capture_output=True, check=True
            ).stdout
            self.assertRegex(si_info, rf"(?m)^Pages:\s+{expected_si_pages}$")
            checksum = subprocess.run(
                ["sha256sum", "-c", "SHA256SUMS.txt"], cwd=output, text=True, capture_output=True, check=False
            )
            self.assertEqual(checksum.returncode, 0, checksum.stderr)
            checksum_paths = [line.split(maxsplit=1)[1] for line in (output / "SHA256SUMS.txt").read_text().splitlines()]
            self.assertEqual(checksum_paths, sorted(package_contract.UPLOAD_PAYLOAD))
            self.assertNotIn("attestation/manuscript.log", checksum_paths)

            registry = json.loads((output / "visual-registry.json").read_text(encoding="utf-8"))
            staged_pdf_hash = hashlib.sha256((output / "manuscript.pdf").read_bytes()).hexdigest()
            main_rows = [row for row in registry if row.get("document") == "main"]
            si_rows = [row for row in registry if row.get("document") == "si"]
            expected_main = [*(f"Fig {n}" for n in range(1, 8))]
            expected_si = [*(f"Fig S{n}" for n in range(1, 10))]
            self.assertEqual([row["item"] for row in main_rows], expected_main)
            self.assertEqual([row["item"] for row in si_rows], expected_si)
            self.assertEqual([row["logical_id"] for row in main_rows], ["F1", "F3", "F4", "F5", "F7", "F8", "F11"])
            self.assertEqual([row["logical_id"] for row in si_rows], [*(f"S{n}" for n in range(1, 10))])
            self.assertNotIn("Fig 8", expected_main)
            self.assertTrue(all(row.get("pdf_sha256") == staged_pdf_hash for row in main_rows))
            staged_si_hash = hashlib.sha256((output / "SI.pdf").read_bytes()).hexdigest()
            self.assertTrue(all(row.get("pdf_sha256") == staged_si_hash for row in si_rows))
            graphical_abstract = next(row for row in registry if row["item"] == "Graphical abstract")
            self.assertEqual(graphical_abstract["document"], "standalone")
            self.assertEqual(
                graphical_abstract["assets"],
                {
                    "A3_graphical_abstract.pdf": hashlib.sha256(
                        (output / "A3_graphical_abstract.pdf").read_bytes()
                    ).hexdigest(),
                    "A3_graphical_abstract.png": hashlib.sha256(
                        (output / "A3_graphical_abstract.png").read_bytes()
                    ).hexdigest(),
                },
            )
            package_spec = json.loads((output / "release-spec.json").read_text(encoding="utf-8"))
            self.assertEqual(package_spec["documents"][0]["sha256"], staged_pdf_hash)
            self.assertEqual(package_spec["documents"][0]["pdf"], "manuscript.pdf")
            self.assertEqual(package_spec["documents"][1]["sha256"], staged_si_hash)
            self.assertEqual(package_spec["documents"][1]["pdf"], "SI.pdf")
            self.assertEqual(package_spec["standalone_assets"], graphical_abstract["assets"])
            spec_graphical_abstract = next(
                item for item in package_spec["visual_items"] if item["id"] == "Graphical abstract"
            )
            self.assertEqual(
                {asset["path"]: asset["sha256"] for asset in spec_graphical_abstract["assets"]},
                graphical_abstract["assets"],
            )
            package_contract.verify_attestation(output)
            attestation_path = output / "attestation" / "build-attestation.json"
            tampered = json.loads(attestation_path.read_text(encoding="utf-8"))
            tampered["document"]["sha256"] = "0" * 64
            attestation_path.write_text(json.dumps(tampered), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "attestation hash mismatch"):
                package_contract.verify_attestation(output)
        self.assertEqual(before, self._figure_digest())

    def test_canonical_upload_checksum_excludes_audit_material_and_attestation_locks_paths(self) -> None:
        flat = ROOT / "submission" / "flat_Proteins"
        checksum_paths = [line.split(maxsplit=1)[1] for line in (flat / "SHA256SUMS.txt").read_text().splitlines()]
        self.assertEqual(checksum_paths, sorted(package_contract.UPLOAD_PAYLOAD))
        self.assertNotIn("manuscript.log", checksum_paths)
        attestation = json.loads((ROOT / "attestation" / "build-attestation.json").read_text(encoding="utf-8"))
        release_spec = json.loads((ROOT / "release-spec.json").read_text(encoding="utf-8"))
        self.assertEqual(release_spec["build_attestation"], "attestation/build-attestation.json")
        # The checked-in attestation must bind the authority PDF to the input
        # closure that produced it.  Comparing the recorded closure digest to
        # the source specification catches a promoted PDF whose figures, theme
        # or renderer have since moved, which a spec-file inequality check
        # cannot: that only reported whether the source file had been touched.
        self.assertEqual(
            attestation["build_inputs_digest"],
            release_spec["build_inputs_digest"],
        )
        self.assertEqual(
            attestation["source_release_spec"]["sha256"],
            hashlib.sha256((ROOT / "release-spec.json").read_bytes()).hexdigest(),
        )
        self.assertEqual(
            attestation["document"]["sha256"],
            hashlib.sha256((flat / "manuscript.pdf").read_bytes()).hexdigest(),
        )
        self.assertEqual(attestation["document"]["pages"], release_spec["documents"]["main"]["pages"])
        self.assertEqual(attestation["si_document"]["pages"], release_spec["documents"]["si"]["pages"])
        self.assertEqual(
            attestation["si_document"]["sha256"],
            hashlib.sha256((flat / "SI.pdf").read_bytes()).hexdigest(),
        )
        self.assertEqual(attestation["package_release_spec"]["path"], "release-spec.json")
        self.assertEqual(
            attestation["package_release_spec"]["sha256"],
            hashlib.sha256((flat / "release-spec.json").read_bytes()).hexdigest(),
        )
        self.assertEqual(attestation["builder"]["path"], "build_submission_pdf.sh")
        self.assertEqual(len(release_spec["package"]["upload_payload"]), len(package_contract.UPLOAD_PAYLOAD))
        self.assertIn("SI.pdf", {entry["path"] for entry in release_spec["package"]["upload_payload"]})
        for entry in release_spec["package"]["upload_payload"]:
            path = flat / entry["path"]
            self.assertFalse(Path(entry["path"]).is_absolute())
            self.assertEqual(entry["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
        for entry in release_spec["standalone_assets"]:
            path = ROOT / entry["path"]
            self.assertFalse(Path(entry["path"]).is_absolute())
            self.assertEqual(entry["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
        all_versioned = release_spec["versioned_paths"]
        self.assertTrue(all((ROOT / path).is_file() and not Path(path).is_absolute() for path in all_versioned))

    def test_source_spec_main_document_hash_binds_the_face_on_disk(self) -> None:
        """documents.main.sha256 must name the PDF that is actually there.

        Added 2026-09-09 after this field was found recording a84b3230, four
        faces behind the promoted 38bbaa0d. It had rotted because nothing wrote
        it and nothing read it: verify_submission_package.py recomputes the
        PACKAGE mirror's sha from measured bytes and reads back only `pages`
        from the source spec, so the source `sha256` was unreachable from both
        directions. "No test binds it" was the cause, not an excuse.

        What this now guarantees: release-spec.json cannot name a main-document
        face that differs from the bytes at documents.main.pdf. It does NOT
        guarantee that face is the promoted one, nor that it is reproducible --
        the attestation and build_inputs_digest checks above cover those.
        """
        release_spec = json.loads((ROOT / "release-spec.json").read_text(encoding="utf-8"))
        main = release_spec["documents"]["main"]
        pdf = ROOT / main["pdf"]
        self.assertTrue(pdf.is_file(), main["pdf"])
        measured = hashlib.sha256(pdf.read_bytes()).hexdigest()
        self.assertEqual(main["sha256"], measured, "documents.main.sha256 is stale")
        info = subprocess.run(
            ["pdfinfo", str(pdf)], text=True, capture_output=True, check=True
        ).stdout
        pages = int(re.search(r"^Pages:\s+(\d+)$", info, re.M).group(1))
        self.assertEqual(main["pages"], pages)
        # pdf_sha256 is a vestigial duplicate of sha256 that nothing reads. It
        # is permitted to stay null, but it may never assert a different face.
        self.assertIn(main.get("pdf_sha256"), (None, measured))

    def test_root_manuscript_pdf_is_byte_identical_to_the_flat_payload(self) -> None:
        """The tracked root PDF must equal the declared deliverable.

        Added 2026-09-09. The root copy was found at 9cb34bd4 while the flat
        payload held 38bbaa0d, and a note in FIGURE-PROVENANCE.md had declared
        the divergence intentional pending "a separately evidenced rebuild".
        That rebuild was performed and reproduced the flat payload byte-for-byte
        from current source, so flat is the reproducible face and root now
        tracks it. Every other capsule holds this invariant; A3 now does too.
        """
        root_pdf = ROOT / "manuscript.pdf"
        flat_pdf = ROOT / "submission/flat_Proteins/manuscript.pdf"
        self.assertTrue(root_pdf.is_file() and flat_pdf.is_file())
        self.assertEqual(
            hashlib.sha256(root_pdf.read_bytes()).hexdigest(),
            hashlib.sha256(flat_pdf.read_bytes()).hexdigest(),
            "manuscript/manuscript.pdf has diverged from the flat_Proteins payload",
        )

    def test_release_spec_authority_paths_are_manuscript_root_relative(self) -> None:
        release_spec = json.loads((ROOT / "release-spec.json").read_text(encoding="utf-8"))
        authority_paths = [
            release_spec["registry"],
            release_spec["build_attestation"],
            *release_spec["versioned_paths"],
            *(document["pdf"] for document in release_spec["documents"].values()),
            *(document["log"] for document in release_spec["documents"].values()),
            *(asset["path"] for asset in release_spec["standalone_assets"]),
            release_spec["package"]["directory"],
            release_spec["package"]["checksum"],
        ]
        for relative_path in authority_paths:
            candidate = Path(relative_path)
            self.assertFalse(candidate.is_absolute(), relative_path)
            self.assertNotIn("..", candidate.parts, relative_path)
            self.assertTrue((ROOT / candidate).resolve().is_relative_to(ROOT.resolve()), relative_path)
        self.assertIn("scripts/tests/test_release_contract.py", release_spec["versioned_paths"])

    def test_hash_bound_build_input_closure_is_nonempty_exact_and_versioned(self) -> None:
        release_spec = json.loads((ROOT / "release-spec.json").read_text(encoding="utf-8"))
        inputs = release_spec["build_inputs"]
        self.assertGreater(len(inputs), 10)
        self.assertEqual(len(inputs), len({row["path"] for row in inputs}))
        self.assertIn("scripts/R/make_four_panel_all.R", {row["path"] for row in inputs})
        self.assertIn("scripts/R/make_f7_20gene_ece.R", {row["path"] for row in inputs})
        self.assertIn("scripts/R/make_a3_graphical_abstract.R", {row["path"] for row in inputs})
        self.assertTrue(all(row["path"] in release_spec["versioned_paths"] for row in inputs))
        self.assertEqual(
            release_spec["build_inputs_digest"],
            package_contract.normalized_build_inputs_digest(inputs),
        )
        with tempfile.TemporaryDirectory(prefix="a3-input-closure-") as temporary:
            stage = Path(temporary) / "clean"
            completed = subprocess.run(
                ["python3", str(STAGE_INPUTS), "--source-root", str(ROOT), "--spec", str(ROOT / "release-spec.json"), "--stage", str(stage)],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(len([path for path in stage.rglob("*") if path.is_file()]), len(inputs) + 1)

    def test_hash_bound_closure_rejects_tampered_renderer_before_publication(self) -> None:
        before = self._figure_digest()
        with tempfile.TemporaryDirectory(prefix="a3-clean-authority-") as temporary:
            clean_root = Path(temporary) / "clean"
            staged = subprocess.run(
                ["python3", str(STAGE_INPUTS), "--source-root", str(ROOT), "--spec", str(ROOT / "release-spec.json"), "--stage", str(clean_root)],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(staged.returncode, 0, staged.stderr)
            renderer = clean_root / "scripts" / "R" / "make_four_panel_all.R"
            renderer.write_text(renderer.read_text(encoding="utf-8") + "\n# tampered\n", encoding="utf-8")
            with tempfile.TemporaryDirectory(prefix="a3-external-publish-") as external:
                output = Path(external) / "published"
                completed = subprocess.run(
                    ["bash", str(clean_root / "scripts" / "build_submission_pdf.sh"), str(output)],
                    text=True, capture_output=True, check=False, timeout=60,
                )
                self.assertNotEqual(completed.returncode, 0)
                self.assertIn("build input hash mismatch", completed.stderr)
                self.assertFalse(output.exists())
        self.assertEqual(before, self._figure_digest())

    def test_complete_release_refuses_a_repeated_publish_without_mutation(self) -> None:
        before = self._figure_digest()
        with tempfile.TemporaryDirectory(prefix="a3-release-repeat-") as temporary:
            output = Path(temporary) / "published"
            first = subprocess.run(
                ["bash", str(FULL_BUILD), str(output)], text=True, capture_output=True, check=False, timeout=420
            )
            self.assertEqual(first.returncode, 0, first.stderr)
            package_digest = hashlib.sha256((output / "SHA256SUMS.txt").read_bytes()).hexdigest()
            second = subprocess.run(
                ["bash", str(FULL_BUILD), str(output)], text=True, capture_output=True, check=False, timeout=60
            )
            self.assertNotEqual(second.returncode, 0)
            self.assertIn("refusing to overwrite output path", second.stderr)
            self.assertEqual(package_digest, hashlib.sha256((output / "SHA256SUMS.txt").read_bytes()).hexdigest())
        self.assertEqual(before, self._figure_digest())

    def test_inline_bibliography_matches_first_citation_order_and_has_metadata(self) -> None:
        def expand(path: Path) -> str:
            text = path.read_text(encoding="utf-8")

            def include(match: re.Match[str]) -> str:
                target = path.parent / f"{match.group(1)}.tex"
                return expand(target)

            return re.sub(r"\\input\{([^}]+)\}", include, text)

        source = expand(ROOT / "manuscript.tex")
        cited = []
        for match in re.finditer(r"\\cite[a-zA-Z*]*\s*(?:\[[^]]*\]\s*)?\{([^}]*)\}", source):
            cited.extend(key.strip() for key in match.group(1).split(",") if key.strip())
        first_cited = list(dict.fromkeys(cited))
        bibliography = re.findall(r"\\bibitem\{([^}]+)\}\s*([^\n]+)", source)
        keys = [key for key, _ in bibliography]
        self.assertEqual(keys, first_cited)
        self.assertEqual(len(keys), len(set(keys)))
        for key, metadata in bibliography:
            self.assertRegex(
                metadata,
                r"(?:doi:|https?://|OpenReview:|\b(?:Nature|Science|Bioinformatics|NeurIPS|ICLR|IEEE|PNAS|PLoS|Nat Genet|Nat Methods|MIT Press|Advances in Neural Information Processing Systems)\b)",
                f"{key} lacks DOI or explicit bibliographic metadata",
            )

        entries = dict(bibliography)
        self.assertIn("36:64331-64379", entries["Notin2023"])
        self.assertIn("Official proceedings record:", entries["Notin2023"])
        self.assertNotIn("10.52202/075280-2810", entries["Notin2023"])
        self.assertIn("23(4):772-784", entries["Dinh2026VESM"])
        self.assertIn("doi:10.1038/s41592-026-03050-9", entries["Dinh2026VESM"])
        self.assertIn("NeurIPS 2021 Datasets and Benchmarks Track", entries["Dallago2021FLIP"])
        self.assertIn("Official proceedings record:", entries["Dallago2021FLIP"])
        self.assertIn("No proceedings DOI is assigned", entries["Dallago2021FLIP"])

    def test_visual_review_source_is_invariant_and_separate_from_published_registry(self) -> None:
        release_spec = json.loads(RELEASE_SPEC.read_text(encoding="utf-8"))
        self.assertEqual(release_spec["registry"], "visual-review-source.json")
        build_paths = {row["path"] for row in release_spec["build_inputs"]}
        self.assertIn("visual-review-source.json", build_paths)
        self.assertNotIn("submission/flat_Proteins/visual-registry.json", build_paths)
        source_registry = json.loads((ROOT / "visual-review-source.json").read_text(encoding="utf-8"))
        self.assertTrue(source_registry)
        self.assertTrue(all("pdf_sha256" not in row and "assets" not in row for row in source_registry))


if __name__ == "__main__":
    unittest.main()
