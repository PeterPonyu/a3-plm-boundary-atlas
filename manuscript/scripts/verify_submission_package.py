#!/usr/bin/env python3
"""Bind, attest, and verify an externally staged A3 upload package.

The checked-in release specification is the authority for the canonical PDF.
An external build deliberately receives new bytes, so it gets a package-local
registry/spec and an audit-only attestation.  The attestation and LaTeX log are
not upload payload and are therefore excluded from ``SHA256SUMS.txt``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path


MAIN_DOCUMENT = "main"
SI_DOCUMENT = "si"
STANDALONE_ITEM = "Graphical abstract"
ASSET_NAMES = ("A3_graphical_abstract.pdf", "A3_graphical_abstract.png")
UPLOAD_PAYLOAD = (
    "manuscript.pdf",
    "SI.pdf",
    "A3_graphical_abstract.pdf",
    "A3_graphical_abstract.png",
    "cover_letter_Proteins.md",
    "submission_checklist.md",
    "MANIFEST.txt",
    "visual-registry.json",
)
MAIN_FIGURE_ITEMS = tuple(f"Fig {n}" for n in range(1, 8))
SI_FIGURE_ITEMS = tuple(f"Fig S{n}" for n in range(1, 10))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalized_build_inputs_digest(inputs: list[dict[str, str]]) -> str:
    material = "".join(
        f"{row['path']}\t{row['sha256'].lower()}\n" for row in sorted(inputs, key=lambda row: row["path"])
    )
    return hashlib.sha256(material.encode("utf-8")).hexdigest()


def page_count(pdf: Path) -> int:
    result = subprocess.run(["pdfinfo", str(pdf)], text=True, capture_output=True, check=True)
    for line in result.stdout.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise ValueError(f"pdfinfo did not report pages for {pdf}")


def extracted_figure_pages(pdf: Path, expected_items: set[str]) -> dict[str, int]:
    """Fail closed unless the final PDF has each authoritative caption once."""
    found: dict[str, int] = {}
    pages = page_count(pdf)
    for page in range(1, pages + 1):
        text = subprocess.run(
            ["pdftotext", "-f", str(page), "-l", str(page), str(pdf), "-"],
            text=True, capture_output=True, check=True,
        ).stdout
        for token in set(re.findall(r"Fig S?[1-9]\.", text)):
            item = token[:-1]
            if item not in expected_items:
                fail(f"final PDF has ghost or unregistered figure caption: {item}")
            if item in found:
                fail(f"final PDF repeats figure caption: {item}")
            found[item] = page
    if set(found) != expected_items:
        fail(f"final PDF figure captions do not match registry: missing={sorted(expected_items-set(found))}")
    return found


def fail(message: str) -> None:
    raise ValueError(message)


def checked_package_path(package: Path, relative: str) -> Path:
    candidate = package / relative
    resolved = candidate.resolve(strict=True)
    if package not in (resolved, *resolved.parents):
        fail(f"package path escapes release root: {relative}")
    if not resolved.is_file():
        fail(f"package path is not a regular file: {relative}")
    return resolved


def verify_attestation(package: Path) -> None:
    """Verify audit evidence without regenerating it; useful for tamper checks."""
    package = package.resolve(strict=True)
    attestation_path = checked_package_path(package, "attestation/build-attestation.json")
    attestation = json.loads(attestation_path.read_text(encoding="utf-8"))
    if attestation.get("capsule_id") != "A3":
        fail("attestation capsule ID is not A3")
    checks = (
        ("package_release_spec", "path"),
        ("document", "path"),
        ("si_document", "path"),
        ("log", "path"),
        ("log_policy", "path"),
    )
    for section, path_key in checks:
        entry = attestation.get(section, {})
        path = checked_package_path(package, entry.get(path_key, ""))
        if digest(path) != entry.get("sha256"):
            fail(f"attestation hash mismatch for {section}")
    package_spec = json.loads(checked_package_path(package, "release-spec.json").read_text(encoding="utf-8"))
    if package_spec.get("build_inputs_digest") != attestation.get("build_inputs_digest"):
        fail("attestation build_inputs_digest does not match package release specification")
    if page_count(checked_package_path(package, attestation["document"]["path"])) != attestation["document"].get("pages"):
        fail("attestation page count does not match manuscript")
    if page_count(checked_package_path(package, attestation["si_document"]["path"])) != attestation["si_document"].get("pages"):
        fail("attestation page count does not match SI")
    payload = attestation.get("upload_payload")
    if set(payload or {}) != set(UPLOAD_PAYLOAD):
        fail("attestation upload payload is incomplete or has unexpected paths")
    for name, expected_hash in payload.items():
        if digest(checked_package_path(package, name)) != expected_hash:
            fail(f"attestation upload payload hash mismatch for {name}")


def bind_and_verify(
    package: Path, source_spec_path: Path, source_registry_path: Path, builder: Path
) -> None:
    package = package.resolve(strict=True)
    source_spec = json.loads(source_spec_path.read_text(encoding="utf-8"))
    source_registry = json.loads(source_registry_path.read_text(encoding="utf-8"))
    if source_spec.get("capsule_id") != "A3":
        fail("source release specification is not for A3")
    build_inputs = source_spec.get("build_inputs")
    if not isinstance(build_inputs, list) or not build_inputs:
        fail("source release specification has no build input closure")
    build_inputs_digest = normalized_build_inputs_digest(build_inputs)
    if source_spec.get("build_inputs_digest") != build_inputs_digest:
        fail("source release specification has an invalid build input digest")

    source_documents = source_spec.get("documents")
    if not isinstance(source_documents, dict) or set(source_documents) != {MAIN_DOCUMENT, SI_DOCUMENT}:
        fail("A3 release specification must define main and si documents")
    source_document = source_documents.get(MAIN_DOCUMENT)
    si_source_document = source_documents.get(SI_DOCUMENT)
    if not isinstance(source_document, dict) or not isinstance(si_source_document, dict):
        fail("A3 release specification has invalid main or si document records")
    source_visual_items = list(source_document.get("visual_items", [])) + list(si_source_document.get("visual_items", []))
    standalone = source_spec.get("standalone_assets", [])
    if not isinstance(standalone, list) or not standalone:
        fail("A3 release specification must declare graphical-abstract assets")
    source_visual_items.append({"id": STANDALONE_ITEM, "document": "standalone"})
    if not isinstance(source_visual_items, list) or not all(isinstance(item, dict) for item in source_visual_items):
        fail("A3 release specification has invalid visual items")

    for name in UPLOAD_PAYLOAD[:-1]:
        checked_package_path(package, name)
    manuscript = checked_package_path(package, "manuscript.pdf")
    si_pdf = checked_package_path(package, "SI.pdf")
    manuscript_hash = digest(manuscript)
    si_hash = digest(si_pdf)
    expected_items = {item["id"] for item in source_visual_items}
    actual_items = [row.get("item") for row in source_registry]
    if not actual_items or len(actual_items) != len(set(actual_items)) or set(actual_items) != expected_items - {STANDALONE_ITEM}:
        fail("source visual registry does not exactly match the non-empty release specification")

    expected_main_items = {row["item"] for row in source_registry if row.get("document") == MAIN_DOCUMENT}
    expected_si_items = {row["item"] for row in source_registry if row.get("document") == SI_DOCUMENT}
    if expected_main_items != set(MAIN_FIGURE_ITEMS):
        fail("source visual registry must declare exactly seven main figures on the main document")
    if expected_si_items != set(SI_FIGURE_ITEMS):
        fail("source visual registry must declare exactly nine supplementary figures on the si document")
    final_pages = extracted_figure_pages(manuscript, expected_main_items)
    si_pages = extracted_figure_pages(si_pdf, expected_si_items)

    standalone_assets: dict[str, str] = {}
    rebound_rows = []
    for row in source_registry:
        row = dict(row)
        if row.get("document") == MAIN_DOCUMENT:
            row["pdf_sha256"] = manuscript_hash
            row["page"] = final_pages[row["item"]]
        elif row.get("document") == SI_DOCUMENT:
            row["pdf_sha256"] = si_hash
            row["page"] = si_pages[row["item"]]
        elif row.get("item") == STANDALONE_ITEM:
            if row.get("document") != "standalone":
                fail("graphical abstract must be a standalone registry item")
            standalone_assets = {name: digest(checked_package_path(package, name)) for name in ASSET_NAMES}
            row["assets"] = standalone_assets
        else:
            fail(f"unexpected registry document for {row.get('item')!r}")
        rebound_rows.append(row)

    if not standalone_assets:
        standalone_assets = {name: digest(checked_package_path(package, name)) for name in ASSET_NAMES}
        rebound_rows.append({"item": STANDALONE_ITEM, "document": "standalone", "page": "standalone", "assets": standalone_assets})

    staged_document = dict(source_document)
    staged_document.update({"id": MAIN_DOCUMENT, "pdf": "manuscript.pdf", "sha256": manuscript_hash, "pages": page_count(manuscript), "log": "attestation/manuscript.log"})
    if staged_document["pages"] != source_document["pages"]:
        fail(
            f"staged manuscript has {staged_document['pages']} pages, "
            f"release specification requires {source_document['pages']}"
        )
    staged_si_document = dict(si_source_document)
    staged_si_document.update({"id": SI_DOCUMENT, "pdf": "SI.pdf", "sha256": si_hash, "pages": page_count(si_pdf), "log": "attestation/SI.log"})
    if staged_si_document["pages"] != si_source_document["pages"]:
        fail(
            f"staged SI has {staged_si_document['pages']} pages, "
            f"release specification requires {si_source_document['pages']}"
        )

    staged_visual_items = []
    for item in source_visual_items:
        item = dict(item)
        if item.get("id") == STANDALONE_ITEM:
            item["assets"] = [{"path": name, "sha256": standalone_assets[name]} for name in ASSET_NAMES]
        staged_visual_items.append(item)

    (package / "visual-registry.json").write_text(json.dumps(rebound_rows, indent=2) + "\n", encoding="utf-8")
    payload_hashes = {name: digest(checked_package_path(package, name)) for name in UPLOAD_PAYLOAD}
    package_spec = {
        "schema_version": 2,
        "capsule_id": "A3",
        "registry": "visual-registry.json",
        "documents": [staged_document, staged_si_document],
        "visual_items": staged_visual_items,
        "standalone_assets": standalone_assets,
        "upload_payload": [{"path": name, "sha256": payload_hashes[name]} for name in UPLOAD_PAYLOAD],
        "build_attestation": "attestation/build-attestation.json",
        "versioned_paths": {"log": "attestation/manuscript.log", "si_log": "attestation/SI.log", "log_policy": "attestation/log-policy.json"},
        "build_inputs_digest": build_inputs_digest,
        "package_checksum": "SHA256SUMS.txt",
    }
    (package / "release-spec.json").write_text(json.dumps(package_spec, indent=2) + "\n", encoding="utf-8")

    log = checked_package_path(package, "attestation/manuscript.log")
    policy = checked_package_path(package, "attestation/log-policy.json")
    attestation = {
        "schema_version": 1,
        "capsule_id": "A3",
        "epoch": "1970-01-01T00:00:00Z",
        "builder": {
            "path": "build_submission_pdf.sh",
            "sha256": digest(builder),
            "command": "bash manuscript/scripts/build_submission_pdf.sh /absolute/nonexistent-output-directory",
        },
        "source_release_spec": {"path": "release-spec.json", "sha256": digest(source_spec_path)},
        "package_release_spec": {"path": "release-spec.json", "sha256": digest(package / "release-spec.json")},
        "document": {"path": "manuscript.pdf", "sha256": manuscript_hash, "pages": staged_document["pages"]},
        "si_document": {"path": "SI.pdf", "sha256": si_hash, "pages": staged_si_document["pages"]},
        "log": {"path": "attestation/manuscript.log", "sha256": digest(log)},
        "log_policy": {"path": "attestation/log-policy.json", "sha256": digest(policy)},
        "upload_payload": payload_hashes,
        "build_inputs_digest": build_inputs_digest,
    }
    attestation_path = package / "attestation" / "build-attestation.json"
    attestation_path.write_text(json.dumps(attestation, indent=2) + "\n", encoding="utf-8")
    verify_attestation(package)

    staged_registry = json.loads((package / "visual-registry.json").read_text(encoding="utf-8"))
    main_rows = [row for row in staged_registry if row.get("document") == MAIN_DOCUMENT]
    if not main_rows or any(row.get("pdf_sha256") != manuscript_hash for row in main_rows):
        fail("staged main registry rows are not bound to staged manuscript.pdf")
    si_rows = [row for row in staged_registry if row.get("document") == SI_DOCUMENT]
    if not si_rows or any(row.get("pdf_sha256") != si_hash for row in si_rows):
        fail("staged SI registry rows are not bound to staged SI.pdf")
    graphical_abstract = next((row for row in staged_registry if row.get("item") == STANDALONE_ITEM), None)
    if graphical_abstract is None or graphical_abstract.get("assets") != standalone_assets:
        fail("staged graphical-abstract registry row is not bound to package assets")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--source-spec", type=Path, required=True)
    parser.add_argument("--source-registry", type=Path, required=True)
    parser.add_argument("--builder", type=Path, required=True)
    args = parser.parse_args()
    try:
        bind_and_verify(args.package, args.source_spec, args.source_registry, args.builder)
    except (OSError, ValueError, subprocess.CalledProcessError, json.JSONDecodeError) as error:
        print(f"A3 staged package validation failed: {error}", file=sys.stderr)
        return 1
    print("A3 staged package registry/spec/attestation binding passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
