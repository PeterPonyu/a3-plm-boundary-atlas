#!/usr/bin/env python3
"""Refresh the public inventory after reviewed publication edits, using apply_patch.

Original score, analysis and build-source bytes are protected by the export
allowlist. Manuscript text/PDFs and original binding metadata may change during
the coordinator's DOI/DAS integration. Hashes never attest scientific readiness.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = "PUBLIC_RELEASE_MANIFEST.json"
CHECKSUMS = "SHA256SUMS.txt"
STATUS = "AUTHOR_REVIEW_SAME_SCORE_INTEGRATED_NOT_SUBMISSION_READY"
TEMP_DIRS = {".git", "__pycache__", ".pytest_cache"}
TEMP_SUFFIXES = {".pyc", ".pyo", ".aux", ".log", ".out", ".fls", ".fdb_latexmk", ".synctex.gz"}
FORBIDDEN_SUFFIXES = {".pt", ".pth", ".safetensors", ".bin", ".zip", ".tar.gz", ".partial", ".incomplete"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def publication_editable(relative: str) -> bool:
    """Only coordinator publication text/PDF and capsule binding-metadata edits."""
    return (
        relative.startswith("manuscript/") and relative.endswith(".tex")
        or relative in {
            "manuscript/manuscript.pdf", "manuscript/SI.pdf",
            "manuscript/release-current.json", "manuscript/release-spec.json",
            "reproducibility/method-contract.json",
        }
    )


def forbidden(relative: str) -> bool:
    parts = Path(relative).parts
    lower = relative.lower()
    return (
        any(part in {"submission", "parts", ".omc", ".omx", ".superpowers"} for part in parts)
        or any(lower.endswith(suffix) for suffix in FORBIDDEN_SUFFIXES)
        or any(token in lower for token in ("cover_letter", "pre-integration", "pre-unified", "prepolish", "backup"))
    )


def inventory(root: Path) -> list[Path]:
    files = []
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root).as_posix()
        if any(part in TEMP_DIRS for part in path.relative_to(root).parts):
            continue
        if path.is_symlink():
            raise ValueError(f"Symlink is not a portable release input: {relative}")
        if not path.is_file() or any(relative.endswith(suffix) for suffix in TEMP_SUFFIXES):
            continue
        if forbidden(relative):
            raise ValueError(f"Excluded/private/weight/backup payload in public tree: {relative}")
        files.append(path)
    return files


def apply_file(path: Path, text: str) -> None:
    relative = path.relative_to(ROOT).as_posix()
    lines = ["*** Begin Patch"]
    if path.exists():
        original = path.read_text(encoding="utf-8")
        if original == text:
            return
        lines.extend([f"*** Update File: {relative}", "@@"])
        lines.extend("-" + line for line in original.splitlines())
        lines.extend("+" + line for line in text.splitlines())
    else:
        lines.append(f"*** Add File: {relative}")
        lines.extend("+" + line for line in text.splitlines())
    lines.append("*** End Patch")
    subprocess.run(["apply_patch"], input="\n".join(lines) + "\n", text=True, cwd=ROOT, check=True)


def main() -> None:
    allowlist = json.loads((ROOT / "provenance/export-allowlist.json").read_text())
    origin = {item["path"]: item for item in allowlist["source_files"]}
    drift = []
    for relative, record in origin.items():
        path = ROOT / relative
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"Original source dependency missing: {relative}")
        actual = sha(path)
        if actual != record["source_sha256"]:
            if not publication_editable(relative):
                raise ValueError(f"Scientific/result/build source changed; review required: {relative}")
            drift.append({"path": relative, "source_sha256": record["source_sha256"], "public_sha256": actual})
    files = []
    for path in inventory(ROOT):
        relative = path.relative_to(ROOT).as_posix()
        if relative in {MANIFEST, CHECKSUMS}:
            continue
        source = origin.get(relative)
        files.append({
            "path": relative, "sha256": sha(path), "bytes": path.stat().st_size,
            "source_sha256": source["source_sha256"] if source else None,
            "role": source["role"] if source else "public_release_documentation_adapter_test_or_audit",
        })
    source_release = json.loads((ROOT / "manuscript/release-current.json").read_text())
    record = {
        "schema_version": "A3-public-research-release-v1",
        "status": STATUS, "scientific_ready_attested": False,
        "repository": "https://github.com/PeterPonyu/a3-plm-boundary-atlas",
        "zenodo": {"version_doi": "10.5281/zenodo.23101189", "publication_status": "version_bound_release; remote state independently verified in publication receipt", "concept_doi": "10.5281/zenodo.21869994"},
        "source_capsule_manifests": {
            "classification": "source-version provenance; not the current public inventory",
            "method": "reproducibility/method-contract.json",
            "manuscript": "manuscript/release-current.json",
            "all_256_locked_and_237_current_build_input_files_included_at_export": True,
        },
        "verification_mode": "read-only standard-library checksum/row-score/OOF/refit-record reconciliation; no new inference or fits",
        "reference_sequence_boundary": allowlist["reference_sequence_payload"],
        "source_publication_edits": drift,
        "source_document_record": source_release["documents"],
        "files": files,
        "file_count_excluding_manifest_and_checksums": len(files),
        "bytes_excluding_manifest_and_checksums": sum(item["bytes"] for item in files),
        "self_hash_policy": "Manifest excludes itself and SHA256SUMS.txt. SHA256SUMS.txt includes this manifest and every listed file, not itself.",
    }
    apply_file(ROOT / MANIFEST, json.dumps(record, indent=2, sort_keys=True) + "\n")
    sums = "".join(f"{item['sha256']}  {item['path']}\n" for item in files)
    sums += f"{sha(ROOT / MANIFEST)}  {MANIFEST}\n"
    apply_file(ROOT / CHECKSUMS, sums)
    print(json.dumps({"files": len(files) + 2, "payload_bytes": record["bytes_excluding_manifest_and_checksums"], "publication_source_edits": len(drift), "scientific_ready_attested": False}))


if __name__ == "__main__":
    main()
