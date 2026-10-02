#!/usr/bin/env python3
"""Stage the exact, hash-bound inputs required by an A3 release build."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalized_digest(inputs: list[dict[str, str]]) -> str:
    """Digest a canonical path/hash sequence, independent of JSON formatting."""
    material = "".join(f"{row['path']}\t{row['sha256'].lower()}\n" for row in sorted(inputs, key=lambda row: row["path"]))
    return hashlib.sha256(material.encode("utf-8")).hexdigest()


def checked_inputs(root: Path, spec_path: Path) -> list[dict[str, str]]:
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    inputs = spec.get("build_inputs")
    if not isinstance(inputs, list) or not inputs:
        raise ValueError("release specification has no build_inputs closure")
    if spec.get("build_inputs_digest") != normalized_digest(inputs):
        raise ValueError("release specification build_inputs_digest is not normalized")
    versioned = set(spec.get("versioned_paths", []))
    seen: set[str] = set()
    checked: list[dict[str, str]] = []
    for row in inputs:
        if not isinstance(row, dict) or set(row) != {"path", "sha256"}:
            raise ValueError("each build input must have exactly path and sha256")
        relative, expected = row["path"], row["sha256"]
        candidate = Path(relative)
        if (not isinstance(relative, str) or not isinstance(expected, str) or len(expected) != 64
                or candidate.is_absolute() or ".." in candidate.parts or relative in seen):
            raise ValueError(f"invalid build input: {relative!r}")
        seen.add(relative)
        if relative not in versioned:
            raise ValueError(f"build input omitted from versioned_paths: {relative}")
        source = root / candidate
        if source.is_symlink() or not source.is_file() or not source.resolve().is_relative_to(root.resolve()):
            raise ValueError(f"build input is not a regular file inside source root: {relative}")
        if digest(source) != expected:
            raise ValueError(f"build input hash mismatch: {relative}")
        checked.append({"path": relative, "sha256": expected})
    return checked


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--stage", type=Path, required=True)
    args = parser.parse_args()
    root = args.source_root.resolve(strict=True)
    spec = args.spec.resolve(strict=True)
    stage = args.stage
    if stage.exists() or stage.is_symlink():
        raise SystemExit("A3 release input stage must not already exist")
    try:
        inputs = checked_inputs(root, spec)
        stage.mkdir(parents=True)
        # The authority itself cannot self-hash into its own closure.  It is
        # copied only after its hash-bound input manifest has been checked.
        shutil.copyfile(spec, stage / "release-spec.json")
        for row in inputs:
            destination = stage / row["path"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(root / row["path"], destination)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        shutil.rmtree(stage, ignore_errors=True)
        raise SystemExit(f"A3 release input staging failed: {error}")
    print(f"A3 staged {len(inputs)} hash-bound release inputs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
