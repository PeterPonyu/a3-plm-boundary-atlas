#!/usr/bin/env python3
"""Portable explicit-input adapter for the unchanged completed A3 scorer.

Not called by public verification. --preflight-only reads identities without
model inference; the full route requires official weights and upstream inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_scorer():
    path = ROOT / "analysis/revision_20260922/rescore_clinvar_candidate.py"
    spec = importlib.util.spec_from_file_location("a3_original_completed_scorer", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load exported scorer")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical-csvs", required=True, type=Path)
    parser.add_argument("--n6c-summary", required=True, type=Path)
    parser.add_argument("--n6c-scorer", required=True, type=Path)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    if not args.preflight_only and (args.checkpoint is None or args.out is None):
        parser.error("Full scoring requires explicit --checkpoint and --out")
    if args.out is not None and args.out.resolve().is_relative_to(ROOT):
        parser.error("New scoring output must be outside the published release tree")
    provenance = json.loads((ROOT / "analysis/revision_20260923/row-scores/provenance.json").read_text())
    binding = json.loads((ROOT / "analysis/revision_20260923/row-scores/run-binding.json").read_text())
    scorer = load_scorer()
    scorer.CLINICAL = args.clinical_csvs.resolve(strict=True)
    scorer.SUMMARY = args.n6c_summary.resolve(strict=True)
    scorer.SCORER = args.n6c_scorer.resolve(strict=True)
    if sha(scorer.SUMMARY) != provenance["N6c_summary_sha256"] or sha(scorer.SCORER) != provenance["N6c_scorer_sha256"]:
        raise ValueError("Upstream N6c summary/scorer bytes do not match recorded run")
    if sha(Path(scorer.__file__)) != binding["scorer_sha256"] or sha(scorer.MANIFEST) != binding["manifest_sha256"]:
        raise ValueError("Exported scientific scorer/manifest identity changed")
    grouped, references, source_hashes, _ = scorer.source_rows()
    if references != provenance["wt_sequence_sha256"] or source_hashes != provenance["clinical_source_sha256"]:
        raise ValueError("Upstream clinical labels/reference/source bytes differ from completed run")
    print(json.dumps({"bound_variants": sum(map(len, grouped.values())), "bound_genes": len(grouped), "source_identity_matches": True, "model_inference_performed": False}, indent=2))
    if args.preflight_only:
        return
    forwarded = [sys.argv[0], "--checkpoint", str(args.checkpoint), "--out", str(args.out), "--device", args.device]
    if args.resume:
        forwarded.append("--resume")
    sys.argv = forwarded
    scorer.main()


if __name__ == "__main__":
    main()
