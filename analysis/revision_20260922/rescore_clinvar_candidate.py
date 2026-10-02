#!/usr/bin/env python3
"""Recover row-level N6c fp32 WT-marginal LLRs, without touching old scores.

The source clinical CSVs and 20-gene manifest must match variant by variant.
The official fair-esm checkpoint must be complete; a partial download is never
loaded. A per-gene rounded AUROC gate compares the new predictions with the
historic N6c summary. Any difference is an explicit failure, not an adjusted
score. This script neither trains a model nor changes the submitted paper.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
import zipfile
from collections import Counter, defaultdict
from importlib.metadata import version
from pathlib import Path

HERE = Path(__file__).resolve().parent
CAPSULE = HERE.parents[1]
SOURCE = CAPSULE.parents[1] / "active/reliability-program/biochemical/biochem-research/research"
N6 = SOURCE / "results/N6_clinvar"
SCORER = SOURCE / "scaffolds/protein_ligand/stage0/featurize/n6c_consistent_rescore.py"
CLINICAL = N6 / "gof_test/clinical_csvs"
MANIFEST = CAPSULE / "analysis/xfit_calibration_20260918/inputs/manifest_gof_650M.csv"
SUMMARY = N6 / "n6c_consistent.json"
DEFAULT_CHECKPOINT = Path.home() / ".cache/torch/hub/checkpoints/esm2_t33_650M_UR50D.pt"
CHECKPOINT_BYTES = 2604537549  # official checkpoint HTTP Content-Length, verify before loading
OUT = HERE / "candidate-v1-row-scores"
CTX = 1024
AA = set("ACDEFGHIKLMNPQRSTVWY")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def require_complete_checkpoint(path: Path) -> int:
    size = path.stat().st_size if path.exists() else 0
    if size != CHECKPOINT_BYTES:
        raise ValueError(f"Incomplete checkpoint: {size}/{CHECKPOINT_BYTES}; no inference attempted")
    return size


def source_rows():
    manifest = list(csv.DictReader(MANIFEST.open(newline="")))
    if len(manifest) != 4677 or len({r["example_id"] for r in manifest}) != len(manifest):
        raise ValueError("Expected 4,677 unique manifest examples")
    grouped = defaultdict(list)
    for row in manifest:
        grouped[row["target_id"]].append(row)
    n6 = json.loads(SUMMARY.read_text())
    if set(grouped) != set(n6["per_gene"]) or len(grouped) != 20:
        raise ValueError("Not the same 20 N6c genes")
    sequence_hashes, clinical_hashes = {}, {}
    for gene, entries in grouped.items():
        path = CLINICAL / f"{gene}.csv"
        clinical_hashes[gene] = digest(path)
        clinical = list(csv.DictReader(path.open(newline="")))
        by_mut = defaultdict(list)
        wildtypes = set()
        for r in clinical:
            mut = re.fullmatch(r"([A-Z])(\d+)([A-Z])", r["mutant"])
            if mut is None:
                raise ValueError(f"Invalid source HGVS amino-acid mutation: {gene} {r['mutant']}")
            wt, pos, alt = mut.group(1), int(mut.group(2)), mut.group(3)
            seq = r["mutated_sequence"]
            if alt not in AA or not 1 <= pos <= len(seq) or seq[pos-1] != alt:
                raise ValueError(f"Invalid source sequence: {gene} {r['mutant']}")
            wildtypes.add(seq[:pos-1] + wt + seq[pos:])
            by_mut[r["mutant"]].append(r)
        if len(wildtypes) != 1:
            raise ValueError(f"Multiple source WT sequences: {gene}")
        wt_seq = next(iter(wildtypes))
        if len(wt_seq) != n6["per_gene"][gene]["wt_len"]:
            raise ValueError(f"Wrong N6c WT length: {gene}")
        sequence_hashes[gene] = hashlib.sha256(wt_seq.encode()).hexdigest()
        if Counter(r["ligand_id"] for r in entries) != Counter(r["mutant"] for r in clinical):
            raise ValueError(f"Manifest/source variant multiplicity mismatch: {gene}")
        for m in entries:
            origin = by_mut[m["ligand_id"]]
            if len(origin) != 1:
                raise ValueError(f"Nonunique clinical source mutation: {gene} {m['ligand_id']}")
            if int(float(m["label_value"])) != int(float(origin[0]["DMS_score_bin"])):
                raise ValueError(f"ClinVar label mismatch: {gene} {m['ligand_id']}")
            mut = re.fullmatch(r"([A-Z])(\d+)([A-Z])", m["ligand_id"])
            if not mut or wt_seq[int(mut.group(2))-1] != mut.group(1):
                raise ValueError(f"Reference differs at variant: {gene} {m['ligand_id']}")
    return grouped, sequence_hashes, clinical_hashes, n6


def score_gene(gene, entries, wt_seq, model, alphabet, bc, device, diagnostics=None):
    import torch

    idx = {a: alphabet.get_idx(a) for a in AA}

    @torch.no_grad()
    def logp(seq):
        _, _, token = bc([("wt", seq)])
        result = model(token.to(device))["logits"][0].float()
        return torch.log_softmax(result, dim=-1).cpu()

    # N6c computes whole-sequence fp32 for every gene. For >1024, it also
    # computes a position-centered <=1024-residue window as OOD sensitivity.
    whole = logp(wt_seq)
    if diagnostics is not None:
        repeat = logp(wt_seq)
        diagnostics[gene] = {"repeat_max_abs_logp_diff": float((whole-repeat).abs().max())}
        if not torch.equal(whole, repeat):
            raise ValueError(f"Same-input scoring is not bitwise repeatable: {gene}")
        del repeat
    if len(wt_seq) <= CTX:
        windows = {}
    else:
        half = CTX // 2 - 1
        windows = {}
        for pos in sorted({int(re.fullmatch(r"([A-Z])(\d+)([A-Z])", m["ligand_id"]).group(2))
                           for m in entries}):
            start = max(0, pos-1-half)
            end = min(len(wt_seq), pos-1+half+1)
            if end-start > CTX:
                raise AssertionError("Window longer than N6c context")
            windows[pos] = (logp(wt_seq[start:end])[pos-start], start, end)
    output = []
    for m in entries:
        mut = re.fullmatch(r"([A-Z])(\d+)([A-Z])", m["ligand_id"])
        w, pos, a = mut.group(1), int(mut.group(2)), mut.group(3)
        whole_llr = float(whole[pos, idx[a]] - whole[pos, idx[w]])
        if windows:
            token, start, end = windows[pos]
            win_llr = float(token[idx[a]]-token[idx[w]])
            selected, mode = win_llr, "windowed"
        else:
            start, end, win_llr = 0, len(wt_seq), None
            selected, mode = whole_llr, "whole"
        output.append(dict(gene=gene, example_id=m["example_id"], mutation=m["ligand_id"],
                           source_file=f"{gene}.csv", label_benign=int(float(m["label_value"])),
                           wt_sha256=hashlib.sha256(wt_seq.encode()).hexdigest(), wt_len=len(wt_seq),
                           window_start_zero_based=start, window_end_exclusive=end,
                           llr_whole_fp32=whole_llr, llr_window_fp32=win_llr,
                           primary_mode=mode, primary_llr_fp32=selected))
    return output


def source_wt(gene: str, wt_hash: str) -> str:
    """Reconstruct from a validated source record and recheck its identity."""
    clinical = list(csv.DictReader((CLINICAL/f"{gene}.csv").open(newline="")))
    mut = re.fullmatch(r"([A-Z])(\d+)([A-Z])", clinical[0]["mutant"])
    pos = int(mut.group(2))
    seq = clinical[0]["mutated_sequence"]
    wt_seq = seq[:pos-1] + mut.group(1) + seq[pos:]
    if hashlib.sha256(wt_seq.encode()).hexdigest() != wt_hash:
        raise ValueError(f"Source sequence changed during run: {gene}")
    return wt_seq


def check_auroc(rows, n6):
    from sklearn.metrics import roc_auc_score

    observed = {}
    for gene, rr in rows.items():
        expected = n6["per_gene"][gene]
        ys = [r["label_benign"] for r in rr]
        if sum(ys) != expected["n_clinvar_benign"] or len(ys)-sum(ys) != expected["n_path"]:
            raise ValueError(f"N6c source counts disagree: {gene}")
        modes = ("whole", "win") if expected["windowed"] else ("whole",)
        for mode in modes:
            scores = [r["llr_whole_fp32" if mode == "whole" else "llr_window_fp32"] for r in rr]
            got = float(roc_auc_score(ys, scores))
            reference = expected.get(f"A_clinvar_{mode}")
            if reference is not None and abs(got-reference) > 0.000051:
                raise ValueError(f"N6c AUROC mismatch {gene} {mode}: {got:.7f} vs {reference:.4f}")
            observed[gene + ":" + mode] = dict(auroc=got, N6c_rounded=reference)
    return observed


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    ap.add_argument("--preflight-only", action="store_true")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()
    grouped, sequence_hashes, clinical_hashes, n6 = source_rows()
    checkpoint_size = args.checkpoint.stat().st_size if args.checkpoint.exists() else 0
    print(json.dumps(dict(n_variants=sum(map(len, grouped.values())), n_genes=len(grouped),
                          checkpoint_bytes=checkpoint_size, expected_checkpoint_bytes=CHECKPOINT_BYTES,
                          source_manifest_sha256=digest(MANIFEST),
                          reference_wt_sha256=sequence_hashes), indent=2), flush=True)
    if args.preflight_only:
        return
    out = args.out
    if out.exists() and not args.resume:
        raise SystemExit(f"Refusing to overwrite candidate output: {out}; use --resume with bound inputs")
    require_complete_checkpoint(args.checkpoint)
    download_provenance = args.checkpoint.parent / "download-provenance.json"
    checkpoint_hash = digest(args.checkpoint)
    if not download_provenance.is_file():
        raise SystemExit("Official download provenance is required for this recovery run")
    downloaded = json.loads(download_provenance.read_text())
    if downloaded["sha256"] != checkpoint_hash or downloaded["bytes"] != CHECKPOINT_BYTES:
        raise SystemExit("Checkpoint does not match its completed-download record")
    with zipfile.ZipFile(args.checkpoint) as archive:
        bad_member = archive.testzip()
        if bad_member is not None:
            raise SystemExit(f"Checkpoint ZIP CRC failed: {bad_member}")
    print("Checkpoint length, SHA-256 binding and all ZIP CRCs verified", flush=True)

    binding = dict(checkpoint_sha256=checkpoint_hash,
                   manifest_sha256=digest(MANIFEST), N6c_summary_sha256=digest(SUMMARY),
                   clinical_source_sha256=clinical_hashes,
                   scorer_sha256=digest(Path(__file__)), context=CTX, device=args.device)
    out.mkdir(parents=True, exist_ok=True)
    binding_file = out / "run-binding.json"
    if binding_file.exists():
        if json.loads(binding_file.read_text()) != binding:
            raise SystemExit("Resume binding changed; use a new output directory")
    else:
        binding_file.write_text(json.dumps(binding, indent=2) + "\n")

    import torch
    import esm
    if args.device == "cuda" and not torch.cuda.is_available():
        raise SystemExit("CUDA requested but not available")
    torch.set_num_threads(4)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    model_data = torch.load(str(args.checkpoint), map_location="cpu", weights_only=False)
    model, alphabet = esm.pretrained.load_model_and_alphabet_core(args.checkpoint.stem, model_data)
    del model_data
    model = model.float().to(args.device).eval()
    print(f"Loaded {sum(p.numel() for p in model.parameters()):,} parameters in fp32", flush=True)
    bc = alphabet.get_batch_converter()
    all_rows = {}
    diagnostics = {}
    for gene, entries in sorted(grouped.items()):
        cache = out / f"{gene}.json"
        if cache.exists():
            saved = json.loads(cache.read_text())
            if {r["example_id"] for r in saved["rows"]} != {r["example_id"] for r in entries}:
                raise ValueError(f"Resumed identity mismatch: {gene}")
            all_rows[gene] = saved["rows"]
            diagnostics[gene] = saved["diagnostics"]
            print(f"reused bound scores {gene}: {len(entries)} clinical variants", flush=True)
            continue
        wt_seq = source_wt(gene, sequence_hashes[gene])
        started = time.monotonic()
        all_rows[gene] = score_gene(gene, entries, wt_seq, model, alphabet, bc, args.device, diagnostics)
        diagnostics[gene]["seconds"] = time.monotonic() - started
        cache.write_text(json.dumps(dict(rows=all_rows[gene], diagnostics=diagnostics[gene]), indent=2) + "\n")
        print(f"scored {gene}: {len(entries)} clinical variants, {diagnostics[gene]['seconds']:.1f}s", flush=True)
    try:
        gates = check_auroc(all_rows, n6)
    except ValueError as exc:
        (out / "AUROC-GATE-FAILED.json").write_text(json.dumps(dict(
            error=str(exc), status="Rows preserved for diagnosis; not manuscript evidence",
            diagnostics=diagnostics), indent=2) + "\n")
        raise
    rows = [row for group in all_rows.values() for row in group]
    if len(rows) != 4677 or len({r["example_id"] for r in rows}) != 4677:
        raise AssertionError("Incomplete scored inventory")
    with (out / "scores.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (out / "provenance.json").write_text(json.dumps(dict(
        status="internal candidate; manuscript HOLD", scorer="N6c fp32 WT-marginal",
        checkpoint_sha256=checkpoint_hash, checkpoint_bytes=checkpoint_size,
        scores_sha256=digest(out / "scores.csv"), checkpoint_zip_crc="all members passed",
        clinical_source_sha256=clinical_hashes, wt_sequence_sha256=sequence_hashes,
        manifest_sha256=digest(MANIFEST), N6c_summary_sha256=digest(SUMMARY),
        N6c_scorer_sha256=digest(SCORER), versions={name: version(name) for name in
                                           ("torch", "fair-esm", "scikit-learn")},
        context=CTX, device=args.device, n_rows=len(rows), auc_gates=gates,
        diagnostics=diagnostics, tf32=False, runtime_gpu=torch.cuda.get_device_name(0)
        if args.device == "cuda" else None), indent=2) + "\n")
    print("saved", out, flush=True)


if __name__ == "__main__":
    main()
