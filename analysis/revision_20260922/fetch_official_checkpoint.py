#!/usr/bin/env python3
"""Resumable, byte-range fetch of the official fair-esm ESM-2 650M checkpoint.

The completed file lives in a new candidate directory; the truncated torch
cache is never overwritten. Each part validates HTTP 206 and Content-Range,
then the assembled file validates byte length before use. The clinical scorer
still verifies model loading and the historical N6c per-gene AUROC separately.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "candidate-v1-checkpoint"
URL = "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt"
LENGTH = 2604537549
PART = 16 * 1024 * 1024
N = (LENGTH + PART - 1) // PART


def expected(i):
    start = PART * i
    end = min(LENGTH, start+PART) - 1
    return start, end, end-start+1


def download(i, retries):
    folder = OUT / "parts"
    start, end, n = expected(i)
    complete = folder / f"{i:04d}.bin"
    if complete.exists() and complete.stat().st_size == n:
        return i
    partial = folder / f"{i:04d}.partial"
    for attempt in range(retries):
        offset = partial.stat().st_size if partial.exists() else 0
        if offset > n:
            raise ValueError(f"Part {i} exceeds its expected size")
        if offset == n:
            partial.rename(complete)
            return i
        request = urllib.request.Request(URL, headers={
            "Range": f"bytes={start+offset}-{end}",
            "User-Agent": "A3-N6c-resumable-checkpoint/1",
        })
        try:
            with urllib.request.urlopen(request, timeout=45) as response:
                if response.status != 206:
                    raise ValueError(f"Part {i}: expected HTTP 206, got {response.status}")
                expected_range = f"bytes {start+offset}-{end}/{LENGTH}"
                if response.headers.get("Content-Range") != expected_range:
                    raise ValueError(f"Part {i}: range mismatch {response.headers.get('Content-Range')}")
                with partial.open("ab") as handle:
                    while chunk := response.read(1 << 20):
                        handle.write(chunk)
            if partial.stat().st_size == n:
                partial.rename(complete)
                return i
        except (OSError, urllib.error.URLError) as exc:
            if attempt == retries-1:
                raise RuntimeError(f"Part {i} incomplete after {retries} attempts") from exc
            time.sleep(min(3+attempt, 15))
    raise RuntimeError(f"Part {i} incomplete after {retries} attempts")


def assemble():
    target = OUT / "esm2_t33_650M_UR50D.pt"
    if target.exists():
        raise SystemExit(f"Refusing to overwrite an existing checkpoint: {target}")
    temporary = OUT / "esm2_t33_650M_UR50D.incomplete"
    if temporary.exists():
        raise SystemExit(f"Refusing to overwrite partial assembly: {temporary}")
    sha = hashlib.sha256()
    with temporary.open("xb") as destination:
        for i in range(N):
            part = OUT / "parts" / f"{i:04d}.bin"
            if part.stat().st_size != expected(i)[2]:
                raise RuntimeError(f"Invalid part size {i}")
            with part.open("rb") as handle:
                while block := handle.read(1 << 20):
                    destination.write(block)
                    sha.update(block)
        destination.flush()
        os.fsync(destination.fileno())
    if temporary.stat().st_size != LENGTH:
        raise AssertionError("Assembled checkpoint has incorrect byte length")
    temporary.rename(target)
    (OUT / "download-provenance.json").write_text(json.dumps(dict(
        official_url=URL, bytes=LENGTH, sha256=sha.hexdigest(),
        range_size=PART, n_parts=N, clinical_score_status="not yet inferred"), indent=2) + "\n")
    print(f"Complete candidate checkpoint {target}, SHA256 {sha.hexdigest()}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--retries", type=int, default=12)
    ap.add_argument("--status", action="store_true")
    args = ap.parse_args()
    if args.workers < 1 or args.retries < 1:
        raise SystemExit("Positive workers and retries required")
    done = [i for i in range(N) if (OUT/"parts"/f"{i:04d}.bin").exists() and
            (OUT/"parts"/f"{i:04d}.bin").stat().st_size == expected(i)[2]]
    print(f"Official ESM-2 checkpoint: {len(done)}/{N} complete ranges", flush=True)
    if args.status:
        return
    if (OUT / "esm2_t33_650M_UR50D.pt").exists():
        raise SystemExit("Candidate checkpoint already exists; validate it before inference")
    (OUT/"parts").mkdir(parents=True, exist_ok=True)
    failures = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        tasks = [pool.submit(download, i, args.retries) for i in range(N) if i not in done]
        for j, task in enumerate(as_completed(tasks), 1):
            try:
                task.result()
            except Exception as exc:
                failures.append(str(exc))
                print(f"Range failed, retaining partial bytes: {exc}", flush=True)
            if j % 8 == 0 or j == len(tasks):
                print(f"Ranges complete: {len(done)+j-len(failures)}/{N}", flush=True)
    if failures:
        raise SystemExit(f"{len(failures)} ranges still incomplete; rerun to resume")
    assemble()


if __name__ == "__main__":
    main()
