#!/usr/bin/env python3
"""Robustness of the locked N4 antibody_interface cell (n=525, rho=0.165).

650M per-mutation scored.csv is not on disk (n4_stats.py expected
~/biochem-data/N4_skempi/scored.csv). This script therefore cannot
re-bootstrap the locked 650M pairs. It:

  1. restates the locked published cell from results.json
  2. recomputes Spearman + bootstrap/jackknife on the surviving 3B
     WT-marginal scores for the same ab ∩ interface mutation set
  3. reports leave-one-PDB-out and 1%/99% ddG-winsor sensitivity

This is a sensitivity analysis, not a replacement of the locked 0.165 cell.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
OUT = HERE / "outputs"
CAP = HERE.parents[1]
LOCKED = CAP / "results-snapshot" / "N4_antibody_ddg" / "results.json"
SCORED_3B = INPUTS / "scored_3b_wt.csv"
SEED = 20260918
N_BOOT = 2000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


def pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    mx = sum(x) / n
    my = sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    dx = math.sqrt(sum((a - mx) ** 2 for a in x))
    dy = math.sqrt(sum((b - my) ** 2 for b in y))
    if dx == 0.0 or dy == 0.0:
        return float("nan")
    return num / (dx * dy)


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(ranks(x), ranks(y))


def percentile(values: list[float], q: float) -> float:
    ordered = [v for v in values if not math.isnan(v)]
    ordered.sort()
    if not ordered:
        return float("nan")
    pos = (len(ordered) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(ordered) - 1)
    w = pos - lo
    return ordered[lo] * (1.0 - w) + ordered[hi] * w


def boot_spearman(x: list[float], y: list[float], seed: int, n_boot: int) -> dict:
    rng = random.Random(seed)
    n = len(x)
    point = spearman(x, y)
    draws = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        draws.append(spearman([x[i] for i in idx], [y[i] for i in idx]))
    return {
        "rho": point,
        "ci_lo": percentile(draws, 0.025),
        "ci_hi": percentile(draws, 0.975),
        "n": n,
        "n_boot": n_boot,
        "seed": seed,
    }


def winsor(values: list[float], lo_q: float, hi_q: float) -> list[float]:
    ordered = sorted(values)
    lo = percentile(ordered, lo_q)
    hi = percentile(ordered, hi_q)
    return [min(hi, max(lo, v)) for v in values]


def load_3b() -> list[dict]:
    rows = []
    with SCORED_3B.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["llr_3b_wt"] == "" or row["ddg"] == "":
                continue
            rows.append(
                {
                    "mutkey": row["mutkey"],
                    "pdb": row["pdb"],
                    "pdbid": row["pdbid"],
                    "llr": float(row["llr_3b_wt"]),
                    "ddg": float(row["ddg"]),
                    "ab": row["ab"] == "True",
                    "interface": row["interface"] == "True",
                }
            )
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    locked = json.loads(LOCKED.read_text(encoding="utf-8"))
    cell = locked["antibody_interface"]
    rows = load_3b()
    abint = [r for r in rows if r["ab"] and r["interface"]]
    x = [r["llr"] for r in abint]
    y = [-r["ddg"] for r in abint]
    primary = boot_spearman(x, y, seed=SEED, n_boot=N_BOOT)
    original_seed0 = boot_spearman(x, y, seed=0, n_boot=1000)
    y_win = [-d for d in winsor([r["ddg"] for r in abint], 0.01, 0.99)]
    winsorized = boot_spearman(x, y_win, seed=SEED, n_boot=N_BOOT)

    by_pdb: dict[str, list[int]] = defaultdict(list)
    for i, row in enumerate(abint):
        by_pdb[row["pdbid"]].append(i)
    jack = []
    for pdb, idx_drop in sorted(by_pdb.items()):
        keep = [i for i in range(len(abint)) if i not in set(idx_drop)]
        if len(keep) < 10:
            continue
        rho = spearman([x[i] for i in keep], [y[i] for i in keep])
        jack.append(
            {
                "pdbid": pdb,
                "n_dropped": len(idx_drop),
                "n_kept": len(keep),
                "rho": rho,
                "delta_vs_full": rho - primary["rho"],
            }
        )
    jack.sort(key=lambda r: abs(r["delta_vs_full"]), reverse=True)

    payload = {
        "analysis": "A3-ABINT-ROBUST-20260918",
        "locked_650M_cell": {
            "source": str(LOCKED.relative_to(CAP)),
            "source_sha256": sha256(LOCKED),
            "n": cell["n"],
            "rho": cell["rho"],
            "rounded_3dp": round(cell["rho"], 3),
            "ci_lo": cell["ci_lo"],
            "ci_hi": cell["ci_hi"],
            "printed": "0.165 [0.075, 0.249]",
            "historical_bootstrap": "n4_stats.py seed=0 B=1000 mutation-level percentile",
        },
        "650M_row_status": {
            "scored_csv_expected": "~/biochem-data/N4_skempi/scored.csv",
            "present": False,
            "consequence": "cannot independently re-bootstrap the locked 650M n=525 pairs from raw LLR",
        },
        "sensitivity_3B": {
            "file": "scored_3b_wt.csv",
            "sha256": sha256(SCORED_3B),
            "n_mapped": len(rows),
            "n_antibody_interface": len(abint),
            "matches_locked_n": len(abint) == 525,
            "model": "ESM-2 3B WT-marginal (not the locked 650M score)",
            "bootstrap_seed_20260918_B2000": primary,
            "bootstrap_seed_0_B1000_like_n4_stats": original_seed0,
            "ddg_winsor_1_99": winsorized,
            "n_pdbs": len(by_pdb),
            "leave_one_pdb_max_abs_delta": jack[0] if jack else None,
            "leave_one_pdb_rho_min": min(r["rho"] for r in jack) if jack else None,
            "leave_one_pdb_rho_max": max(r["rho"] for r in jack) if jack else None,
        },
    }
    (OUT / "antibody_interface_robustness.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    with (OUT / "antibody_interface_loo_pdb.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["pdbid", "n_dropped", "n_kept", "rho", "delta_vs_full"],
        )
        writer.writeheader()
        for row in jack:
            writer.writerow({k: (round(v, 6) if isinstance(v, float) else v) for k, v in row.items()})
    print("=== A3-ABINT-ROBUST-20260918 ===")
    print(
        f"locked 650M n={cell['n']} rho={cell['rho']:.6f} "
        f"[{cell['ci_lo']:.6f}, {cell['ci_hi']:.6f}] -> printed 0.165 [0.075, 0.249]"
    )
    print("650M scored.csv ABSENT; cannot re-bootstrap locked pairs")
    print(
        f"3B ab∩interface n={len(abint)} rho={primary['rho']:.4f} "
        f"[{primary['ci_lo']:.4f}, {primary['ci_hi']:.4f}] seed={SEED} B={N_BOOT}"
    )
    if jack:
        print(
            f"LOO-PDB n_pdb={len(by_pdb)} rho range "
            f"[{payload['sensitivity_3B']['leave_one_pdb_rho_min']:.4f}, "
            f"{payload['sensitivity_3B']['leave_one_pdb_rho_max']:.4f}] "
            f"worst {jack[0]['pdbid']} delta={jack[0]['delta_vs_full']:.4f}"
        )


if __name__ == "__main__":
    main()
