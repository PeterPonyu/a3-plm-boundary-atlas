"""Dependency-free Platt and isotonic calibrators (A3-XFIT-CAL-20260918).

Copied from the project Stage-0 calibrators so this capsule analysis does not
import outside the owned tree. Gradient-descent Platt and PAV isotonic are
the same algorithms; this file is not a restoration of any historical fit.
"""
from __future__ import annotations

import math
from bisect import bisect_right
from typing import Callable, Sequence

EPS = 1e-6

Calibrator = Callable[[Sequence[float]], list[float]]


def clip01(p: float) -> float:
    return min(1.0 - EPS, max(EPS, p))


def logit(p: float) -> float:
    p = clip01(p)
    return math.log(p / (1.0 - p))


def sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def fit_platt(
    p_cal: Sequence[float], y_cal: Sequence[float], iters: int = 4000, lr: float = 0.05
) -> Calibrator:
    x = [logit(p) for p in p_cal]
    a, b = 1.0, 0.0
    n = len(x)
    for _ in range(iters):
        ga = gb = 0.0
        for xi, ti in zip(x, y_cal):
            qi = sigmoid(a * xi + b)
            err = qi - ti
            ga += err * xi
            gb += err
        a -= lr * ga / n
        b -= lr * gb / n

    def apply(probs: Sequence[float]) -> list[float]:
        return [sigmoid(a * logit(p) + b) for p in probs]

    apply.coeff = (a, b)  # type: ignore[attr-defined]
    return apply


def fit_isotonic(p_cal: Sequence[float], y_cal: Sequence[float]) -> Calibrator:
    if len(p_cal) != len(y_cal) or not p_cal:
        raise ValueError("isotonic fit requires nonempty equal-length scores and labels")
    pairs = sorted(
        zip((float(p) for p in p_cal), (float(t) for t in y_cal)), key=lambda z: z[0]
    )
    xs: list[float] = []
    aggregated: list[tuple[float, float]] = []
    for p, t in pairs:
        if not math.isfinite(p) or not math.isfinite(t) or t not in (0., 1.):
            raise ValueError("isotonic scores must be finite and labels binary")
        if xs and p == xs[-1]:
            s, c = aggregated[-1]
            aggregated[-1] = (s + t, c + 1.)
        else:
            xs.append(p)
            aggregated.append((t, 1.))
    blocks: list[list[float]] = []
    for s, c in aggregated:
        blocks.append([s, c, s / c, 1.])
        while len(blocks) >= 2 and blocks[-2][2] >= blocks[-1][2]:
            s2, c2, _, unique2 = blocks.pop()
            s1, c1, _, unique1 = blocks.pop()
            s, c = s1 + s2, c1 + c2
            blocks.append([s, c, s / c, unique1 + unique2])
    fitted: list[float] = []
    for _, _, v, unique in blocks:
        # A block may contain several training observations at one x; expand
        # over distinct x values, never over row count.
        fitted.extend([v] * int(unique))

    if len(fitted) != len(xs):
        raise AssertionError("isotonic fit did not assign one value per unique score")

    def apply(probs: Sequence[float]) -> list[float]:
        out: list[float] = []
        for p in probs:
            p = float(p)
            idx = min(max(bisect_right(xs, p) - 1, 0), len(fitted) - 1)
            out.append(clip01(fitted[idx]))
        return out

    return apply


def ece(probs: Sequence[float], labels: Sequence[int], n_bins: int = 10) -> float:
    """Equal-width ECE on [0,1]. Bin b is [b/n, (b+1)/n) except last [0.9, 1]."""
    bins: list[list[tuple[float, int]]] = [[] for _ in range(n_bins)]
    for p, y in zip(probs, labels):
        b = min(n_bins - 1, int(p * n_bins))
        bins[b].append((p, y))
    n = len(probs)
    if n == 0:
        return float("nan")
    total = 0.0
    for block in bins:
        if not block:
            continue
        conf = sum(p for p, _ in block) / len(block)
        acc = sum(y for _, y in block) / len(block)
        total += (len(block) / n) * abs(conf - acc)
    return total


def reliability_points(
    probs: Sequence[float], labels: Sequence[int], n_bins: int = 10
) -> list[dict]:
    bins: list[list[tuple[float, int]]] = [[] for _ in range(n_bins)]
    for p, y in zip(probs, labels):
        b = min(n_bins - 1, int(p * n_bins))
        bins[b].append((p, y))
    out = []
    for bi, block in enumerate(bins):
        if not block:
            continue
        conf = sum(p for p, _ in block) / len(block)
        acc = sum(y for _, y in block) / len(block)
        out.append(
            {
                "bin": bi,
                "bin_lo": bi / n_bins,
                "bin_hi": (bi + 1) / n_bins,
                "bin_mid": (bi + 0.5) / n_bins,
                "mean_pred": conf,
                "obs_rate": acc,
                "n": len(block),
            }
        )
    return out
