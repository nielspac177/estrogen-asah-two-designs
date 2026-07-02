"""Fixed- and random-effects pooling of per-stratum estimates (ADR-0002).

Estimates are combined on the log-OR scale by inverse-variance weighting;
random effects use the DerSimonian-Laird tau^2.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class Pooled:
    estimate: float          # OR
    ci_low: float
    ci_high: float
    tau2: float
    q: float
    k: int
    model: str


def _combine(log_ests: list[float], ses: list[float], random: bool) -> Pooled:
    k = len(log_ests)
    w = [1.0 / (se * se) for se in ses]
    fixed = sum(wi * li for wi, li in zip(w, log_ests, strict=True)) / sum(w)
    q = sum(wi * (li - fixed) ** 2 for wi, li in zip(w, log_ests, strict=True))
    tau2 = 0.0
    if random and k > 1:
        c = sum(w) - sum(wi * wi for wi in w) / sum(w)
        tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
        w = [1.0 / (se * se + tau2) for se in ses]
    pooled = sum(wi * li for wi, li in zip(w, log_ests, strict=True)) / sum(w)
    se_pooled = math.sqrt(1.0 / sum(w))
    return Pooled(
        estimate=math.exp(pooled),
        ci_low=math.exp(pooled - 1.96 * se_pooled),
        ci_high=math.exp(pooled + 1.96 * se_pooled),
        tau2=tau2,
        q=q,
        k=k,
        model="random" if random else "fixed",
    )


def pool_or(estimates: list[dict], random: bool = False) -> Pooled:
    """Pool a list of dicts with keys 'or', 'ci_low', 'ci_high' (per stratum)."""
    log_ests, ses = [], []
    for e in estimates:
        log_or = math.log(e["or"])
        se = (math.log(e["ci_high"]) - math.log(e["ci_low"])) / (2 * 1.96)
        log_ests.append(log_or)
        ses.append(se)
    return _combine(log_ests, ses, random)
