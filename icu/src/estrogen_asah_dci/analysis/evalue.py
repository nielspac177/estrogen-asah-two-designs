"""E-values for unmeasured confounding (VanderWeele & Ding 2017).

An E-value is the minimum strength of association (on the risk-ratio scale) that
an unmeasured confounder would need with both exposure and outcome to fully
explain away an observed effect.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


def or_to_rr(odds_ratio: float, baseline_risk: float) -> float:
    """Approximate a risk ratio from an odds ratio given the outcome baseline risk."""
    return odds_ratio / (1.0 - baseline_risk + baseline_risk * odds_ratio)


def _evalue_point(rr: float) -> float:
    r = rr if rr >= 1 else 1.0 / rr
    return r + math.sqrt(r * (r - 1.0))


@dataclass
class EValue:
    point: float
    ci_bound: float | None  # E-value for the CI limit closest to the null (1.0)


def evalue_rr(rr: float, lo: float | None = None, hi: float | None = None) -> EValue:
    """E-value for a risk ratio and (optionally) the CI limit nearest the null."""
    point = _evalue_point(rr)
    ci_bound = None
    if lo is not None and hi is not None:
        if lo <= 1.0 <= hi:
            ci_bound = 1.0  # CI crosses null -> confounder of RR 1 suffices
        else:
            nearest = lo if rr > 1 else hi
            ci_bound = _evalue_point(nearest)
    return EValue(point=point, ci_bound=ci_bound)


def evalue_or(odds_ratio: float, baseline_risk: float,
              lo: float | None = None, hi: float | None = None) -> EValue:
    """E-value for an odds ratio, converting to RR at the given baseline risk."""
    rr = or_to_rr(odds_ratio, baseline_risk)
    rr_lo = or_to_rr(lo, baseline_risk) if lo is not None else None
    rr_hi = or_to_rr(hi, baseline_risk) if hi is not None else None
    return evalue_rr(rr, rr_lo, rr_hi)
