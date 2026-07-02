"""Two-sample summary-data Mendelian-randomization estimators.

Inputs are harmonized per-SNP effect sizes on a common effect allele:
``beta_exp, se_exp`` (SNP->exposure) and ``beta_out, se_out`` (SNP->outcome).
Implemented from first principles for transparency: inverse-variance weighted
(IVW), MR-Egger (with a directional-pleiotropy intercept test), and the weighted
median. Estimates are on the log-odds scale when the outcome GWAS is binary, so
``exp(estimate)`` is the causal odds ratio per unit exposure.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass
class MRResult:
    method: str
    estimate: float          # causal effect (log-OR per unit exposure)
    se: float
    ci_low: float
    ci_high: float
    p: float
    n_snps: int

    @property
    def odds_ratio(self) -> float:
        return float(np.exp(self.estimate))

    def as_dict(self) -> dict:
        return {
            "method": self.method, "estimate": self.estimate, "se": self.se,
            "ci_low": self.ci_low, "ci_high": self.ci_high, "p": self.p,
            "or": self.odds_ratio, "or_ci_low": float(np.exp(self.ci_low)),
            "or_ci_high": float(np.exp(self.ci_high)), "n_snps": self.n_snps,
        }


def _ci(est: float, se: float, method: str, n: int, extra_dof: int = 0) -> MRResult:
    z = est / se if se > 0 else np.nan
    p = float(2 * stats.norm.sf(abs(z))) if se > 0 else np.nan
    return MRResult(method, float(est), float(se),
                    float(est - 1.96 * se), float(est + 1.96 * se), p, n)


def ivw(bx, sx, by, sy) -> MRResult:
    """Inverse-variance weighted (fixed-effect) estimator."""
    bx, sx, by, sy = map(np.asarray, (bx, sx, by, sy))
    w = 1.0 / sy**2
    est = np.sum(w * bx * by) / np.sum(w * bx**2)
    se = np.sqrt(1.0 / np.sum(w * bx**2))
    return _ci(est, se, "IVW", len(bx))


def mr_egger(bx, sx, by, sy) -> tuple[MRResult, dict]:
    """MR-Egger: slope = pleiotropy-robust causal estimate; intercept tests pleiotropy."""
    bx, sx, by, sy = map(np.asarray, (bx, sx, by, sy))
    # orient so all exposure effects are positive
    sign = np.sign(bx)
    bx_o, by_o = bx * sign, by * sign
    w = 1.0 / sy**2
    X = np.column_stack([np.ones_like(bx_o), bx_o])
    W = np.diag(w)
    xtwx = X.T @ W @ X
    beta = np.linalg.solve(xtwx, X.T @ W @ by_o)
    resid = by_o - X @ beta
    dof = max(len(bx_o) - 2, 1)
    sigma2 = (resid @ (w * resid)) / dof
    cov = sigma2 * np.linalg.inv(xtwx)
    intercept, slope = beta[0], beta[1]
    se_slope = np.sqrt(cov[1, 1])
    se_int = np.sqrt(cov[0, 0])
    p_int = float(2 * stats.t.sf(abs(intercept / se_int), dof)) if se_int > 0 else np.nan
    res = _ci(slope, se_slope, "MR-Egger", len(bx_o))
    return res, {"intercept": float(intercept), "intercept_se": float(se_int),
                 "intercept_p": p_int}


def weighted_median(bx, sx, by, sy, n_boot: int = 1000, seed: int = 0) -> MRResult:
    """Weighted-median estimator (consistent if >50% of weight is from valid instruments)."""
    bx, sx, by, sy = map(np.asarray, (bx, sx, by, sy))
    ratios = by / bx
    weights = (bx**2) / (sy**2)  # inverse-variance of the Wald ratio (approx)

    def wm(r, w):
        order = np.argsort(r)
        r, w = r[order], w[order]
        cw = np.cumsum(w) - 0.5 * w
        cw /= np.sum(w)
        k = np.searchsorted(cw, 0.5)
        k = min(max(k, 1), len(r) - 1)
        # linear interpolation between k-1 and k
        return r[k - 1] + (r[k] - r[k - 1]) * (0.5 - cw[k - 1]) / (cw[k] - cw[k - 1])

    est = wm(ratios, weights)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        bxi = rng.normal(bx, sx)
        byi = rng.normal(by, sy)
        boots.append(wm(byi / bxi, (bxi**2) / (sy**2)))
    se = float(np.std(boots))
    return _ci(est, se, "Weighted median", len(bx))


def cochran_q(bx, sx, by, sy) -> dict:
    """Cochran's Q heterogeneity around the IVW estimate (instrument invalidity signal)."""
    bx, sx, by, sy = map(np.asarray, (bx, sx, by, sy))
    b = ivw(bx, sx, by, sy).estimate
    w = 1.0 / sy**2
    q = float(np.sum(w * (by - b * bx) ** 2))
    dof = len(bx) - 1
    i2 = max(0.0, (q - dof) / q) if q > 0 else 0.0
    return {"Q": q, "df": dof, "p": float(stats.chi2.sf(q, dof)), "i2": i2}
