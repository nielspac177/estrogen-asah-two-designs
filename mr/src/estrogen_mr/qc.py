"""MR quality-control, sensitivity, power, and equivalence testing.

Everything here runs on summary statistics alone (no reference panel):
random-effects IVW, Steiger directionality, leave-one-out, a global MR-PRESSO
heterogeneity/outlier test, power / minimum-detectable-effect, and TOST
equivalence against a smallest-effect-of-interest.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from .methods import MRResult, ivw


def ivw_random(bx, sx, by, sy) -> MRResult:
    """IVW with multiplicative random-effects SE (inflate by sqrt(Q/df) if >1)."""
    bx, sx, by, sy = map(np.asarray, (bx, sx, by, sy))
    fe = ivw(bx, sx, by, sy)
    w = 1.0 / sy**2
    q = float(np.sum(w * (by - fe.estimate * bx) ** 2))
    df = len(bx) - 1
    phi = max(1.0, (q / df) ** 0.5)
    se = fe.se * phi
    est = fe.estimate
    z = est / se
    return MRResult("IVW (random effects)", est, se, est - 1.96 * se, est + 1.96 * se,
                    float(2 * stats.norm.sf(abs(z))), len(bx))


def _r2_from_z(beta, se, n):
    f = (np.asarray(beta) / np.asarray(se)) ** 2
    return f / (f + n - 2)


def steiger(bx, sx, by, sy, n_exp, n_out) -> dict:
    """Steiger directionality: keep SNPs where exposure r^2 > outcome r^2."""
    r2x = _r2_from_z(bx, sx, n_exp)
    r2y = _r2_from_z(by, sy, n_out)
    correct = r2x > r2y
    return {"n": len(bx), "n_correct": int(correct.sum()),
            "frac_correct": float(correct.mean()), "keep_mask": correct}


def leave_one_out(bx, sx, by, sy) -> list[dict]:
    bx, sx, by, sy = map(np.asarray, (bx, sx, by, sy))
    out = []
    for i in range(len(bx)):
        m = np.arange(len(bx)) != i
        r = ivw(bx[m], sx[m], by[m], sy[m])
        out.append({"dropped": i, "or": r.odds_ratio, "p": r.p})
    return out


def mr_presso_global(bx, sx, by, sy, n_sim: int = 2000, seed: int = 0) -> dict:
    """Global MR-PRESSO-style test: observed RSS vs its null distribution."""
    bx, sx, by, sy = map(np.asarray, (bx, sx, by, sy))
    w = 1.0 / sy**2

    def rss(bx_, by_):
        b = np.sum(w * bx_ * by_) / np.sum(w * bx_**2)
        return np.sum(w * (by_ - b * bx_) ** 2)

    obs = rss(bx, by)
    bhat = np.sum(w * bx * by) / np.sum(w * bx**2)
    rng = np.random.default_rng(seed)
    sims = np.empty(n_sim)
    for k in range(n_sim):
        by_sim = rng.normal(bhat * bx, sy)
        sims[k] = rss(bx, by_sim)
    p = float((np.sum(sims >= obs) + 1) / (n_sim + 1))
    return {"rss_obs": float(obs), "p_global": p}


def power_mde(se: float, sd_units: float = 4.0, alpha: float = 0.05,
              power: float = 0.80, effects_or=(0.80, 0.85, 0.90)) -> dict:
    """Minimum detectable effect + power, on the estimate's native scale and per SD."""
    za, zb = stats.norm.ppf(1 - alpha / 2), stats.norm.ppf(power)
    mde_log = (za + zb) * se
    se_sd = se * sd_units

    def pw(or_target, se_):  # power to detect a given OR (two-sided)
        theta = abs(np.log(or_target))
        return float(stats.norm.cdf(theta / se_ - za))

    return {
        "se_per_unit": se, "se_per_sd": se_sd, "sd_units": sd_units,
        "mde_or_native_protect": float(np.exp(-mde_log)),
        "mde_or_native_harm": float(np.exp(mde_log)),
        "mde_or_perSD_protect": float(np.exp(-(za + zb) * se_sd)),
        "mde_or_perSD_harm": float(np.exp((za + zb) * se_sd)),
        "power_perSD_at": {f"OR={o}": pw(o, se_sd) for o in effects_or},
    }


def tost(estimate: float, se: float, sesoi_or: float = 0.90) -> dict:
    """Two one-sided tests of equivalence to the null within +/- log(SESOI)."""
    delta = abs(np.log(sesoi_or))
    # H0a: theta <= -delta (real protection) ; reject if estimate sufficiently above -delta
    z_lower = (estimate - (-delta)) / se
    p_lower = float(stats.norm.sf(z_lower))          # small -> reject strong protection
    # H0b: theta >= +delta (real harm) ; reject if estimate sufficiently below +delta
    z_upper = (delta - estimate) / se
    p_upper = float(stats.norm.sf(z_upper))          # small -> reject strong harm
    return {"sesoi_or": sesoi_or, "delta_log": delta,
            "p_reject_strong_protection": p_lower,
            "p_reject_strong_harm": p_upper,
            "equivalent_to_null": bool(p_lower < 0.05 and p_upper < 0.05)}
