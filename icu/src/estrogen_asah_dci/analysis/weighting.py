"""Propensity-score weighting and covariate balance (7.4).

Estimates a propensity for postmenopausal status (women only) and forms IPW
weights (ATE or ATT). Balance is summarized with standardized mean differences
(target |SMD| < 0.1).
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

# NB: age is deliberately EXCLUDED — menopausal stratum is defined by an age
# threshold, so age is collinear with the exposure (no overlap). Balance is
# assessed on the remaining measured confounders; residual age/menopause
# confounding is addressed by the age-restricted sensitivity analysis, not IPW.
DEFAULT_COVARIATES = ["htn", "smoking", "diabetes"]


def propensity_weights(df: pd.DataFrame, covariates: list[str] | None = None,
                       estimand: str = "ATE") -> pd.DataFrame:
    """Return the women pre/post subset with columns ``ps`` and ``ipw``."""
    covariates = covariates or DEFAULT_COVARIATES
    d = df[df["menopausal_stratum"].isin(["premenopausal", "postmenopausal"])].copy()
    d["t"] = (d["menopausal_stratum"] == "postmenopausal").astype(float)
    for c in covariates:
        d[c] = pd.to_numeric(d[c].astype("float") if str(d[c].dtype) in ("boolean", "bool")
                             else d[c], errors="coerce")
    d = d.dropna(subset=["t", *covariates])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = smf.logit("t ~ " + " + ".join(covariates), d).fit(disp=0, maxiter=100)
    d["ps"] = res.predict(d).clip(1e-3, 1 - 1e-3)
    if estimand == "ATT":
        d["ipw"] = np.where(d["t"] == 1, 1.0, d["ps"] / (1 - d["ps"]))
    else:  # ATE
        d["ipw"] = np.where(d["t"] == 1, 1 / d["ps"], 1 / (1 - d["ps"]))
    return d


def standardized_diff(d: pd.DataFrame, var: str, weights: str | None = None) -> float:
    """|SMD| for ``var`` between t=1 and t=0 groups, optionally weighted."""
    x = pd.to_numeric(d[var].astype("float") if str(d[var].dtype) in ("boolean", "bool")
                      else d[var], errors="coerce")
    t = d["t"]
    w = d[weights] if weights else pd.Series(1.0, index=d.index)

    def wm(mask):
        ww = w[mask]
        return np.average(x[mask], weights=ww)

    def wv(mask, m):
        ww = w[mask]
        return np.average((x[mask] - m) ** 2, weights=ww)

    m1, m0 = wm(t == 1), wm(t == 0)
    v1, v0 = wv(t == 1, m1), wv(t == 0, m0)
    denom = np.sqrt((v1 + v0) / 2) or np.nan
    return float(abs(m1 - m0) / denom)


def balance_table(df: pd.DataFrame, covariates: list[str] | None = None,
                  estimand: str = "ATE") -> pd.DataFrame:
    """|SMD| per covariate, unweighted vs IPW-weighted."""
    covariates = covariates or DEFAULT_COVARIATES
    d = propensity_weights(df, covariates, estimand)
    rows = [
        {
            "covariate": c,
            "smd_unweighted": standardized_diff(d, c, None),
            "smd_weighted": standardized_diff(d, c, "ipw"),
        }
        for c in covariates
    ]
    return pd.DataFrame(rows)
