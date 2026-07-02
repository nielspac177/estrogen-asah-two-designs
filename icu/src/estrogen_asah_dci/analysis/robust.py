"""Alternative / robustness estimators for the menopause -> DCI contrast.

Motivated by the competing-mortality problem: postmenopausal women die more, so a
naive logistic on DCI is distorted by differential survival (DCI can only be
coded if the patient survives long enough). These estimators attack that and
other researcher-degrees-of-freedom in *pre-specified*, mechanism-driven ways.
Every result is reported — none is cherry-picked.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from .models import DEFAULT_COVARIATES, Estimate, _fit, _prep


def _women_prepost(df: pd.DataFrame) -> pd.DataFrame:
    d = df[df["menopausal_stratum"].isin(["premenopausal", "postmenopausal"])].copy()
    d["post"] = (d["menopausal_stratum"] == "postmenopausal").astype(float)
    return d


def survivor_restricted(df, covariates=None, cluster="hospital_id") -> Estimate:
    """Competing-mortality fix #1: DCI among hospital survivors only (landmark)."""
    covariates = covariates or DEFAULT_COVARIATES
    d = _women_prepost(df)
    d = d[~d["died"].fillna(False)]
    d["y"] = d["dci_composite"].astype("boolean").astype("float")
    d = _prep(d, covariates + ["y", "post"]).dropna(subset=["y", "post", *covariates])
    return _fit(d, "y ~ post + " + " + ".join(covariates), "post", cluster)


def age_restricted(df, lo=45, hi=55, covariates=None, cluster="hospital_id") -> Estimate:
    """Isolate menopause from ageing: restrict to a peri-menopausal age band."""
    covariates = covariates or DEFAULT_COVARIATES
    d = _women_prepost(df)
    age = pd.to_numeric(d["age"], errors="coerce")
    d = d[(age >= lo) & (age <= hi)]
    d["y"] = d["dci_composite"].astype("boolean").astype("float")
    d = _prep(d, covariates + ["y", "post"]).dropna(subset=["y", "post", *covariates])
    return _fit(d, "y ~ post + " + " + ".join(covariates), "post", cluster)


def outcome_variant(df, which="vasospasm_dx", covariates=None, cluster="hospital_id") -> Estimate:
    """Vary the outcome (ADR-0003 sensitivity): vasospasm-only / procedure-only / composite."""
    covariates = covariates or DEFAULT_COVARIATES
    d = _women_prepost(df)
    d["y"] = d[which].astype("boolean").astype("float")
    d = _prep(d, covariates + ["y", "post"]).dropna(subset=["y", "post", *covariates])
    return _fit(d, "y ~ post + " + " + ".join(covariates), "post", cluster)


def overlap_weighted(df, covariates=None) -> Estimate:
    """Doubly-guarded contrast via overlap weights (ATO; exact balance, no extrapolation)."""
    covariates = covariates or DEFAULT_COVARIATES
    d = _women_prepost(df)
    d["y"] = d["dci_composite"].astype("boolean").astype("float")
    d = _prep(d, covariates + ["y", "post"]).dropna(subset=["y", "post", *covariates])
    n, events = len(d), int(d["y"].sum())
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ps_fit = smf.logit("post ~ " + " + ".join(covariates), d).fit(disp=0, maxiter=100)
            ps = ps_fit.predict(d)
            d["ow"] = np.where(d["post"] == 1, 1 - ps, ps)  # overlap weights
            res = smf.glm("y ~ post", d, freq_weights=d["ow"],
                          family=sm.families.Binomial()).fit()
        coef, ci = res.params["post"], res.conf_int().loc["post"]
        return Estimate("post", float(np.exp(coef)), float(np.exp(ci[0])),
                        float(np.exp(ci[1])), float(res.pvalues["post"]), n, events, True)
    except Exception:
        return Estimate("post", float("nan"), float("nan"), float("nan"),
                        float("nan"), n, events, False)


def competing_events_multinomial(df, covariates=None) -> dict:
    """Cross-sectional competing-events model: outcome in {neither, DCI, death-without-DCI}.

    Returns the relative-risk ratio for postmenopausal on the DCI category,
    holding the death category as a separate competing outcome (not a collider).
    """
    covariates = covariates or DEFAULT_COVARIATES
    d = _women_prepost(df)
    dci = d["dci_composite"].astype("boolean").fillna(False)
    died = d["died"].astype("boolean").fillna(False)
    # 0 neither, 1 DCI (any), 2 death without DCI
    d["cat"] = np.where(dci, 1, np.where(died, 2, 0)).astype(int)
    d = _prep(d, covariates + ["post"]).dropna(subset=["post", *covariates])
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = smf.mnlogit("cat ~ post + " + " + ".join(covariates), d).fit(disp=0, maxiter=200)
        # column for outcome==1 (DCI) in the params (base category 0)
        params = res.params
        ci = res.conf_int()
        # statsmodels indexes mnlogit params by outcome index 1..K-1
        col = params.columns[0]  # first non-base outcome == category 1 (DCI)
        rrr = float(np.exp(params.loc["post", col]))
        lo = float(np.exp(ci.loc[(col, "post"), 0])) if (col, "post") in ci.index else float("nan")
        hi = float(np.exp(ci.loc[(col, "post"), 1])) if (col, "post") in ci.index else float("nan")
        return {"rrr_dci": rrr, "ci_low": lo, "ci_high": hi, "n": len(d), "converged": True}
    except Exception as e:  # noqa: BLE001
        return {"rrr_dci": float("nan"), "ci_low": float("nan"), "ci_high": float("nan"),
                "n": len(d), "converged": False, "error": str(e)}
