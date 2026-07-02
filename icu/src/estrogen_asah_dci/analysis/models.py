"""Regression models for the DCI analysis.

Primary: logistic regression of the DCI composite on menopausal stratum
(postmenopausal vs premenopausal, women), adjusted, with cluster-robust standard
errors by hospital (ADR-0002). Fits are wrapped so non-convergence/separation on
small cohorts returns NaN estimates rather than raising — the pipeline always
completes and reports what could be estimated.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

# AUDIT FIX (ADR-0005): age is a deterministic driver of the age-defined menopause
# exposure, so co-adjusting a smooth age term with the menopause step is a
# functional-form artifact (false-precision CI) and is prohibited. The outcome
# model uses the same measured confounders as the propensity model (weighting.py).
DEFAULT_COVARIATES = ["htn", "smoking", "diabetes"]


@dataclass
class Estimate:
    term: str
    odds_ratio: float
    ci_low: float
    ci_high: float
    p: float
    n: int
    events: int
    converged: bool

    def as_dict(self) -> dict:
        return {
            "term": self.term, "or": self.odds_ratio,
            "ci_low": self.ci_low, "ci_high": self.ci_high,
            "p": self.p, "n": self.n, "events": self.events,
            "converged": self.converged,
        }


def _prep(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    d = df.copy()
    for c in cols:
        if str(d[c].dtype) in ("boolean", "bool"):
            d[c] = d[c].astype("float")
        else:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


def _fit(d: pd.DataFrame, formula: str, term: str, cluster: str | None) -> Estimate:
    y = d["y"]
    n, events = len(d), int(np.nansum(y.values))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = smf.logit(formula, d)
            if cluster and cluster in d and d[cluster].nunique() > 1:
                res = model.fit(disp=0, maxiter=100,
                                cov_type="cluster", cov_kwds={"groups": d[cluster]})
            else:
                res = model.fit(disp=0, maxiter=100)
        if term not in res.params:
            raise ValueError(f"term {term} absent")
        coef = res.params[term]
        ci = res.conf_int().loc[term]
        return Estimate(term, float(np.exp(coef)), float(np.exp(ci[0])),
                        float(np.exp(ci[1])), float(res.pvalues[term]), n, events, True)
    except Exception:
        return Estimate(term, float("nan"), float("nan"), float("nan"),
                        float("nan"), n, events, False)


def primary_logistic(df: pd.DataFrame, covariates: list[str] | None = None,
                     cluster: str = "hospital_id") -> Estimate:
    """Postmenopausal vs premenopausal (women) on the DCI composite (7.3)."""
    covariates = covariates or DEFAULT_COVARIATES
    d = df[df["menopausal_stratum"].isin(["premenopausal", "postmenopausal"])].copy()
    d["post"] = (d["menopausal_stratum"] == "postmenopausal").astype(float)
    d["y"] = d["dci_composite"].astype("boolean").astype("float")
    d = _prep(d, covariates + ["y", "post"])
    d = d.dropna(subset=["y", "post", *covariates])
    formula = "y ~ post + " + " + ".join(covariates)
    return _fit(d, formula, "post", cluster)


def sex_difference(df: pd.DataFrame, covariates: list[str] | None = None,
                   cluster: str = "hospital_id") -> Estimate:
    """Secondary: male vs female on the DCI composite (7.7)."""
    covariates = covariates or DEFAULT_COVARIATES
    d = df.copy()
    d["male"] = (df["sex"] == "M").astype(float)
    d["y"] = d["dci_composite"].astype("boolean").astype("float")
    d = _prep(d, covariates + ["y", "male"])
    d = d.dropna(subset=["y", "male", *covariates])
    formula = "y ~ male + " + " + ".join(covariates)
    return _fit(d, formula, "male", cluster)


def exposure_logistic(df: pd.DataFrame, exposure: str = "hrt_exposure",
                      covariates: list[str] | None = None,
                      cluster: str = "hospital_id") -> Estimate:
    """Exploratory: an exposure flag (e.g. HRT) on the DCI composite (7.7)."""
    covariates = covariates or DEFAULT_COVARIATES
    d = df.copy()
    d["exp"] = d[exposure].astype("boolean").astype("float")
    d["y"] = d["dci_composite"].astype("boolean").astype("float")
    d = _prep(d, covariates + ["y", "exp"])
    d = d.dropna(subset=["y", "exp", *covariates])
    formula = "y ~ exp + " + " + ".join(covariates)
    return _fit(d, formula, "exp", cluster)


def by_source(df: pd.DataFrame, **kwargs) -> dict[str, Estimate]:
    """Primary estimate within each data source, for meta-analytic pooling."""
    return {
        src: primary_logistic(sub, **kwargs)
        for src, sub in df.groupby("source")
    }


def crude_dci_rates(df: pd.DataFrame) -> pd.DataFrame:
    """DCI rate by menopausal stratum (7.2)."""
    d = df.copy()
    d["dci"] = d["dci_composite"].astype("boolean").astype("float")
    g = d.groupby("menopausal_stratum")["dci"]
    return pd.DataFrame({"n": g.size(), "events": g.sum(), "rate": g.mean()}).reset_index()
