"""Honest multiverse: specification curve + age×sex spline difference-in-differences.

These are the two analyses the adversarial audit endorsed (docs/adr/0005). The
specification curve reports the menopause OR across ALL defensible forks so the
reader sees where the primary sits and how many forks reach significance in each
direction (the antidote to fork-selection). The age×sex spline DiD is the one
design that structurally confronts age–menopause collinearity by borrowing men as
a chronological-ageing reference.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import patsy
import statsmodels.formula.api as smf


def _logit_or(d: pd.DataFrame, covariates: list[str]) -> dict:
    d = d.copy()
    d["y"] = d["y"].astype("float")
    for c in covariates:
        d[c] = pd.to_numeric(
            d[c].astype("float") if str(d[c].dtype) in ("boolean", "bool") else d[c],
            errors="coerce",
        )
    d = d.dropna(subset=["y", "post", *covariates])
    formula = "y ~ post" + ("" if not covariates else " + " + " + ".join(covariates))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = smf.logit(formula, d).fit(disp=0, maxiter=100)
        return {"or": float(np.exp(res.params["post"])), "p": float(res.pvalues["post"]),
                "n": len(d), "events": int(d["y"].sum()), "converged": True}
    except Exception:
        return {"or": float("nan"), "p": float("nan"), "n": len(d),
                "events": int(d["y"].sum()), "converged": False}


def spec_curve(df: pd.DataFrame) -> pd.DataFrame:
    """Post-vs-pre OR across cutoffs × covariate sets × outcomes × sources."""
    cutoffs = [45, 48, 50, 51, 52, 55]
    covsets = {"unadj": [], "clinical": ["htn", "smoking", "diabetes"]}
    outcomes = ["dci_composite", "vasospasm_dx"]
    sources = {"all": None, "mimic_iv": "mimic_iv", "eicu": "eicu"}
    rows = []
    women = df[df["sex"] == "F"].copy()
    women["_age"] = pd.to_numeric(women["age"], errors="coerce")
    for cut in cutoffs:
        for cname, cov in covsets.items():
            for out in outcomes:
                for sname, src in sources.items():
                    d = women if src is None else women[women["source"] == src]
                    d = d.dropna(subset=["_age"]).copy()
                    d["post"] = (d["_age"] >= cut).astype(float)
                    d["y"] = d[out].astype("boolean")
                    r = _logit_or(d, cov)
                    rows.append({"cutoff": cut, "adjust": cname, "outcome": out,
                                 "source": sname, **r})
    return pd.DataFrame(rows)


def spec_curve_summary(curve: pd.DataFrame) -> dict:
    """Where does the primary sit; how many forks are significant, and which way."""
    conv = curve[curve["converged"]]
    sig = conv[conv["p"] < 0.05]
    primary = conv[(conv.cutoff == 51) & (conv.adjust == "clinical")
                   & (conv.outcome == "dci_composite") & (conv.source == "all")]
    return {
        "n_specs": int(len(curve)),
        "n_converged": int(len(conv)),
        "median_or": float(conv["or"].median()),
        "primary_or": float(primary["or"].iloc[0]) if len(primary) else float("nan"),
        "n_sig": int(len(sig)),
        "n_sig_hypothesis_consistent_OR_gt_1": int((sig["or"] > 1).sum()),
        "n_sig_opposite_OR_lt_1": int((sig["or"] < 1).sum()),
    }


def age_sex_did(df: pd.DataFrame, young: int = 45, old: int = 55, n_boot: int = 500,
                seed: int = 0) -> dict:
    """Female-minus-male change in log-odds of DCI from `young` to `old` age.

    A menopause-attributable signal is a female-specific upward inflection after
    ~50 that is NOT shared by men (who age without an abrupt oestrogen transition).
    Estimated on the log-odds scale; CI by nonparametric bootstrap.
    """
    d = df.copy()
    d["y"] = d["dci_composite"].astype("boolean").astype("float")
    d["age"] = pd.to_numeric(d["age"], errors="coerce")
    d["female"] = (d["sex"] == "F").astype(float)
    d = d.dropna(subset=["y", "age", "female"])
    formula = "y ~ bs(age, df=4) * female"

    def did_once(data) -> float:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = smf.logit(formula, data).fit(disp=0, maxiter=100)
        di = res.model.data.design_info
        grid = pd.DataFrame({"age": [young, old, young, old],
                             "female": [1.0, 1.0, 0.0, 0.0]})
        X = patsy.build_design_matrices([di], grid)[0]
        lp = np.asarray(X) @ res.params.values
        return float((lp[1] - lp[0]) - (lp[3] - lp[2]))  # (F_old-F_young)-(M_old-M_young)

    try:
        point = did_once(d)
    except Exception as e:  # noqa: BLE001
        return {"converged": False, "error": str(e), "n": len(d)}

    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        samp = d.iloc[rng.integers(0, len(d), len(d))]
        try:
            boots.append(did_once(samp))
        except Exception:
            continue
    lo, hi = (np.nanpercentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan))
    return {
        "did_logodds": point, "or_ratio": float(np.exp(point)),
        "ci_low": float(np.exp(lo)), "ci_high": float(np.exp(hi)),
        "n": len(d), "n_boot": len(boots), "converged": True,
        "interpretation": "OR-ratio >1 = female-specific rise in DCI odds after menopause "
                          "(hypothesis-consistent); ~1 or <1 = no menopause inflection.",
    }
