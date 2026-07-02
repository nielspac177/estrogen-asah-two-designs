"""Descriptive tables (7.1)."""

from __future__ import annotations

import pandas as pd


def _rate(s: pd.Series) -> float:
    return float(s.astype("boolean").astype("float").mean())


def table_one(df: pd.DataFrame, strata: str = "menopausal_stratum") -> pd.DataFrame:
    """Baseline characteristics by stratum (one column per stratum + overall)."""
    def summarize(sub: pd.DataFrame) -> dict:
        return {
            "n": len(sub),
            "age_mean": round(float(pd.to_numeric(sub["age"], errors="coerce").mean()), 1),
            "female_pct": round(100 * (sub["sex"] == "F").mean(), 1),
            "htn_pct": round(100 * _rate(sub["htn"]), 1),
            "smoking_pct": round(100 * _rate(sub["smoking"]), 1),
            "diabetes_pct": round(100 * _rate(sub["diabetes"]), 1),
            "apache_mean": round(
                float(pd.to_numeric(sub["apache_score"], errors="coerce").mean()), 1
            ),
            "dci_pct": round(100 * _rate(sub["dci_composite"]), 1),
            "mortality_pct": round(100 * _rate(sub["died"]), 1),
        }

    cols = {"overall": summarize(df)}
    for name, sub in df.groupby(strata):
        cols[str(name)] = summarize(sub)
    return pd.DataFrame(cols)
