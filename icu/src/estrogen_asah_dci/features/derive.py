"""Derived analysis features (P6): menopausal stratum + DCI composite.

These are added to the harmonized cohort; they are not part of the raw
extraction schema (ADR-0001).
"""

from __future__ import annotations

import pandas as pd

MENOPAUSE_AGE = 51  # ≈ median age at natural menopause


def menopausal_stratum(df: pd.DataFrame, age_cut: int = MENOPAUSE_AGE) -> pd.Series:
    """premenopausal (F, age<cut) / postmenopausal (F, age>=cut) / male / unknown."""
    def label(row):
        sex, age = row["sex"], row["age"]
        if sex == "M":
            return "male"
        if sex == "F" and pd.notna(age):
            return "premenopausal" if age < age_cut else "postmenopausal"
        return pd.NA

    return df.apply(label, axis=1).astype("string")


def dci_composite(df: pd.DataFrame) -> pd.Series:
    """Primary outcome: vasospasm dx OR rescue procedure OR delayed infarction (ADR-0003)."""
    return (
        df["vasospasm_dx"].fillna(False)
        | df["dci_procedure"].fillna(False)
        | df["delayed_infarction"].fillna(False)
    ).astype("boolean")


def add_features(df: pd.DataFrame, age_cut: int = MENOPAUSE_AGE) -> pd.DataFrame:
    """Return a copy with menopausal_stratum, dci_composite, and analysis helpers."""
    out = df.copy()
    out["menopausal_stratum"] = menopausal_stratum(out, age_cut)
    out["dci_composite"] = dci_composite(out)
    out["female"] = (out["sex"] == "F").astype("boolean")
    # centered age for models (within full cohort)
    if out["age"].notna().any():
        out["age_c"] = (out["age"] - out["age"].mean()).astype("Float64")
    else:
        out["age_c"] = pd.NA
    return out


def missingness(df: pd.DataFrame) -> pd.DataFrame:
    """Per-column missing fraction — a small data-quality report."""
    frac = df.isna().mean().sort_values(ascending=False)
    return frac.rename("missing_fraction").reset_index(names="column")
