"""The harmonized cohort schema shared by both sources.

Every extractor (MIMIC-IV, eICU) must emit a DataFrame with exactly these
columns and dtypes — one row per patient-admission. Features (menopausal
stratum, DCI composite) are added downstream in ``features/`` and are NOT part
of the raw extraction schema.
"""

from __future__ import annotations

import pandas as pd

# column -> pandas dtype for the raw harmonized cohort (pre-features)
COHORT_SCHEMA: dict[str, str] = {
    "source": "string",            # "mimic_iv" | "eicu"
    "patient_id": "string",        # unique within source
    "hospital_id": "string",       # for clustering; "mimic_iv" is a single site
    "age": "Float64",              # years (eICU ">89" -> 90)
    "sex": "string",               # "M" | "F"
    # outcome components
    "vasospasm_dx": "boolean",
    "dci_procedure": "boolean",
    "delayed_infarction": "boolean",
    # exposure (secondary) + covariates
    "hrt_exposure": "boolean",
    "aneurysm_secured": "boolean",
    "treatment_modality": "string",   # "clip" | "coil" | <NA>
    "htn": "boolean",
    "smoking": "boolean",
    "diabetes": "boolean",
    "severity_gcs": "Float64",        # lowest GCS if available
    "apache_score": "Float64",        # eICU APACHE IV; <NA> for MIMIC
    # secondary outcomes
    "died": "boolean",                # in-hospital mortality
    "poor_disposition": "boolean",    # death or hospice
    "icu_los_days": "Float64",
}

COHORT_COLUMNS: list[str] = list(COHORT_SCHEMA)


def empty_cohort() -> pd.DataFrame:
    """An empty cohort frame with the correct columns/dtypes."""
    return pd.DataFrame({c: pd.Series(dtype=t) for c, t in COHORT_SCHEMA.items()})


def coerce_cohort(df: pd.DataFrame) -> pd.DataFrame:
    """Reorder to canonical columns and coerce dtypes; raise on missing columns."""
    missing = set(COHORT_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"cohort missing columns: {sorted(missing)}")
    out = df[COHORT_COLUMNS].copy()
    for col, dtype in COHORT_SCHEMA.items():
        out[col] = out[col].astype(dtype)
    return out


def validate_cohort(df: pd.DataFrame) -> pd.DataFrame:
    """Assert schema + basic plausibility; return the frame unchanged."""
    coerce_cohort(df)  # column/dtype check
    if len(df):
        assert df["sex"].dropna().isin(["M", "F"]).all(), "sex must be M/F"
        assert df["source"].dropna().isin(["mimic_iv", "eicu"]).all()
        ages = df["age"].dropna()
        assert (ages >= 0).all() and (ages <= 120).all(), "age out of range"
        assert df["patient_id"].is_unique, "patient_id must be unique within cohort"
    return df
