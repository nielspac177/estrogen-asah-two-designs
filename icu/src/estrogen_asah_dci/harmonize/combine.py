"""Pool the per-source cohorts into one harmonized frame (P5).

Sources are concatenated (never merged on patient), validated, and kept
distinguishable via ``source`` / ``hospital_id`` for stratified analysis
(ADR-0002).
"""

from __future__ import annotations

import pandas as pd

from .common_schema import COHORT_COLUMNS, coerce_cohort, validate_cohort


def combine(*cohorts: pd.DataFrame) -> pd.DataFrame:
    """Concatenate per-source cohorts, coerce to the schema, validate."""
    frames = [c for c in cohorts if c is not None and len(c)]
    if not frames:
        from .common_schema import empty_cohort
        return empty_cohort()
    combined = pd.concat([coerce_cohort(c) for c in frames], ignore_index=True)
    combined = combined[COHORT_COLUMNS]
    return validate_cohort(combined)
