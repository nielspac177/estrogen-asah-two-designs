"""End-to-end cohort assembly shared by the numbered scripts.

If ``config/paths.yaml`` resolves to real, present PhysioNet data, the cohort is
extracted from MIMIC-IV + eICU. Otherwise the pipeline falls back to a clearly
labelled SYNTHETIC cohort (seeded simulator) so the full pipeline — and CI — runs
with zero credentialed data.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import REPO_ROOT, load_paths
from .features.derive import add_features

OUTPUTS = REPO_ROOT / "outputs"


def _real_data_available() -> bool:
    if not (REPO_ROOT / "config" / "paths.yaml").exists():
        return False
    try:
        p = load_paths()
        return p.mimic_iv_db.exists() and p.eicu_dir.exists()
    except Exception:
        return False


def build_featured_cohort() -> tuple[pd.DataFrame, dict]:
    """Return (featured cohort, meta). Uses real data if configured, else synthetic."""
    from . import synthetic as syn

    if _real_data_available():
        from .extract import eicu, mimic_iv
        from .harmonize.combine import combine

        p = load_paths()
        con = mimic_iv.connect(str(p.mimic_iv_db), read_only=True)
        m = mimic_iv.build_cohort(con)
        d = eicu.build_cohort(str(p.eicu_dir))
        cohort = combine(m, d)
        meta = {"mode": "real", "n": len(cohort),
                "n_mimic": len(m), "n_eicu": len(d)}
    else:
        cohort = syn.simulate_cohort(n=1500, seed=0)
        meta = {"mode": "synthetic", "n": len(cohort)}

    return add_features(cohort), meta


def ensure_outputs() -> Path:
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    return OUTPUTS
