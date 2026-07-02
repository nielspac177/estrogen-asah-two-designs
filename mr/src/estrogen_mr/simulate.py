"""Synthetic instruments with a known causal effect — for estimator validation.

Generates two-sample summary statistics under a specified true causal log-OR,
optionally with a fraction of pleiotropic (invalid) instruments so that MR-Egger
and the weighted median can be checked for robustness where IVW is biased.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def simulate_instruments(n_snps: int = 40, theta: float = 0.0, seed: int = 0,
                         pleiotropy_frac: float = 0.0, pleiotropy_mean: float = 0.1,
                         se_exp: float = 0.02, se_out: float = 0.03) -> pd.DataFrame:
    """Return harmonized per-SNP arrays (beta_exp/se_exp/beta_out/se_out).

    True model: beta_out = theta * beta_exp + alpha (pleiotropy) + noise.
    """
    rng = np.random.default_rng(seed)
    gamma = rng.uniform(0.05, 0.25, n_snps)           # true SNP->exposure effects
    beta_exp = rng.normal(gamma, se_exp)
    alpha = np.zeros(n_snps)
    if pleiotropy_frac > 0:
        k = int(round(pleiotropy_frac * n_snps))
        idx = rng.choice(n_snps, k, replace=False)
        alpha[idx] = rng.normal(pleiotropy_mean, 0.02, k)  # directional pleiotropy
    beta_out = rng.normal(theta * gamma + alpha, se_out)
    return pd.DataFrame({
        "SNP": [f"rs{i}" for i in range(n_snps)],
        "beta_exp": beta_exp, "se_exp": np.full(n_snps, se_exp),
        "beta_out": beta_out, "se_out": np.full(n_snps, se_out),
    })
