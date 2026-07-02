"""Multivariable IVW recovers known direct effects."""

import numpy as np

from estrogen_mr.mvmr import mvmr_ivw


def test_mvmr_recovers_two_effects():
    rng = np.random.default_rng(0)
    n = 300
    Bx = rng.normal(0, 0.1, (n, 2))          # two exposures' SNP effects
    true = np.array([0.30, -0.20])
    se_out = np.full(n, 0.02)
    by = Bx @ true + rng.normal(0, se_out)
    res = mvmr_ivw(Bx, by, se_out, ["x1", "x2"])
    assert res.n_snps == n
    d = res.as_dict()["exposures"]
    import math
    assert abs(math.log(d["x1"]["or"]) - 0.30) < 0.05
    assert abs(math.log(d["x2"]["or"]) - (-0.20)) < 0.05
