"""MR estimators recover a known causal effect; robust methods resist pleiotropy."""

import math

from estrogen_mr.methods import cochran_q, ivw, mr_egger, weighted_median
from estrogen_mr.simulate import simulate_instruments

TRUE = math.log(1.5)  # true causal log-OR


def _d(**kw):
    return simulate_instruments(n_snps=80, theta=TRUE, seed=1, **kw)


def _arrs(df):
    return df["beta_exp"], df["se_exp"], df["beta_out"], df["se_out"]


def test_ivw_recovers_truth_no_pleiotropy():
    r = ivw(*_arrs(_d()))
    assert abs(r.estimate - TRUE) < 0.08
    assert r.ci_low < TRUE < r.ci_high
    assert r.n_snps == 80


def test_weighted_median_recovers_truth():
    r = weighted_median(*_arrs(_d()), n_boot=300, seed=2)
    assert abs(r.estimate - TRUE) < 0.12


def test_egger_intercept_null_when_no_pleiotropy():
    _, info = mr_egger(*_arrs(_d()))
    assert abs(info["intercept"]) < 0.02
    assert info["intercept_p"] > 0.05


def test_egger_detects_directional_pleiotropy():
    df = _d(pleiotropy_frac=0.4, pleiotropy_mean=0.15)
    ivw_est = ivw(*_arrs(df)).estimate
    egg, info = mr_egger(*_arrs(df))
    # IVW is biased upward by directional pleiotropy; Egger intercept flags it
    assert ivw_est > TRUE
    assert info["intercept"] > 0
    assert egg.estimate < ivw_est  # Egger slope corrects toward truth


def test_cochran_q_flags_heterogeneity():
    clean = cochran_q(*_arrs(_d()))
    dirty = cochran_q(*_arrs(_d(pleiotropy_frac=0.4, pleiotropy_mean=0.15)))
    assert dirty["Q"] > clean["Q"]
    assert dirty["i2"] >= clean["i2"]
