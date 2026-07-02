"""P7 — E-value/meta exact math; models recover a known effect on simulated data."""

import math

import numpy as np

from estrogen_asah_dci import synthetic as syn
from estrogen_asah_dci.analysis import models
from estrogen_asah_dci.analysis.evalue import evalue_or, evalue_rr, or_to_rr
from estrogen_asah_dci.analysis.meta import pool_or
from estrogen_asah_dci.analysis.weighting import balance_table
from estrogen_asah_dci.features.derive import add_features


# ---- E-value (exact) ----
def test_evalue_point_known_values():
    assert evalue_rr(2.0).point == 2 + math.sqrt(2 * 1)          # 3.414...
    # protective RR mirrors: E-value(0.5) == E-value(2.0)
    assert math.isclose(evalue_rr(0.5).point, evalue_rr(2.0).point)


def test_evalue_ci_crossing_null_is_one():
    assert evalue_rr(1.3, lo=0.8, hi=2.1).ci_bound == 1.0


def test_or_to_rr_rare_outcome_approx_equal_or():
    assert math.isclose(or_to_rr(2.0, 0.001), 2.0, rel_tol=1e-2)


def test_evalue_or_protective():
    ev = evalue_or(0.67, baseline_risk=0.2, lo=0.42, hi=1.05)
    assert ev.point > 1
    assert ev.ci_bound == 1.0  # CI crosses null


# ---- meta pooling ----
def test_pool_two_identical_estimates():
    e = {"or": 0.6, "ci_low": 0.4, "ci_high": 0.9}
    pooled = pool_or([e, e], random=False)
    assert math.isclose(pooled.estimate, 0.6, rel_tol=1e-6)
    assert pooled.k == 2
    assert pooled.ci_high - pooled.ci_low < (0.9 - 0.4)  # narrower than a single study


# ---- models recover a known effect ----
def test_primary_recovers_known_or():
    df = add_features(syn.simulate_cohort(n=6000, seed=1, log_or_post=math.log(0.5)))
    est = models.primary_logistic(df)
    assert est.converged
    assert 0.35 < est.odds_ratio < 0.72          # recovers ~0.5 (protective)
    assert est.ci_low < est.odds_ratio < est.ci_high


def test_crude_rates_cover_strata():
    df = add_features(syn.simulate_cohort(n=1000, seed=2))
    rates = models.crude_dci_rates(df)
    assert set(rates["menopausal_stratum"]) >= {"premenopausal", "postmenopausal", "male"}
    assert (rates["rate"].between(0, 1)).all()


def test_by_source_and_pool():
    df = add_features(syn.simulate_cohort(n=6000, seed=3, log_or_post=math.log(0.6)))
    ests = {k: v for k, v in models.by_source(df).items() if v.converged}
    pooled = pool_or([e.as_dict() for e in ests.values()], random=True)
    assert 0.4 < pooled.estimate < 0.85


def test_balance_table_shape_and_finite():
    df = add_features(syn.simulate_cohort(n=2000, seed=4))
    bt = balance_table(df)
    assert list(bt["covariate"]) == ["htn", "smoking", "diabetes"]
    assert np.isfinite(bt["smd_weighted"]).all()
    assert (bt["smd_weighted"] < 0.3).all()


def test_models_dont_crash_on_tiny_cohort():
    # the 18-row real-synthetic cohort may not converge, but must not raise
    con = syn.load_mimic_duckdb()
    from estrogen_asah_dci.extract import mimic_iv
    df = add_features(mimic_iv.build_cohort(con))
    est = models.primary_logistic(df)
    assert est.n >= 0 and est.term == "post"
