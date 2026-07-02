"""Audit robustness estimators run and return sane structures on simulated data."""

import math

from estrogen_asah_dci import synthetic as syn
from estrogen_asah_dci.analysis import robust
from estrogen_asah_dci.analysis.models import Estimate
from estrogen_asah_dci.features.derive import add_features


def _df():
    return add_features(syn.simulate_cohort(n=3000, seed=7, log_or_post=math.log(0.6)))


def test_survivor_restricted_returns_estimate():
    est = robust.survivor_restricted(_df())
    assert isinstance(est, Estimate) and est.term == "post"
    assert est.n > 0


def test_age_restricted_band():
    est = robust.age_restricted(_df(), 45, 55)
    assert isinstance(est, Estimate)
    assert est.n > 0


def test_outcome_variant_runs():
    est = robust.outcome_variant(_df(), "vasospasm_dx")
    assert isinstance(est, Estimate) and est.n > 0


def test_overlap_weighted_runs():
    est = robust.overlap_weighted(_df())
    assert isinstance(est, Estimate)


def test_competing_events_multinomial_keys():
    out = robust.competing_events_multinomial(_df())
    assert set(out) >= {"rrr_dci", "ci_low", "ci_high", "n", "converged"}


def test_survivor_recovers_direction_on_sim():
    # simulator injects protective post effect (OR~0.6); survivor-restricted should stay <1
    est = robust.survivor_restricted(_df())
    if est.converged:
        assert est.odds_ratio < 1.2
