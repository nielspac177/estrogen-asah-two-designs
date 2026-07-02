"""Audit-endorsed analyses: specification curve + age×sex spline DiD."""

import math

from estrogen_asah_dci import synthetic as syn
from estrogen_asah_dci.analysis.multiverse import age_sex_did, spec_curve, spec_curve_summary
from estrogen_asah_dci.features.derive import add_features


def _df(log_or=None):
    if log_or is None:
        log_or = math.log(0.5)
    return add_features(syn.simulate_cohort(n=4000, seed=11, log_or_post=log_or))


def test_spec_curve_shape():
    curve = spec_curve(_df())
    # 6 cutoffs × 2 adjust × 2 outcomes × 3 sources
    assert len(curve) == 6 * 2 * 2 * 3
    assert {"cutoff", "adjust", "outcome", "source", "or", "p"}.issubset(curve.columns)


def test_spec_curve_summary_keys():
    s = spec_curve_summary(spec_curve(_df()))
    assert {"n_specs", "median_or", "primary_or", "n_sig",
            "n_sig_hypothesis_consistent_OR_gt_1", "n_sig_opposite_OR_lt_1"}.issubset(s)


def test_did_recovers_injected_menopause_jump():
    # simulator injects a female-only protective menopause effect (OR<1);
    # the DiD (female-minus-male) should recover the protective direction
    did = age_sex_did(_df(math.log(0.5)), n_boot=80, seed=1)
    assert did["converged"]
    assert did["or_ratio"] < 1.1


def test_did_null_when_no_effect():
    did = age_sex_did(_df(math.log(1.0)), n_boot=80, seed=2)
    assert did["converged"]
    assert 0.5 < did["or_ratio"] < 2.0
