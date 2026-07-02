"""QC / power / equivalence sanity checks."""

import math

from estrogen_mr import qc
from estrogen_mr.methods import ivw
from estrogen_mr.simulate import simulate_instruments


def _arrs(df):
    return df["beta_exp"], df["se_exp"], df["beta_out"], df["se_out"]


def test_random_effects_ci_not_narrower_than_fixed():
    df = simulate_instruments(80, theta=math.log(1.3), seed=1, pleiotropy_frac=0.3)
    fe = ivw(*_arrs(df))
    re = qc.ivw_random(*_arrs(df))
    assert re.se >= fe.se - 1e-9
    assert math.isclose(re.estimate, fe.estimate, rel_tol=1e-9)


def test_steiger_mostly_correct_for_strong_instruments():
    df = simulate_instruments(60, theta=math.log(1.4), seed=2)
    st = qc.steiger(*_arrs(df), n_exp=200000, n_out=17000)
    assert st["frac_correct"] > 0.5   # exposure r^2 exceeds outcome r^2 for most valid instruments


def test_power_mde_monotone_and_bounds():
    p = qc.power_mde(se=0.03, sd_units=4.0)
    assert p["mde_or_perSD_protect"] < 1 < p["mde_or_perSD_harm"]
    pw = p["power_perSD_at"]
    assert pw["OR=0.8"] > pw["OR=0.9"]        # larger effects easier to detect


def test_tost_excludes_strong_protection_when_estimate_near_null():
    # estimate ~ +0.11 per SD, se ~0.11 -> reject strong protection, not strong harm
    eq = qc.tost(0.114, 0.115, sesoi_or=0.90)
    assert eq["p_reject_strong_protection"] < 0.05
    assert eq["p_reject_strong_harm"] > 0.05
    assert eq["equivalent_to_null"] is False


def test_leave_one_out_length():
    df = simulate_instruments(20, theta=0.0, seed=3)
    assert len(qc.leave_one_out(*_arrs(df))) == 20
