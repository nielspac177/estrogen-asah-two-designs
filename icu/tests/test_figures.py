"""P8 — figures render to non-empty PNGs."""

import pandas as pd

from estrogen_asah_dci.viz import figures


def test_forest_plot(tmp_path):
    rows = [
        {"label": "post vs pre", "or": 0.68, "ci_low": 0.44, "ci_high": 1.04},
        {"label": "male vs female", "or": 1.4, "ci_low": 1.1, "ci_high": 1.8},
        {"label": "nan drop", "or": float("nan"), "ci_low": 1, "ci_high": 2},
    ]
    p = figures.forest_plot(rows, tmp_path)
    assert p.exists() and p.stat().st_size > 1000


def test_rates_plot(tmp_path):
    rates = pd.DataFrame({
        "menopausal_stratum": ["premenopausal", "postmenopausal", "male"],
        "n": [259, 678, 563], "events": [74, 162, 243],
        "rate": [0.286, 0.239, 0.432],
    })
    p = figures.rates_plot(rates, tmp_path)
    assert p.exists() and p.stat().st_size > 1000


def test_love_plot(tmp_path):
    bal = pd.DataFrame({
        "covariate": ["htn", "smoking", "diabetes"],
        "smd_unweighted": [0.09, 0.02, 0.01], "smd_weighted": [0.001, 0.004, 0.004],
    })
    p = figures.love_plot(bal, tmp_path)
    assert p.exists() and p.stat().st_size > 1000


def test_cohort_flow(tmp_path):
    p = figures.cohort_flow([("aSAH", 1500), ("Women", 937), ("Pre", 259), ("Post", 678)], tmp_path)
    assert p.exists() and p.stat().st_size > 1000
