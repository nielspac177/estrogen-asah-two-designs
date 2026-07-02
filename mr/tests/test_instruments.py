"""Distance clumping and position-based harmonization."""

import pandas as pd

from estrogen_mr.instruments import distance_clump, f_statistic, harmonize_by_position


def test_distance_clump_keeps_lowest_p_per_window():
    df = pd.DataFrame({
        "CHR": [1, 1, 1, 2],
        "POS": [1_000_000, 1_500_000, 5_000_000, 1_000_000],
        "Pval": [1e-9, 1e-12, 1e-8, 1e-10],
    })
    kept = distance_clump(df, window=1_000_000)
    # chr1: 1.0M & 1.5M within 1Mb -> keep lower-p (1.5M); 5.0M separate; chr2 separate
    assert set(kept["POS"]) == {1_500_000, 5_000_000, 1_000_000}
    assert 1_000_000 in kept[kept["CHR"] == 2]["POS"].values
    assert 1_000_000 not in kept[kept["CHR"] == 1]["POS"].values


def test_harmonize_by_position_flips_and_drops():
    inst = pd.DataFrame({
        "key": ["1:100", "1:200", "1:300"],
        "Effect_Allele": ["A", "A", "C"], "Other_Allele": ["G", "G", "G"],
        "EAF": [0.3, 0.3, 0.5], "Effect": [0.1, 0.2, 0.3], "SE": [0.01, 0.01, 0.01],
    })
    out = pd.DataFrame({
        "key": ["1:100", "1:200", "1:300"],
        # 200 swapped; 300 palindromic C/G @ EAF .5
        "oa_eff": ["A", "G", "C"], "oa_non": ["G", "A", "G"],
        "beta_out": [0.5, 0.5, 0.5], "se_out": [0.02, 0.02, 0.02],
    })
    h = harmonize_by_position(inst, out).set_index("key")
    assert h.loc["1:100", "by"] == 0.5      # aligned
    assert h.loc["1:200", "by"] == -0.5     # swapped -> flipped
    assert "1:300" not in h.index           # ambiguous palindrome dropped


def test_f_statistic():
    assert f_statistic([0.1, 0.2], [0.01, 0.01]) > 100
