"""Allele harmonization: sign flips on swap, drops ambiguous/incompatible SNPs."""

import pandas as pd

from estrogen_mr.harmonize import harmonize


def _exp():
    return pd.DataFrame({
        "SNP": ["rs1", "rs2", "rs3", "rs4"],
        "effect_allele": ["A", "A", "A", "A"],
        "other_allele": ["G", "G", "T", "G"],
        "beta": [0.1, 0.2, 0.3, 0.4],
        "se": [0.01, 0.01, 0.01, 0.01],
        "eaf": [0.3, 0.3, 0.5, 0.3],
    })


def test_swap_flips_sign_and_alignment():
    out = pd.DataFrame({
        "SNP": ["rs1", "rs2"],
        "effect_allele": ["A", "G"],   # rs1 aligned, rs2 swapped
        "other_allele": ["G", "A"],
        "beta": [0.5, 0.5],
        "se": [0.02, 0.02],
    })
    h = harmonize(_exp(), out).set_index("SNP")
    assert h.loc["rs1", "beta_out"] == 0.5      # aligned, unchanged
    assert h.loc["rs2", "beta_out"] == -0.5     # swapped -> sign flipped


def test_palindromic_ambiguous_dropped():
    # rs3 is A/T (palindromic) with eaf 0.5 -> ambiguous -> dropped
    out = pd.DataFrame({
        "SNP": ["rs3"], "effect_allele": ["A"], "other_allele": ["T"],
        "beta": [0.3], "se": [0.02],
    })
    assert harmonize(_exp(), out).empty


def test_incompatible_alleles_dropped():
    # C/A is neither the same nor the strand-complement of rs1's A/G -> incompatible
    out = pd.DataFrame({
        "SNP": ["rs1"], "effect_allele": ["C"], "other_allele": ["A"],
        "beta": [0.3], "se": [0.02],
    })
    assert harmonize(_exp(), out).empty


def test_strand_flip_kept():
    # C/T IS the strand-complement of A/G (same variant, opposite strand) -> kept
    out = pd.DataFrame({
        "SNP": ["rs1"], "effect_allele": ["C"], "other_allele": ["T"],
        "beta": [0.3], "se": [0.02],
    })
    assert len(harmonize(_exp(), out)) == 1
