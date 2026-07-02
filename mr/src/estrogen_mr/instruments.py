"""Instrument selection and position-based harmonization for two-sample MR.

Note on clumping: proper LD clumping needs a reference panel (e.g. 1000G) + PLINK.
When unavailable we use *distance* clumping (keep the most significant SNP, drop
others within ±window on the same chromosome) — a conservative approximation that
tends to retain independent signals but should be replaced by r²-clumping for a
final analysis. This limitation is reported alongside results.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

_PALINDROMIC = ({"A", "T"}, {"C", "G"})


def distance_clump(df: pd.DataFrame, pos="POS", chrom="CHR", pval="Pval",
                   window: int = 1_000_000) -> pd.DataFrame:
    """Greedy distance clumping: lowest-p SNP kept, neighbours within ±window dropped."""
    keep = []
    for _, g in df.sort_values(pval).groupby(chrom):
        chosen: list[int] = []
        for idx, r in g.iterrows():
            if all(abs(r[pos] - p) > window for p in chosen):
                chosen.append(r[pos])
                keep.append(idx)
    return df.loc[keep]


def harmonize_by_position(inst: pd.DataFrame, out: pd.DataFrame,
                          eaf_ambiguous: float = 0.42) -> pd.DataFrame:
    """Align outcome effects to exposure effect allele, keyed on CHR:POS.

    ``inst`` needs Effect_Allele/Other_Allele/EAF/Effect/SE + a ``key`` column;
    ``out`` needs oa_eff/oa_non/beta_out/se_out + a ``key`` column.
    Returns per-SNP bx/sx/by/sy for the MR estimators.
    """
    m = inst.merge(out, on="key")
    rows = []
    for _, r in m.iterrows():
        ea, oa = str(r["Effect_Allele"]).upper(), str(r["Other_Allele"]).upper()
        oe, on_ = str(r["oa_eff"]).upper(), str(r["oa_non"]).upper()
        b = r["beta_out"]
        if {ea, oa} != {oe, on_}:
            continue  # position collision / incompatible alleles
        if (oe, on_) == (oa, ea):
            b = -b  # allele swap
        if {ea, oa} in _PALINDROMIC:
            eaf = float(r.get("EAF", np.nan))
            if not np.isnan(eaf) and abs(eaf - 0.5) < (0.5 - eaf_ambiguous):
                continue  # ambiguous palindrome
        rows.append({"key": r["key"], "bx": r["Effect"], "sx": r["SE"],
                     "by": b, "sy": r["se_out"]})
    return pd.DataFrame(rows)


def f_statistic(bx, sx) -> float:
    bx, sx = np.asarray(bx), np.asarray(sx)
    return float(np.mean((bx / sx) ** 2))
