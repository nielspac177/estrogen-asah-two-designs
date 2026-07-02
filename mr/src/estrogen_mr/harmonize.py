"""Harmonize exposure and outcome GWAS to a common effect allele.

Aligns the outcome effect estimate to the exposure's effect allele (flipping the
sign when alleles are swapped) and drops strand-ambiguous palindromic SNPs with
intermediate frequency, which cannot be aligned reliably.
"""

from __future__ import annotations

import pandas as pd

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
_PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}

REQUIRED = ["SNP", "effect_allele", "other_allele", "beta", "se"]


def _palindromic(a1: str, a2: str) -> bool:
    return frozenset((a1.upper(), a2.upper())) in _PALINDROMIC


def harmonize(exp: pd.DataFrame, out: pd.DataFrame,
              eaf_ambiguous: float = 0.42) -> pd.DataFrame:
    """Merge on SNP; return columns beta_exp/se_exp/beta_out/se_out aligned to exposure allele.

    ``exp`` must carry ``eaf`` (effect-allele frequency) to screen palindromes.
    """
    for name, df in (("exposure", exp), ("outcome", out)):
        missing = set(REQUIRED) - set(df.columns)
        if missing:
            raise ValueError(f"{name} missing columns: {sorted(missing)}")

    m = exp.merge(out, on="SNP", suffixes=("_exp", "_out"))
    rows = []
    for _, r in m.iterrows():
        e1, e2 = r["effect_allele_exp"].upper(), r["other_allele_exp"].upper()
        o1, o2 = r["effect_allele_out"].upper(), r["other_allele_out"].upper()
        beta_out, se_out = r["beta_out"], r["se_out"]

        # same alleles, possibly swapped or on the complementary strand
        if {e1, e2} == {o1, o2}:
            if (o1, o2) == (e2, e1):
                beta_out = -beta_out          # allele swap -> flip sign
        elif {_COMPLEMENT[o1], _COMPLEMENT[o2]} == {e1, e2}:
            if (_COMPLEMENT[o1], _COMPLEMENT[o2]) == (e2, e1):
                beta_out = -beta_out          # strand flip + swap
        else:
            continue  # incompatible alleles -> drop

        # drop ambiguous palindromes at intermediate frequency
        eaf = r.get("eaf", r.get("eaf_exp"))
        ambiguous = eaf is not None and abs(float(eaf) - 0.5) < (0.5 - eaf_ambiguous)
        if _palindromic(e1, e2) and ambiguous:
            continue

        rows.append({"SNP": r["SNP"], "beta_exp": r["beta_exp"], "se_exp": r["se_exp"],
                     "beta_out": beta_out, "se_out": se_out})
    return pd.DataFrame(rows)
