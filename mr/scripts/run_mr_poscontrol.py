"""Positive control: age at natural menopause -> breast cancer, SAME pipeline.

Later menopause is an established risk factor for breast cancer (~OR 1.05 per
year). Recovering this with our instruments/harmonisation/estimators validates
the build, matching, and power, and makes the aSAH null interpretable. If it
fails, the aSAH estimate is an artifact. (Michailidou 2017, GCST004988.)
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from estrogen_mr.instruments import distance_clump, f_statistic, harmonize_by_position
from estrogen_mr.methods import cochran_q, ivw, mr_egger, weighted_median

ROOT = Path(__file__).resolve().parents[1]
GWAS = ROOT / "data" / "gwas"
OUT = ROOT / "outputs"
BC = GWAS / "bc_build37.tsv.gz"   # GRCh37 build, matches ANM/Bakker


def main() -> None:
    OUT.mkdir(exist_ok=True)
    exp = pd.read_csv(GWAS / "anm_ruth2021.txt.gz", sep="\t")
    exp["Pval"] = pd.to_numeric(exp["Pval"], errors="coerce")
    exp = exp[exp["Pval"] < 5e-8].dropna(subset=["Pval", "Effect", "SE"])
    inst = distance_clump(exp).copy()
    inst["key"] = inst["CHR"].astype(str) + ":" + inst["POS"].astype(str)
    keys = set(inst["key"])

    # Use the ORIGINAL (author-submitted, GRCh37) columns — NOT the hm_* harmonised
    # columns, which are lifted to GRCh38 and would not match GRCh37 instruments.
    cols = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error"]
    hits = []
    for ch in pd.read_csv(BC, sep="\t", usecols=cols, chunksize=1_000_000, low_memory=False):
        ch = ch.dropna(subset=["chromosome", "base_pair_location", "beta", "standard_error"])
        chrom = pd.to_numeric(ch["chromosome"], errors="coerce")
        ch = ch[chrom.notna()]
        ch["key"] = (chrom[chrom.notna()].astype(int).astype(str) + ":"
                     + pd.to_numeric(ch["base_pair_location"]).astype(int).astype(str))
        hits.append(ch[ch["key"].isin(keys)])
    out = pd.concat(hits, ignore_index=True).rename(columns={
        "effect_allele": "oa_eff", "other_allele": "oa_non",
        "beta": "beta_out", "standard_error": "se_out"})

    h = harmonize_by_position(inst, out[["key", "oa_eff", "oa_non", "beta_out", "se_out"]])
    bx, sx, by, sy = (h[c].values for c in ("bx", "sx", "by", "sy"))

    egg, egg_i = mr_egger(bx, sx, by, sy)
    res = {
        "exposure": "age_at_natural_menopause_Ruth2021",
        "outcome": "breast_cancer_Michailidou2017 (76192 cases / 63082 controls)",
        "n_harmonized": int(len(h)), "mean_F": round(f_statistic(bx, sx), 1),
        "IVW": ivw(bx, sx, by, sy).as_dict(),
        "MR_Egger": egg.as_dict(), "egger_intercept": egg_i,
        "weighted_median": weighted_median(bx, sx, by, sy, n_boot=1000).as_dict(),
        "cochran_q": cochran_q(bx, sx, by, sy),
        "expected": "OR ~1.05 per year later menopause (established); direction/scale validation",
    }
    (OUT / "mr_poscontrol_breastcancer.json").write_text(json.dumps(res, indent=2))
    for m in ("IVW", "MR_Egger", "weighted_median"):
        r = res[m]
        print(f"  {m:16s} OR/yr={r['or']:.3f} ({r['or_ci_low']:.3f}-{r['or_ci_high']:.3f})  p={r['p']:.3g}")
    print(f"  n={len(h)}, mean F={res['mean_F']}, expected ~1.05/yr")


if __name__ == "__main__":
    main()
