"""Primary MR: age at natural menopause (Ruth 2021) -> aneurysmal SAH (Bakker 2020).

The untested gap (ADR-0001). Requires the GWAS in data/gwas/ (see sources.py);
they are git-ignored. Distance-clumping is an LD-clumping approximation (see
instruments.py) — flagged in outputs. Genome build assumed GRCh37/hg19 for both.
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


def main() -> None:
    OUT.mkdir(exist_ok=True)
    exp = pd.read_csv(GWAS / "anm_ruth2021.txt.gz", sep="\t")
    exp["Pval"] = pd.to_numeric(exp["Pval"], errors="coerce")
    exp = exp[exp["Pval"] < 5e-8].dropna(subset=["Pval", "Effect", "SE"])
    inst = distance_clump(exp).copy()
    inst["key"] = inst["CHR"].astype(str) + ":" + inst["POS"].astype(str)

    out = pd.read_csv(GWAS / "bakker_sah_euro_noUKBB.txt.gz", sep="\t")
    out["key"] = out["CHR"].astype(str) + ":" + out["BP"].astype(str)
    out = out.rename(columns={"A_EFF": "oa_eff", "A_NONEFF": "oa_non",
                              "BETA": "beta_out", "SE": "se_out"})

    h = harmonize_by_position(inst, out[["key", "oa_eff", "oa_non", "beta_out", "se_out"]])
    bx, sx, by, sy = h["bx"].values, h["sx"].values, h["by"].values, h["sy"].values

    egger, egg_info = mr_egger(bx, sx, by, sy)
    q = cochran_q(bx, sx, by, sy)
    results = {
        "exposure": "age_at_natural_menopause_Ruth2021",
        "outcome": "aSAH_European_UKBexcluded_Bakker2020",
        "n_instruments_clumped": int(len(inst)),
        "n_harmonized": int(len(h)),
        "mean_F": round(f_statistic(bx, sx), 1),
        "IVW": ivw(bx, sx, by, sy).as_dict(),
        "MR_Egger": egger.as_dict(),
        "weighted_median": weighted_median(bx, sx, by, sy, n_boot=1000).as_dict(),
        "egger_intercept": egg_info,
        "cochran_q": q,
        "clumping": "distance +/-1Mb (LD-clump approximation; see instruments.py)",
        "note": "OR per 1 year later menopause. OR<1 => more lifetime estrogen lowers aSAH risk.",
    }
    (OUT / "mr_menopause_aSAH.json").write_text(json.dumps(results, indent=2))

    print(f"instruments (clumped): {results['n_instruments_clumped']}, "
          f"harmonized: {results['n_harmonized']}, mean F={results['mean_F']}")
    for m in ("IVW", "MR_Egger", "weighted_median"):
        r = results[m]
        ci = f"({r['or_ci_low']:.3f}-{r['or_ci_high']:.3f})"
        print(f"  {m:16s} OR={r['or']:.3f} {ci}  p={r['p']:.3g}")
    print(f"  Egger intercept={egg_info['intercept']:.4f} (p={egg_info['intercept_p']:.3g}); "
          f"Q p={q['p']:.3g}, I2={q['i2']:.2f}")


if __name__ == "__main__":
    main()
