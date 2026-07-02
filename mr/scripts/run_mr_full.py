"""Full-QC MR: age at natural menopause -> aSAH, with sensitivity, power, equivalence.

Implements the adversarial-audit checklist (docs/audit/adversarial_synthesis.md):
random-effects IVW, Steiger, leave-one-out, MR-PRESSO global, power/MDE, TOST.
Reports on the per-SD scale and frames the result as a bound, not "no effect".
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from estrogen_mr import qc
from estrogen_mr.instruments import distance_clump, f_statistic, harmonize_by_position
from estrogen_mr.methods import ivw, mr_egger, weighted_median

ROOT = Path(__file__).resolve().parents[1]
GWAS = ROOT / "data" / "gwas"
OUT = ROOT / "outputs"
N_EXP = 201323      # Ruth 2021 ANM
N_OUT = 17019       # Bakker aSAH-Euro-noUKBB, max effective N
ANM_SD_YEARS = 4.0  # SD of age at natural menopause (sensitivity 3.5-4.5)


def harmonized():
    exp = pd.read_csv(GWAS / "anm_ruth2021.txt.gz", sep="\t")
    exp["Pval"] = pd.to_numeric(exp["Pval"], errors="coerce")
    exp = exp[exp["Pval"] < 5e-8].dropna(subset=["Pval", "Effect", "SE"])
    inst = distance_clump(exp).copy()
    inst["key"] = inst["CHR"].astype(str) + ":" + inst["POS"].astype(str)
    out = pd.read_csv(GWAS / "bakker_sah_euro_noUKBB.txt.gz", sep="\t")
    out["key"] = out["CHR"].astype(str) + ":" + out["BP"].astype(str)
    out = out.rename(columns={"A_EFF": "oa_eff", "A_NONEFF": "oa_non",
                              "BETA": "beta_out", "SE": "se_out"})
    return inst, harmonize_by_position(inst, out[["key", "oa_eff", "oa_non", "beta_out", "se_out"]])


def main() -> None:
    OUT.mkdir(exist_ok=True)
    inst, h = harmonized()
    bx, sx, by, sy = (h[c].values for c in ("bx", "sx", "by", "sy"))

    ivw_fe = ivw(bx, sx, by, sy)
    ivw_re = qc.ivw_random(bx, sx, by, sy)
    egg, egg_i = mr_egger(bx, sx, by, sy)
    wm = weighted_median(bx, sx, by, sy, n_boot=1000)
    st = qc.steiger(bx, sx, by, sy, N_EXP, N_OUT)
    presso = qc.mr_presso_global(bx, sx, by, sy)
    loo = qc.leave_one_out(bx, sx, by, sy)
    pw = qc.power_mde(ivw_re.se, sd_units=ANM_SD_YEARS)
    # equivalence test on the per-SD scale (SESOI is defined per SD)
    eq = qc.tost(ivw_re.estimate * ANM_SD_YEARS, ivw_re.se * ANM_SD_YEARS, sesoi_or=0.90)

    loo_or = [x["or"] for x in loo]
    results = {
        "exposure": "age_at_natural_menopause_Ruth2021",
        "outcome": "aSAH_Euro_UKBexcluded_Bakker2020 (Neff<=17019)",
        "n_instruments_clumped": int(len(inst)), "n_harmonized": int(len(h)),
        "mean_F": round(f_statistic(bx, sx), 1),
        "primary_IVW_random": ivw_re.as_dict(),
        "IVW_fixed": ivw_fe.as_dict(), "MR_Egger": egg.as_dict(),
        "egger_intercept": egg_i, "weighted_median": wm.as_dict(),
        "steiger": {k: v for k, v in st.items() if k != "keep_mask"},
        "mr_presso_global": presso,
        "leave_one_out_or_range": [round(min(loo_or), 3), round(max(loo_or), 3)],
        "power": pw, "equivalence_TOST": eq,
        "clumping": "distance +/-1Mb (LD-clump approximation)",
    }
    (OUT / "mr_menopause_aSAH_full.json").write_text(json.dumps(results, indent=2))

    r = ivw_re
    print(f"instruments {len(h)}  mean F={results['mean_F']}  "
          f"Steiger correct {st['n_correct']}/{st['n']}")
    print(f"PRIMARY IVW (random effects): OR/yr {r.odds_ratio:.3f} "
          f"({r.ci_low and __import__('math').exp(r.ci_low):.3f}-"
          f"{__import__('math').exp(r.ci_high):.3f}), p={r.p:.3g}")
    print(f"  per SD ({ANM_SD_YEARS}y): OR {pw['mde_or_perSD_protect']:.2f} MDE(protect) / "
          f"{pw['mde_or_perSD_harm']:.2f} MDE(harm) at 80% power")
    print(f"  power per-SD: {pw['power_perSD_at']}")
    print(f"  MR-Egger intercept p={egg_i['intercept_p']:.3g}; PRESSO global p={presso['p_global']:.3g}")
    print(f"  leave-one-out OR range: {results['leave_one_out_or_range']}")
    print(f"  TOST vs SESOI OR 0.90/SD: reject strong protection p="
          f"{eq['p_reject_strong_protection']:.3g}; reject strong harm p={eq['p_reject_strong_harm']:.3g}")


if __name__ == "__main__":
    main()
