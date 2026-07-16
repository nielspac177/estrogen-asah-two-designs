"""Formation-vs-rupture MR: age at natural menopause (Ruth 2021) against the
rupture-stratified Bakker 2020 outcomes (aSAH, unruptured IA, IA-combined).

Reuses the exact validated pipeline (distance-clumped ANM instruments, IVW random
effects, MR-Egger, weighted median, Steiger, MR-PRESSO, leave-one-out, power/MDE,
TOST) applied to each outcome so the three disease stages are directly comparable.
Outcome effective N is read per-SNP from the Bakker N column (max over instruments).

The unruptured-IA file is the only Bakker uIA stratum and includes UK Biobank, so
it overlaps the UKB-derived ANM exposure. With strong instruments (mean F~98) the
overlap bias is small, and it biases toward the confounded observational estimate
(which leaned protective, OR<1); a null/non-protective uIA estimate is therefore
conservative. This is flagged in the output.
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
N_EXP = 201323       # Ruth 2021 ANM
ANM_SD_YEARS = 4.0   # SD of age at natural menopause

OUTCOMES = [
    {"key": "aSAH_euro_noUKBB", "file": "bakker_sah_euro_noUKBB.txt.gz",
     "label": "Aneurysmal SAH (ruptured), European, UKB-excluded", "overlap": False},
    {"key": "uIA_euro", "file": "bakker_uia_euro.txt.gz",
     "label": "Unruptured IA (formation), European", "overlap": True},
    {"key": "IA_combined_euro_noUKBB", "file": "bakker_ia_stage1_noUKBB.txt.gz",
     "label": "IA combined (ruptured+unruptured), European, UKB-excluded", "overlap": False},
]


def anm_instruments() -> pd.DataFrame:
    exp = pd.read_csv(GWAS / "anm_ruth2021.txt.gz", sep="\t")
    exp["Pval"] = pd.to_numeric(exp["Pval"], errors="coerce")
    exp = exp[exp["Pval"] < 5e-8].dropna(subset=["Pval", "Effect", "SE"])
    inst = distance_clump(exp).copy()
    inst["key"] = inst["CHR"].astype(str) + ":" + inst["POS"].astype(str)
    return inst


def load_outcome(fname: str) -> pd.DataFrame:
    out = pd.read_csv(GWAS / fname, sep="\t")
    out["key"] = out["CHR"].astype(str) + ":" + out["BP"].astype(str)
    return out.rename(columns={"A_EFF": "oa_eff", "A_NONEFF": "oa_non",
                               "BETA": "beta_out", "SE": "se_out", "Neff": "N"})


def run_one(inst: pd.DataFrame, spec: dict) -> dict:
    out = load_outcome(spec["file"])
    cols = ["key", "oa_eff", "oa_non", "beta_out", "se_out", "N"]
    h = harmonize_by_position(inst, out[cols[:-1]])
    # per-outcome effective N: max Neff over the harmonized instrument positions
    n_out = int(out.merge(h[["key"]], on="key")["N"].max())
    bx, sx, by, sy = (h[c].values for c in ("bx", "sx", "by", "sy"))

    ivw_re = qc.ivw_random(bx, sx, by, sy)
    egg, egg_i = mr_egger(bx, sx, by, sy)
    wm = weighted_median(bx, sx, by, sy, n_boot=1000)
    st = qc.steiger(bx, sx, by, sy, N_EXP, n_out)
    presso = qc.mr_presso_global(bx, sx, by, sy)
    loo = [x["or"] for x in qc.leave_one_out(bx, sx, by, sy)]
    pw = qc.power_mde(ivw_re.se, sd_units=ANM_SD_YEARS)
    eq = qc.tost(ivw_re.estimate * ANM_SD_YEARS, ivw_re.se * ANM_SD_YEARS, sesoi_or=0.90)

    return {
        "outcome_key": spec["key"], "outcome_label": spec["label"],
        "exposure": "age_at_natural_menopause_Ruth2021",
        "sample_overlap_with_exposure": spec["overlap"],
        "n_harmonized": int(len(h)), "outcome_Neff_max": n_out,
        "mean_F": round(f_statistic(bx, sx), 1),
        "primary_IVW_random": ivw_re.as_dict(),
        "IVW_fixed": ivw(bx, sx, by, sy).as_dict(),
        "MR_Egger": egg.as_dict(), "egger_intercept": egg_i,
        "weighted_median": wm.as_dict(),
        "steiger": {k: v for k, v in st.items() if k != "keep_mask"},
        "mr_presso_global": presso,
        "leave_one_out_or_range": [round(min(loo), 3), round(max(loo), 3)],
        "power": pw, "equivalence_TOST": eq,
        "clumping": "distance +/-1Mb (LD-clump approximation)",
    }


def main() -> None:
    OUT.mkdir(exist_ok=True)
    inst = anm_instruments()
    print(f"ANM instruments (clumped): {len(inst)}")
    results = {}
    for spec in OUTCOMES:
        r = run_one(inst, spec)
        results[spec["key"]] = r
        ci = r["primary_IVW_random"]
        print(f"\n{spec['label']}")
        print(f"  harmonized {r['n_harmonized']}  Neff(max) {r['outcome_Neff_max']}  "
              f"F={r['mean_F']}  overlap={spec['overlap']}")
        print(f"  IVW-RE  OR/yr {ci['or']:.3f} ({ci['or_ci_low']:.3f}-{ci['or_ci_high']:.3f}) "
              f"p={ci['p']:.3g}")
        print(f"  Egger {r['MR_Egger']['or']:.3f} (int p={r['egger_intercept']['intercept_p']:.3g})  "
              f"WM {r['weighted_median']['or']:.3f}  PRESSO p={r['mr_presso_global']['p_global']:.3g}")
        print(f"  LOO {r['leave_one_out_or_range']}  Steiger {r['steiger']['n_correct']}/{r['steiger']['n']}")
        print(f"  TOST reject strong protection p={r['equivalence_TOST']['p_reject_strong_protection']:.3g}")
    (OUT / "mr_outcomes_stratified.json").write_text(json.dumps(results, indent=2))
    print("\nsaved outputs/mr_outcomes_stratified.json")


if __name__ == "__main__":
    main()
