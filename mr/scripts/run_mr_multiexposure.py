"""Multi-exposure MR: sex-hormone exposures -> aSAH, one pipeline, one forest.

Adds SHBG (women, Ruth 2020 GCST90012107) and total testosterone
(GCST90012113) alongside age at natural menopause, all against the same Bakker
2020 aSAH outcome. The SHBG-in-women estimate directly tests the published
direction conflict (Molenberg 2022: higher SHBG -> more aSAH, OR 1.18 in women;
Tan/Wu 2025: higher SHBG -> less). Distance clumping (LD-clump approximation).
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from estrogen_mr import qc
from estrogen_mr.instruments import distance_clump, f_statistic, harmonize_by_position
from estrogen_mr.methods import ivw, mr_egger, weighted_median

ROOT = Path(__file__).resolve().parents[1]
GWAS = ROOT / "data" / "gwas"
OUT = ROOT / "outputs"
N_OUT = 17019

# name -> (file, is_gwas_catalog_format, n_exp, prior_note)
EXPOSURES = {
    "age_at_menopause": ("anm_ruth2021.txt.gz", False, 201323, "ANM (Ruth 2021)"),
    "SHBG_women": ("shbg_female.tsv.gz", True, 189473, "SHBG women (Ruth 2020)"),
    "total_testosterone": ("testosterone_total.tsv.gz", True, 194453, "Total T (Ruth 2020)"),
}


def load_sig_instruments(path: Path, catalog: bool) -> pd.DataFrame:
    """awk-prefilter to genome-wide-significant rows, then load + normalize columns."""
    if catalog:  # chromosome,base_pair_location,effect_allele,other_allele,eaf,beta,se,p_value
        pcol, ren = 9, {"chromosome": "CHR", "base_pair_location": "POS",
                        "effect_allele": "Effect_Allele", "other_allele": "Other_Allele",
                        "effect_allele_frequency": "EAF", "beta": "Effect", "standard_error": "SE"}
    else:  # ANM: SNP,CHR,POS,Effect_Allele,Other_Allele,EAF,Effect,SE,Pval,N
        pcol, ren = 9, {}
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf:
        subprocess.run(["bash", "-c", f"gzcat '{path}' | awk 'NR==1 || (${pcol}!=\"NA\" && ${pcol}+0<5e-8)'"],
                       stdout=tf, check=True)
        tmp = tf.name
    df = pd.read_csv(tmp, sep="\t").rename(columns=ren)
    Path(tmp).unlink(missing_ok=True)
    df["Pval"] = pd.to_numeric(df.get("Pval", df.get("p_value")), errors="coerce")
    df = df.dropna(subset=["Effect", "SE", "CHR", "POS"])
    inst = distance_clump(df).copy()
    inst["key"] = inst["CHR"].astype(int).astype(str) + ":" + inst["POS"].astype(int).astype(str)
    return inst


def outcome_frame() -> pd.DataFrame:
    out = pd.read_csv(GWAS / "bakker_sah_euro_noUKBB.txt.gz", sep="\t")
    out["key"] = out["CHR"].astype(str) + ":" + out["BP"].astype(str)
    return out.rename(columns={"A_EFF": "oa_eff", "A_NONEFF": "oa_non",
                               "BETA": "beta_out", "SE": "se_out"})


def main() -> None:
    OUT.mkdir(exist_ok=True)
    out = outcome_frame()[["key", "oa_eff", "oa_non", "beta_out", "se_out"]]
    rows = []
    for name, (fname, catalog, n_exp, note) in EXPOSURES.items():
        inst = load_sig_instruments(GWAS / fname, catalog)
        h = harmonize_by_position(inst, out)
        if len(h) < 5:
            rows.append({"exposure": name, "note": note, "n": len(h), "converged": False})
            continue
        bx, sx, by, sy = (h[c].values for c in ("bx", "sx", "by", "sy"))
        re = qc.ivw_random(bx, sx, by, sy)
        egg, egg_i = mr_egger(bx, sx, by, sy)
        wm = weighted_median(bx, sx, by, sy, n_boot=500)
        st = qc.steiger(bx, sx, by, sy, n_exp, N_OUT)
        loo = [x["or"] for x in qc.leave_one_out(bx, sx, by, sy)]
        rows.append({
            "exposure": name, "note": note, "n_clumped": int(len(inst)), "n": int(len(h)),
            "mean_F": round(f_statistic(bx, sx), 1),
            "IVW_random": re.as_dict(), "IVW_fixed": ivw(bx, sx, by, sy).as_dict(),
            "MR_Egger": egg.as_dict(), "egger_intercept_p": egg_i["intercept_p"],
            "weighted_median": wm.as_dict(),
            "steiger_frac_correct": round(st["frac_correct"], 2),
            "loo_or_range": [round(min(loo), 3), round(max(loo), 3)], "converged": True,
        })
    (OUT / "mr_multiexposure_aSAH.json").write_text(json.dumps(rows, indent=2))

    print(f"{'exposure':22s} {'n':>4s} {'F':>6s}  {'IVW-RE OR (95% CI)':22s} {'p':>7s}  Egger_p  Steiger")
    for r in rows:
        if not r["converged"]:
            print(f"{r['exposure']:22s} n={r['n']} (too few instruments)")
            continue
        v = r["IVW_random"]
        ci = f"{v['or']:.2f} ({v['or_ci_low']:.2f}-{v['or_ci_high']:.2f})"
        print(f"{r['exposure']:22s} {r['n']:>4d} {r['mean_F']:>6.0f}  {ci:22s} {v['p']:>7.3g}  "
              f"{r['egger_intercept_p']:.2f}     {r['steiger_frac_correct']}")


if __name__ == "__main__":
    main()
