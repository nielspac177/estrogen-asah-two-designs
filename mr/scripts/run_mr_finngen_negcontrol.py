"""Negative-control outcome for the ANM MR: acute appendicitis (FinnGen R11
K11_APPENDACUT), an outcome with no plausible causal link to age at natural
menopause. Pairs with the breast-cancer positive control: the positive control
shows the ANM instrument DETECTS a true hormonal effect, the negative control
shows it does NOT spuriously associate where no effect is expected, bounding
residual pleiotropy/confounding. Harmonized by rsID (FinnGen is GRCh38).

Reuses the clumped instrument table written by run_mr_finngen.py
(outputs/anm_instruments_rsid.tsv) and the same streaming rsID extract method.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from estrogen_mr.instruments import f_statistic, harmonize_by_position
from estrogen_mr.methods import ivw, mr_egger, weighted_median

ROOT = Path(__file__).resolve().parents[1]
GWAS = ROOT / "data" / "gwas"
OUT = ROOT / "outputs"
EXTRACT = GWAS / "finngen_R11_K11_APPENDACUT_extract.tsv"


def main() -> None:
    inst = pd.read_csv(OUT / "anm_instruments_rsid.tsv", sep="\t")
    inst["key"] = inst["rsid"]
    df = pd.read_csv(EXTRACT, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    df = df.rename(columns={"rsids": "key", "alt": "oa_eff", "ref": "oa_non",
                            "beta": "beta_out", "sebeta": "se_out"})
    out = df[df["key"].isin(set(inst["rsid"]))][["key", "oa_eff", "oa_non", "beta_out", "se_out"]]
    h = harmonize_by_position(inst, out)
    bx, sx, by, sy = (h[c].values for c in ("bx", "sx", "by", "sy"))

    ivw_re = ivw(bx, sx, by, sy)
    egg, egg_i = mr_egger(bx, sx, by, sy)
    wm = weighted_median(bx, sx, by, sy, n_boot=1000)
    res = {
        "outcome_label": "Acute appendicitis (FinnGen R11 K11_APPENDACUT) [negative control]",
        "exposure": "age_at_natural_menopause_Ruth2021",
        "rationale": "no plausible causal effect of menopausal timing; complements the "
                     "breast-cancer positive control to bound pleiotropy/confounding",
        "n_harmonized": int(len(h)), "mean_F": round(f_statistic(bx, sx), 1),
        "IVW": ivw_re.as_dict(), "MR_Egger": egg.as_dict(), "egger_intercept": egg_i,
        "weighted_median": wm.as_dict(),
    }
    (OUT / "mr_negcontrol_appendicitis.json").write_text(json.dumps(res, indent=2))
    v = ivw_re.as_dict()
    print(f"Negative control (appendicitis): harmonized {len(h)}, F={res['mean_F']}")
    print(f"  IVW OR/yr {v['or']:.3f} ({v['or_ci_low']:.3f}-{v['or_ci_high']:.3f}) p={v['p']:.3g}")
    print(f"  Egger {res['MR_Egger']['or']:.3f} (int p={egg_i['intercept_p']:.3g})  "
          f"WM {res['weighted_median']['or']:.3f}")
    print("  Expectation: null (OR~1). saved outputs/mr_negcontrol_appendicitis.json")


if __name__ == "__main__":
    main()
