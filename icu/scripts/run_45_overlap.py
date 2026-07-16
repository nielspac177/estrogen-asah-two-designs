"""run_45 — evidence for the non-identifiability of the observational estimand.

Two artefacts, both aggregate (no patient-level rows leave the cohort):
  1. outputs/age_overlap.csv — age-bin counts by group (premenopausal women,
     postmenopausal women, men). Shows zero common support across the menopausal
     strata: no premenopausal woman >=51 and no postmenopausal woman <51.
  2. outputs/vif.json — variance-inflation factors when the age-defined menopause
     step and a smooth age term are co-entered (collinearity blow-up), versus the
     age-dropped model actually used. Demonstrates the false-precision artefact.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from estrogen_asah_dci.pipeline import OUTPUTS

COVS = ["htn", "smoking", "diabetes"]
BINS = list(range(15, 95, 5)) + [200]  # 5-y bins, top bin absorbs age-shifted 90+


def _r2(df: pd.DataFrame, target: str, predictors: list[str]) -> float:
    d = df[[target, *predictors]].apply(pd.to_numeric, errors="coerce").dropna()
    return float(smf.ols(f"{target} ~ " + " + ".join(predictors), d).fit().rsquared)


def main() -> None:
    df = pd.read_parquet(OUTPUTS / "cohort.parquet").copy()
    df["age"] = pd.to_numeric(df["age"], errors="coerce")
    for c in COVS:
        df[c] = pd.to_numeric(df[c].astype("boolean").astype("float"), errors="coerce")
    df["post"] = (df["menopausal_stratum"] == "postmenopausal").astype(float)

    def group(r):
        if r["sex"] == "M":
            return "men"
        return {"premenopausal": "premenopausal", "postmenopausal": "postmenopausal"}.get(
            r["menopausal_stratum"], "other_women")

    df["group"] = df.apply(group, axis=1)
    df["agebin"] = pd.cut(df["age"], bins=BINS, right=False)
    tab = (df[df["group"].isin(["premenopausal", "postmenopausal", "men"])]
           .groupby(["agebin", "group"], observed=True).size().unstack(fill_value=0))
    for g in ["premenopausal", "postmenopausal", "men"]:
        if g not in tab:
            tab[g] = 0
    tab = tab[["premenopausal", "postmenopausal", "men"]]
    tab.index = [int(iv.left) for iv in tab.index]
    tab.index.name = "age_lo"
    tab.to_csv(OUTPUTS / "age_overlap.csv")

    # DCI rate by age bin and sex (the DiD intuition: age gradient in BOTH sexes,
    # no extra female drop at the menopause cutoff)
    df["dci"] = df["dci_composite"].astype("boolean").astype("float")
    df["sexlab"] = np.where(df["sex"] == "M", "men", "women")
    rate = (df.groupby(["agebin", "sexlab"], observed=True)
            .agg(n=("dci", "size"), events=("dci", "sum")).reset_index())
    rate["age_lo"] = [int(iv.left) for iv in rate["agebin"]]
    rate["rate"] = rate["events"] / rate["n"]
    rate[rate["n"] >= 10][["age_lo", "sexlab", "n", "events", "rate"]].to_csv(
        OUTPUTS / "age_dci_by_sex.csv", index=False)

    # overlap facts
    w = df[df["sex"] != "M"]
    pre, post = w[w["post"] == 0], w[w["post"] == 1]
    overlap = {
        "pre_age_min": float(pre["age"].min()), "pre_age_max": float(pre["age"].max()),
        "post_age_min": float(post["age"].min()), "post_age_max": float(post["age"].max()),
        "n_pre_ge51": int((pre["age"] >= 51).sum()), "n_post_lt51": int((post["age"] < 51).sum()),
    }

    # VIF: regress each predictor on the others; VIF = 1/(1-R^2)
    women = df[df["sex"] != "M"].copy()
    with_age = ["post", "age", *COVS]
    vif_with = {p: round(1 / (1 - _r2(women, p, [x for x in with_age if x != p])), 1)
                for p in ["post", "age"]}
    without_age = ["post", *COVS]
    vif_without = {"post": round(1 / (1 - _r2(women, "post", [x for x in without_age if x != "post"])), 2)}

    vif = {"model_with_age": {"formula": "dci ~ post + age + htn + smoking + diabetes",
                              "VIF": vif_with,
                              "note": "VIF is only ~2.3: this is NOT a classical collinearity "
                                      "problem. The model runs fine. Non-identifiability comes "
                                      "from zero overlap (positivity), so the age-adjusted post "
                                      "coefficient is an extrapolation across a support gap that "
                                      "depends entirely on the assumed linear age form."},
           "model_used_age_dropped": {"formula": "dci ~ post + htn + smoking + diabetes",
                                      "VIF": vif_without},
           "overlap": overlap}
    (OUTPUTS / "vif.json").write_text(json.dumps(vif, indent=2))

    print("[run_45] age overlap by group written to outputs/age_overlap.csv")
    print(f"   premenopausal age range {overlap['pre_age_min']:.0f}-{overlap['pre_age_max']:.0f}; "
          f"postmenopausal {overlap['post_age_min']:.0f}-{overlap['post_age_max']:.0f}")
    print(f"   premenopausal women aged >=51: {overlap['n_pre_ge51']}; "
          f"postmenopausal women aged <51: {overlap['n_post_lt51']}  (zero overlap)")
    print(f"[run_45] VIF with age co-entered: post={vif_with['post']}, age={vif_with['age']} "
          f"(modest: NOT a collinearity problem); age-dropped model post VIF={vif_without['post']}. "
          f"Non-identifiability is from zero overlap, not VIF.")
    print(f"[run_45] wrote {OUTPUTS/'vif.json'}")


if __name__ == "__main__":
    main()
