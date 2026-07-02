"""run_70 — specification curve + age×sex spline DiD (audit-endorsed analyses)."""

from __future__ import annotations

import json

import pandas as pd

from estrogen_asah_dci.analysis.multiverse import age_sex_did, spec_curve, spec_curve_summary
from estrogen_asah_dci.pipeline import OUTPUTS


def main() -> None:
    cohort = pd.read_parquet(OUTPUTS / "cohort.parquet")

    curve = spec_curve(cohort)
    curve.to_csv(OUTPUTS / "spec_curve.csv", index=False)
    summary = spec_curve_summary(curve)
    did = age_sex_did(cohort)

    (OUTPUTS / "results_multiverse.json").write_text(
        json.dumps({"spec_curve_summary": summary, "age_sex_did": did}, indent=2)
    )

    print("[run_70] SPECIFICATION CURVE")
    for k, v in summary.items():
        print(f"    {k}: {v}")
    print("\n[run_70] AGE×SEX SPLINE DiD (female-minus-male DCI age effect, 55 vs 45)")
    if did.get("converged"):
        print(f"    OR-ratio = {did['or_ratio']:.2f} ({did['ci_low']:.2f}-{did['ci_high']:.2f}), "
              f"n={did['n']}, boots={did['n_boot']}")
        print(f"    {did['interpretation']}")
    else:
        print("    did not converge:", did.get("error"))


if __name__ == "__main__":
    main()
