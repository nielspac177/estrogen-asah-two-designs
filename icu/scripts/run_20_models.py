"""run_20 — primary + secondary models, source meta, E-value (P7.3-7.7)."""

from __future__ import annotations

import json

import pandas as pd

from estrogen_asah_dci.analysis import models
from estrogen_asah_dci.analysis.evalue import evalue_or
from estrogen_asah_dci.analysis.meta import pool_or
from estrogen_asah_dci.pipeline import OUTPUTS


def main() -> None:
    cohort = pd.read_parquet(OUTPUTS / "cohort.parquet")

    primary = models.primary_logistic(cohort)
    sex = models.sex_difference(cohort)
    hrt = models.exposure_logistic(cohort, "hrt_exposure")
    per_source = {k: v for k, v in models.by_source(cohort).items() if v.converged}

    results = {
        "primary_post_vs_pre": primary.as_dict(),
        "secondary_sex": sex.as_dict(),
        "exploratory_hrt": hrt.as_dict(),
        "by_source": {k: v.as_dict() for k, v in per_source.items()},
    }

    if len(per_source) >= 2:
        pooled = pool_or([v.as_dict() for v in per_source.values()], random=True)
        results["meta_random"] = pooled.__dict__

    if primary.converged:
        base = float(cohort["dci_composite"].astype("boolean").astype("float").mean())
        ev = evalue_or(primary.odds_ratio, base, primary.ci_low, primary.ci_high)
        results["evalue"] = {"point": ev.point, "ci_bound": ev.ci_bound, "baseline_risk": base}

    (OUTPUTS / "results.json").write_text(json.dumps(results, indent=2))
    print("[run_20] primary (post vs pre):",
          f"OR={primary.odds_ratio:.2f} ({primary.ci_low:.2f}-{primary.ci_high:.2f})"
          if primary.converged else "did not converge")
    if "evalue" in results:
        print(f"[run_20] E-value point={results['evalue']['point']:.2f}")
    print(f"[run_20] wrote {OUTPUTS/'results.json'}")


if __name__ == "__main__":
    main()
