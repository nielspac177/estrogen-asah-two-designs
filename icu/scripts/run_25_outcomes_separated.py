"""run_25 — delayed cerebral ischaemia and in-hospital mortality as SEPARATE
cause-specific outcomes (competing risks), not a composite.

The audit (ADR-0005) showed death-inclusive composites and survivor restrictions
either push the odds ratio toward >1 mechanically or open a Severity->Death<-DCI
collider. The honest alternative is to report each cause-specific outcome on its
own and read them jointly: DCI and death compete, and because postmenopausal
women die more (age), their DCI denominator is truncated. Estimates use the same
adjusted, hospital-cluster-robust model as the primary.
"""

from __future__ import annotations

import json

import pandas as pd

from estrogen_asah_dci.analysis import models
from estrogen_asah_dci.pipeline import OUTPUTS

OUTCOMES = {"dci": "dci_composite", "mortality": "died"}


def main() -> None:
    cohort = pd.read_parquet(OUTPUTS / "cohort.parquet")

    results = {"contrast_post_vs_pre": {}, "contrast_male_vs_female": {},
               "crude_rates_by_stratum": {}}
    for name, col in OUTCOMES.items():
        results["contrast_post_vs_pre"][name] = models.menopause_on(cohort, col).as_dict()
        results["contrast_male_vs_female"][name] = models.sex_on(cohort, col).as_dict()
        results["crude_rates_by_stratum"][name] = (
            models.rates_by_stratum(cohort, col).to_dict("records"))

    # competing-risk note: joint DCI/death pattern by stratum
    note = []
    for name in OUTCOMES:
        for r in results["crude_rates_by_stratum"][name]:
            note.append(f"{r['menopausal_stratum']} {name}: {r['rate']*100:.1f}% "
                        f"({int(r['events'])}/{int(r['n'])})")
    results["competing_risk_note"] = (
        "DCI and in-hospital death are competing events; postmenopausal women have "
        "higher mortality (age), which truncates their DCI observation window. Read "
        "the two cause-specific estimates jointly, not as independent effects.")

    (OUTPUTS / "results_outcomes_separated.json").write_text(json.dumps(results, indent=2))

    print("[run_25] Cause-specific outcomes (post vs pre, adjusted, cluster-robust):")
    for name in OUTCOMES:
        e = results["contrast_post_vs_pre"][name]
        tag = (f"OR {e['or']:.2f} ({e['ci_low']:.2f}-{e['ci_high']:.2f}) p={e['p']:.3g}"
               if e["converged"] else "did not converge")
        print(f"   {name:9s}: {tag}   events {e['events']}/{e['n']}")
    print("[run_25] Male vs female:")
    for name in OUTCOMES:
        e = results["contrast_male_vs_female"][name]
        print(f"   {name:9s}: OR {e['or']:.2f} ({e['ci_low']:.2f}-{e['ci_high']:.2f}) p={e['p']:.3g}")
    print("[run_25] Crude rates by stratum:")
    for line in note:
        print("   " + line)
    print(f"[run_25] wrote {OUTPUTS/'results_outcomes_separated.json'}")


if __name__ == "__main__":
    main()
