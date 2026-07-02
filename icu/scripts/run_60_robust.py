"""run_60 — pre-specified robustness / competing-risk sweep (audit response).

Reports ALL results (no cherry-picking). See docs/adr/0004-robustness-sap.md.
"""

from __future__ import annotations

import json

import pandas as pd

from estrogen_asah_dci.analysis import models, robust
from estrogen_asah_dci.pipeline import OUTPUTS


def _row(name, est):
    return {"analysis": name, "or": est.odds_ratio, "ci_low": est.ci_low,
            "ci_high": est.ci_high, "p": est.p, "n": est.n,
            "events": est.events, "converged": est.converged}


def main() -> None:
    cohort = pd.read_parquet(OUTPUTS / "cohort.parquet")

    results = [
        _row("primary (logistic)", models.primary_logistic(cohort)),
        _row("survivor-restricted", robust.survivor_restricted(cohort)),
        _row("age-restricted 45-55", robust.age_restricted(cohort, 45, 55)),
        _row("overlap-weighted (ATO)", robust.overlap_weighted(cohort)),
        _row("outcome: vasospasm-only", robust.outcome_variant(cohort, "vasospasm_dx")),
        _row("outcome: DCI-procedure-only", robust.outcome_variant(cohort, "dci_procedure")),
        _row("outcome: delayed-infarction", robust.outcome_variant(cohort, "delayed_infarction")),
    ]
    mnl = robust.competing_events_multinomial(cohort)

    out = {"post_vs_pre_variants": results, "competing_events_multinomial": mnl}
    (OUTPUTS / "results_robust.json").write_text(json.dumps(out, indent=2))

    print(f"{'analysis':28s} {'OR (95% CI)':22s} {'p':>6s} {'n':>5s} {'events':>7s}")
    for r in results:
        ci = (f"{r['or']:.2f} ({r['ci_low']:.2f}-{r['ci_high']:.2f})"
              if r["converged"] else "did not converge")
        p = f"{r['p']:.3f}" if r["converged"] else "-"
        print(f"{r['analysis']:28s} {ci:22s} {p:>6s} {r['n']:>5d} {r['events']:>7d}")
    if mnl["converged"]:
        print(f"\ncompeting-events multinomial — postmenopausal RRR for DCI (vs neither, "
              f"death as separate competing outcome): "
              f"{mnl['rrr_dci']:.2f} ({mnl['ci_low']:.2f}-{mnl['ci_high']:.2f}), n={mnl['n']}")


if __name__ == "__main__":
    main()
