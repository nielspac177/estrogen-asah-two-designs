"""run_40 — render publication figures from analysis outputs (P8)."""

from __future__ import annotations

import json

import pandas as pd

from estrogen_asah_dci.pipeline import OUTPUTS
from estrogen_asah_dci.viz import figures


def forest_rows(results: dict) -> list[dict]:
    rows = []
    mapping = [
        ("Postmenopausal vs premenopausal", "primary_post_vs_pre"),
        ("Male vs female", "secondary_sex"),
        ("HRT exposure", "exploratory_hrt"),
    ]
    for label, key in mapping:
        r = results.get(key)
        if r and r.get("converged"):
            rows.append({"label": label, "or": r["or"],
                         "ci_low": r["ci_low"], "ci_high": r["ci_high"]})
    for src, r in results.get("by_source", {}).items():
        if r.get("converged"):
            rows.append({"label": f"  post vs pre — {src}", "or": r["or"],
                         "ci_low": r["ci_low"], "ci_high": r["ci_high"]})
    m = results.get("meta_random")
    if m:
        rows.append({"label": "Pooled (random effects)", "or": m["estimate"],
                     "ci_low": m["ci_low"], "ci_high": m["ci_high"]})
    return rows


def main() -> None:
    fig_dir = OUTPUTS / "figures"
    cohort = pd.read_parquet(OUTPUTS / "cohort.parquet")
    results = json.loads((OUTPUTS / "results.json").read_text())
    rates = pd.read_csv(OUTPUTS / "crude_rates.csv")
    balance = pd.read_csv(OUTPUTS / "balance.csv")

    women = cohort[cohort["sex"] == "F"]
    steps = [
        ("aSAH cohort", len(cohort)),
        ("Women", len(women)),
        ("Premenopausal (<51)", int((women["menopausal_stratum"] == "premenopausal").sum())),
        ("Postmenopausal (≥51)", int((women["menopausal_stratum"] == "postmenopausal").sum())),
    ]

    paths = [
        figures.cohort_flow(steps, fig_dir),
        figures.rates_plot(rates, fig_dir),
        figures.love_plot(balance, fig_dir),
    ]
    rows = forest_rows(results)
    if rows:
        paths.append(figures.forest_plot(rows, fig_dir))
    print("[run_40] wrote:", *[p.name for p in paths])


if __name__ == "__main__":
    main()
