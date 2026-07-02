"""run_30 — propensity-score covariate balance (P7.4)."""

from __future__ import annotations

import pandas as pd

from estrogen_asah_dci.analysis.weighting import balance_table
from estrogen_asah_dci.pipeline import OUTPUTS


def main() -> None:
    cohort = pd.read_parquet(OUTPUTS / "cohort.parquet")
    bt = balance_table(cohort)
    bt.to_csv(OUTPUTS / "balance.csv", index=False)
    print("[run_30] SMD balance (pre vs post, women):")
    print(bt.to_string(index=False))


if __name__ == "__main__":
    main()
