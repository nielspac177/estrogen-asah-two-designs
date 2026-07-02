"""run_10 — Table 1 and crude DCI rates (P7.1-7.2)."""

from __future__ import annotations

import pandas as pd

from estrogen_asah_dci.analysis.models import crude_dci_rates
from estrogen_asah_dci.pipeline import OUTPUTS
from estrogen_asah_dci.report.tables import table_one


def main() -> None:
    cohort = pd.read_parquet(OUTPUTS / "cohort.parquet")
    t1 = table_one(cohort)
    t1.to_csv(OUTPUTS / "table1.csv")
    rates = crude_dci_rates(cohort)
    rates.to_csv(OUTPUTS / "crude_rates.csv", index=False)
    print("[run_10] Table 1:")
    print(t1.to_string())
    print("\n[run_10] Crude DCI rates by stratum:")
    print(rates.to_string(index=False))


if __name__ == "__main__":
    main()
