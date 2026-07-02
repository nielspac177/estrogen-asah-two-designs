"""run_00 — assemble the featured cohort and write it to outputs/ (P5-P6)."""

from __future__ import annotations

import json

from estrogen_asah_dci.pipeline import build_featured_cohort, ensure_outputs


def main() -> None:
    out = ensure_outputs()
    cohort, meta = build_featured_cohort()
    banner = "REAL DATA" if meta["mode"] == "real" else "SYNTHETIC DATA (no PhysioNet access)"
    print(f"[run_00] {banner} — cohort n={len(cohort)}")
    cohort.to_parquet(out / "cohort.parquet", index=False)
    (out / "cohort_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[run_00] wrote {out/'cohort.parquet'}")


if __name__ == "__main__":
    main()
