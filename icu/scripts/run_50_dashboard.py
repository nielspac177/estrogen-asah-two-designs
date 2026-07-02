"""run_50 — build the self-contained results dashboard (P9)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "dashboard"))
from build import build_dashboard  # noqa: E402

from estrogen_asah_dci.config import REPO_ROOT  # noqa: E402
from estrogen_asah_dci.pipeline import OUTPUTS  # noqa: E402


def main() -> None:
    out = build_dashboard(OUTPUTS, REPO_ROOT / "dashboard" / "dist")
    print(f"[run_50] wrote {out} ({out.stat().st_size // 1024} KB, self-contained)")


if __name__ == "__main__":
    main()
