"""Configuration loading: local data paths and repo layout.

Data paths live in ``config/paths.yaml`` (gitignored). The repo ships
``config/paths.example.yaml`` as a template. Nothing here reads patient data;
it only resolves *where* the (local, credentialed) data lives.
"""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import BaseModel

REPO_ROOT = Path(__file__).resolve().parents[2]


class DataPaths(BaseModel):
    """Resolved locations of the credentialed source databases (local only)."""

    mimic_iv_db: Path
    eicu_dir: Path
    outputs_dir: Path = REPO_ROOT / "outputs"

    def exists(self) -> dict[str, bool]:
        """Report which configured paths are actually present on this machine."""
        return {
            "mimic_iv_db": self.mimic_iv_db.exists(),
            "eicu_dir": self.eicu_dir.exists(),
        }


def load_paths(path: str | Path | None = None) -> DataPaths:
    """Load data paths from a YAML file.

    Resolution order when ``path`` is None:
    1. ``config/paths.yaml`` (real, gitignored) if present
    2. ``config/paths.example.yaml`` (template) as a fallback
    """
    if path is None:
        real = REPO_ROOT / "config" / "paths.yaml"
        example = REPO_ROOT / "config" / "paths.example.yaml"
        path = real if real.exists() else example
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh) or {}
    return DataPaths(**raw)
