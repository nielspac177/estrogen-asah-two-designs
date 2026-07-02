"""Versioned clinical codelists and matcher helpers.

Each codelist is a YAML file in ``config/codelists/``. Definitions are data, not
code, so a change to a phenotype is a reviewable diff (see ADR-0001). The matcher
helpers here are reused by both the MIMIC-IV and eICU extractors so the phenotype
logic is identical across sources.
"""

from __future__ import annotations

import re

import yaml
from pydantic import BaseModel, ConfigDict

from .config import REPO_ROOT

CODELIST_DIR = REPO_ROOT / "config" / "codelists"


class Codelist(BaseModel):
    """A single phenotype definition. Extra keys are allowed and preserved."""

    model_config = ConfigDict(extra="allow")

    name: str
    description: str

    def field(self, key: str) -> list[str]:
        """Return a list-valued field (or [] if absent)."""
        return list(getattr(self, key, None) or [])


def load_codelist(name: str) -> Codelist:
    """Load ``config/codelists/{name}.yaml``."""
    path = CODELIST_DIR / f"{name}.yaml"
    with path.open("r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh) or {}
    return Codelist(**raw)


def available_codelists() -> list[str]:
    """Names of all codelists shipped in the repo.

    Skips dotfiles (e.g. macOS AppleDouble ``._*.yaml`` sidecars that appear on
    exFAT/network volumes) so they are never parsed as codelists.
    """
    return sorted(
        p.stem for p in CODELIST_DIR.glob("*.yaml") if not p.name.startswith(".")
    )


def load_all() -> dict[str, Codelist]:
    return {n: load_codelist(n) for n in available_codelists()}


# ---- matcher helpers (shared by extractors) ----

def _norm(code: str) -> str:
    return str(code).replace(".", "").strip().upper()


def code_matches(code: str, prefixes: list[str]) -> bool:
    """True if ``code`` starts with any prefix (ICD codes compared dot-insensitively)."""
    if not prefixes:
        return False
    c = _norm(code)
    return any(c.startswith(_norm(p)) for p in prefixes)


def any_pattern(text: str | None, patterns: list[str]) -> bool:
    """True if any case-insensitive regex in ``patterns`` matches ``text``."""
    if not text or not patterns:
        return False
    return any(re.search(p, text, re.IGNORECASE) for p in patterns)
