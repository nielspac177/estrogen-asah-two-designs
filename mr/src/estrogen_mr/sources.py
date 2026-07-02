"""Public GWAS summary-statistic sources for the MR (verified available; see docs/adr/0001).

Outcome uses the UK-Biobank-EXCLUDED Bakker files to avoid sample overlap with the
UKB-derived exposures (menopause/menarche/SHBG/testosterone), which would bias a
two-sample MR toward the confounded (observational) estimate.
"""

from __future__ import annotations

OUTCOMES = {
    "IA_cross_ancestry": {
        "study": "Bakker et al. 2020, Nat Genet",
        "phenotype": "Intracranial aneurysm (IA), cross-ancestry",
        "n_cases": 10754, "n_controls": 306882,
        "source": "Figshare 11303372",
        "note": "Use UKB-excluded + rupture-stratified files (aSAH-only vs uIA-only).",
    },
    "aSAH_only_euro_ukb_excluded": {
        "study": "Bakker et al. 2020",
        "phenotype": "Aneurysmal SAH (ruptured), European, UKB-excluded",
        "source": "Figshare 11303372 (European aSAH-only, UKB-excluded)",
        "note": "Primary outcome for the rupture-stratified MR.",
    },
}

EXPOSURES = {
    "age_at_natural_menopause": {
        "study": "Ruth et al. 2021, Nature (ReproGen)",
        "n": 201323,
        "source": "reprogen.org (reprogen_ANM_201K_170621.txt.gz)",
        "note": "CONTINUOUS age at natural menopause: the untested gap. AF removed.",
    },
    "age_at_menarche": {
        "study": "Day et al. 2017, Nat Genet (ReproGen)",
        "n": 370000, "source": "reprogen.org",
    },
    "SHBG": {
        "study": "Ruth et al. 2020, Nat Med",
        "source": "GWAS Catalog GCST90012109 / GCST90012111",
        "note": "Direction-of-effect conflict to resolve (Molenberg 2022 vs JCN 2025).",
    },
    "bioavailable_testosterone": {
        "study": "Ruth et al. 2020, Nat Med",
        "source": "GWAS Catalog GCST90012103",
    },
    "total_testosterone": {
        "study": "Ruth et al. 2020, Nat Med",
        "source": "GWAS Catalog GCST90012113",
    },
    "estradiol": {
        "study": "e.g. Schmitz / UKB estradiol GWAS",
        "note": "Weak instruments (few genome-wide-significant SNPs); report with caution.",
    },
}


def manifest() -> dict:
    return {"outcomes": OUTCOMES, "exposures": EXPOSURES}
