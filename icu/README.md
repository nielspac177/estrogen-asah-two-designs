# Estrogen state & delayed cerebral ischemia after aneurysmal SAH

A reproducible, multi-database ICU pilot testing whether **endogenous estrogen
state (menopausal status)** is associated with **delayed cerebral ischemia (DCI)
/ cerebral vasospasm** after **aneurysmal subarachnoid hemorrhage (aSAH)**, using
**MIMIC-IV v3.1** and **eICU-CRD v2.0**.

> **Status:** hypothesis-generating pilot. The pooled cohort (~1,000 aSAH
> admissions) is modest and underpowered for finely-stratified DCI analysis.
> The contributions are (1) the first human test of the estrogen→DCI hypothesis,
> (2) a fully reproducible pipeline, and (3) motivation + protocol for a
> definitive longitudinal registry study.

## Why

Animal SAH models show estradiol reduces vasospasm/DCI (ovariectomized mice
rupture 47% vs 7%); human cerebral endothelium expresses ERβ/GPER1. This
protective estrogen→DCI effect **has never been tested in humans**, while the
adjacent causal question (hormone→aSAH Mendelian randomization) is already
published. ICU databases can't emulate a chronic-exposure target trial, but they
*can* test estrogen **state** against **acute** DCI outcomes, their strength.

## Design (summary)

- **Cohort:** adults with aSAH (nontraumatic SAH + aneurysm/securing procedure).
- **Primary exposure:** menopausal state proxy in women (pre <51 vs post ≥51);
  men retained for a sex-difference secondary. Secondary: HRT at admission.
- **Primary outcome:** DCI/vasospasm composite (coded vasospasm ∪ rescue cerebral
  angioplasty/IA vasodilator ∪ delayed cerebral infarction).
- **Analysis:** multilevel logistic regression + propensity/IPW within women,
  source-stratified pooling, E-value sensitivity. See `docs/adr/`.

## Data access (not included here)

This repo contains **code only**. MIMIC-IV and eICU-CRD require PhysioNet
credentialing under the Data Use Agreement:
- MIMIC-IV: https://physionet.org/content/mimiciv/
- eICU-CRD: https://physionet.org/content/eicu-crd/

Copy `config/paths.example.yaml` → `config/paths.yaml` (gitignored) and point it
at your local copies.

## Reproducibility

The full pipeline runs on a **committed synthetic fixture** (`data/synthetic/`,
fake data) with no PhysioNet access, so anyone, and CI, can execute it end to
end:

```bash
uv sync --extra dev
uv run pytest            # unit + synthetic-fixture integration tests
make all                 # runs scripts/run_*.py on the synthetic fixture
```

To run on real data, set `config/paths.yaml` and re-run `make all`.

> **exFAT / network volumes:** if the repo lives on a filesystem that cannot host
> a Python venv (e.g. exFAT), put the environment elsewhere:
> `export UV_LINK_MODE=copy UV_PROJECT_ENVIRONMENT="$HOME/.venvs/estrogen-asah-dci"`
> before `uv sync` / `make`.

## Layout

```
src/estrogen_asah_dci/   extract · harmonize · features · analysis · viz · report
config/                  paths.example.yaml + versioned codelists
scripts/                 numbered end-to-end pipeline (run_00 … run_50)
dashboard/               interactive aggregate-results app (no PHI)
data/synthetic/          committed fake cohort for CI / reproducibility
docs/adr/                architecture decision records
paper/                   manuscript, figures, tables (both arms), one place
tests/                   unit + integration
```

## License

Code: MIT (`LICENSE`). No patient-level data is included or redistributable.
