# Estrogen and aneurysmal subarachnoid haemorrhage: companion code

Reproducible code for a two-design study of whether estrogen protects against
aneurysmal subarachnoid haemorrhage (aSAH). This repository is the **companion
code** for the paper: it holds the whole analysis, a step-by-step reproduction
guide, a line-by-line walkthrough of the code, and every figure and table. The
manuscript text itself is not included here.

**Finding.** No measurable protection. A large protective effect is excluded; a
null or a small harmful effect cannot be. The two designs fail in unrelated ways
and reach the same answer.

## Study design

![Study design](figures/study_methods.png)

Two arms, chosen because their weaknesses do not overlap:

- **`icu/`, observational.** Pooled MIMIC-IV + eICU aSAH cohort (n = 1,771),
  menopausal state vs delayed cerebral ischaemia. The design cannot separate
  estrogen from age, which is shown, not hidden, by a specification curve and an
  age-by-sex difference-in-differences model. Adjusted OR 0.86 (0.58–1.28).
- **`mr/`, genetic.** Two-sample Mendelian randomization of genetically predicted
  age at natural menopause on aSAH. IVW OR 1.03 per year (0.97–1.09); a positive
  control (menopause → breast cancer) recovers the known effect and validates the
  pipeline. Multivariable MR and real r² LD clumping confirm the null.

## What is here

| Path | Contents |
|---|---|
| [`REPRODUCE.md`](REPRODUCE.md) | Step-by-step reproduction, including building the MIMIC-IV DuckDB |
| [`INTERNALS.md`](INTERNALS.md) / [`INTERNALS.pdf`](INTERNALS.pdf) | Line-by-line walkthrough of the core analysis code (markdown + rendered PDF) |
| [`FIGURES_AND_TABLES.md`](FIGURES_AND_TABLES.md) | Every figure and table |
| `icu/` | Arm 1 code (60 tests), runs on a committed synthetic fixture |
| `mr/` | Arm 2 code (18 tests), runs on public GWAS |
| `figures/` | Figure generators + rendered PNGs |

## Quick start

```bash
# Arm 1 (ICU) — no credentialed data needed for the tests / synthetic pipeline
cd icu && uv sync --extra dev && uv run pytest && make all

# Arm 2 (MR) — estimator + QC tests need no data
cd mr && uv sync --extra dev && uv run pytest
```

See [`REPRODUCE.md`](REPRODUCE.md) to run on the real MIMIC-IV, eICU, and GWAS data.

## Related repositories

- **[nielspac177/mimic-iv-duckdb](https://github.com/nielspac177/mimic-iv-duckdb)** , 
  build the MIMIC-IV DuckDB database used by the ICU arm.
- **[nielspac177/estrogen-asah-dci](https://github.com/nielspac177/estrogen-asah-dci)** , 
  the ICU arm's development history (step-by-step commits).
- **[nielspac177/estrogen-aneurysm-mr](https://github.com/nielspac177/estrogen-aneurysm-mr)** , 
  the genetic arm's development history.

## Data

No patient-level data or GWAS files are committed. MIMIC-IV and eICU-CRD require
PhysioNet credentialing; the GWAS are public (see `mr/docs/gwas_access.md`).

## Overview and abstract

The full workflow and the graphical abstract:

![Results overview](figures/study_workflow.png)

![Graphical abstract](figures/graphical_abstract.png)

## License

Code: MIT. Uses only de-identified credentialed data and public GWAS summary
statistics under their respective terms.
