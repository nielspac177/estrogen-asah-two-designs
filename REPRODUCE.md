# Reproduce

Two independent arms. The ICU arm needs credentialed MIMIC-IV + eICU; the genetic
arm needs only public GWAS. Everything here runs without either on synthetic data.

## 0. Environment

Install [uv](https://docs.astral.sh/uv/). Each arm is its own project.

```bash
cd icu && uv sync --extra dev   # then, for the other arm:
cd ../mr && uv sync --extra dev
```

On exFAT or network volumes a Python environment cannot live next to the code:

```bash
export UV_LINK_MODE=copy
export UV_PROJECT_ENVIRONMENT="$HOME/.venvs/estrogen-icu"   # per arm
```

Every step below can first be checked with `uv run pytest` (no data needed).

## 1. Arm 1, observational (ICU)

### 1a. Build the MIMIC-IV DuckDB

The ICU extractor reads MIMIC-IV from a DuckDB database. Build it once with the
companion tool:

- **[nielspac177/mimic-iv-duckdb](https://github.com/nielspac177/mimic-iv-duckdb)** , 
  follow its README to load the PhysioNet MIMIC-IV v3.1 CSVs into `mimic4.db`
  (schemas `mimiciv_hosp`, `mimiciv_icu`).

You also need the eICU-CRD v2.0 CSVs (the gzipped files as distributed by
PhysioNet); the extractor reads them by name, no database required.

### 1b. Point the config at your local data

```bash
cd icu
cp config/paths.example.yaml config/paths.yaml   # gitignored
# edit config/paths.yaml:
#   mimic_iv_db: /path/to/mimic4.db
#   eicu_dir:    /path/to/eicu-collaborative-research-database-2.0
```

### 1c. Run the pipeline

```bash
make all          # runs scripts/run_00 .. run_70 in order; real data if paths.yaml is set,
                  # otherwise a labelled synthetic cohort
```

Outputs (git-ignored) land in `icu/outputs/`: `cohort.parquet`, `table1.csv`,
`results.json` (primary + secondary models, E-value), `results_robust.json`
(competing-risk / sensitivity sweep), `results_multiverse.json` (specification
curve + age×sex difference-in-differences), `balance.csv`, figures, and the
`dashboard/dist/index.html`.

Key numbers on real data: 1,771 aSAH admissions; primary adjusted OR 0.86
(0.58–1.28); specification curve 41/71 significant forks, all anti-protective;
age×sex difference-in-differences 1.04 (0.62–1.76).

## 2. Arm 2, genetic (Mendelian randomization)

### 2a. Download the GWAS

Public summary statistics only. `mr/docs/gwas_access.md` gives exact URLs and the
build/format gotchas. In short, into `mr/data/gwas/`:

- Outcome: Bakker 2020 intracranial-aneurysm GWAS (Figshare 11303372), the
  `SAH-only_European_excludingUKBB` file.
- Exposures: age at natural menopause (Ruth 2021, ReproGen); SHBG and bioavailable
  testosterone in women (Ruth 2020, GWAS Catalog GCST90012107 / GCST90012102,
  Build37 files).
- Positive control: breast cancer (Michailidou 2017, GCST004988, the **Build37**
  file, not the GRCh38 harmonised one).

### 2b. Run the analyses

```bash
cd mr
uv run python scripts/run_mr_full.py         # ANM -> aSAH, full QC + power + equivalence
uv run python scripts/run_mr_poscontrol.py   # positive control: ANM -> breast cancer
uv run python scripts/run_mr_multiexposure.py # SHBG + testosterone -> aSAH
uv run python scripts/run_mvmr.py            # multivariable MR (SHBG + bioT)
uv run python scripts/run_mr_clumping_sensitivity.py
```

For gold-standard r² LD clumping (`scripts/run_mr_ld_clump.py`) you also need
PLINK 2 and a 1000 Genomes European reference in `mr/data/ldref/` (EUR.bed/bim/fam
+ the `plink2` binary); the script maps the chr:pos instruments to rsIDs via the
reference bim and clumps at r² < 0.001 / 10 Mb.

Key numbers: primary IVW OR 1.03/yr (0.97–1.09); positive control weighted-median
1.055/yr (1.041–1.069); MVMR SHBG 1.06, bioT 1.03; r²-clumped 1.03 (0.98–1.09).

## 3. Regenerate the figures

```bash
python figures/study_methods.py
python figures/study_workflow.py
python figures/mr_forest.py
python figures/graphical_abstract.py
```

The design and abstract figures are self-contained. The forest embeds the final MR
estimates. All rendered PNGs are committed in `figures/` and shown in
[`FIGURES_AND_TABLES.md`](FIGURES_AND_TABLES.md).
