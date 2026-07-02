# STROBE / RECORD-RWE checklist

| Item | Addressed in |
|------|--------------|
| Title/abstract, design + summary | Title; Abstract |
| Background/rationale | Introduction |
| Objectives / hypothesis | Introduction (final paragraph) |
| Study design | Methods → Data sources and study design; ADR-0002 |
| Setting / data sources | Methods → Data sources (MIMIC-IV, eICU) |
| Participants + eligibility (RECORD: codes) | Methods → Participants; `config/codelists/`; ADR-0001 |
| Variables: exposure/outcome/covariates | Methods → Exposure, outcome, covariates; ADR-0003 |
| Data sources/measurement | Methods; codelists |
| Bias | Discussion (competing risk, confounding, immortal time) |
| Study size | Results (n=1,771); acknowledged underpowered |
| Quantitative variables (age threshold) | Methods → Exposure |
| Statistical methods | Methods → Statistical methods |
| Confounding control + E-value | Methods; Results (E-value 1.54) |
| Descriptive data | Results; Table 1 (`outputs/table1.csv`) |
| Outcome data | Results (221/1,771 DCI) |
| Main results (OR + CI) | Results (primary 0.86, 0.58–1.28) |
| Other analyses (sex, HRT, source, sensitivity) | Results |
| Key results | Discussion (para 1) |
| Limitations | Discussion (para 2) |
| Interpretation | Discussion |
| Generalisability | Discussion |
| Data/code availability (RECORD) | Data and code availability |

RECORD-specific: complete codelists are version-controlled (`config/codelists/`);
cohort selection is reproducible (`scripts/run_00_cohort.py`); a synthetic fixture
reproduces the pipeline without credentialed data.
