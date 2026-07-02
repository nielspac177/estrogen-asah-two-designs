# Table 1. Estimand comparison across the two designs

| | Arm 1, Observational (ICU) | Arm 2, Genetic (Mendelian randomization) |
|---|---|---|
| **Exposure** | Menopausal state (age-based proxy: <51 vs ≥51) | Genetically predicted age at natural menopause (Ruth 2021 instruments) |
| **Outcome node** | Delayed cerebral ischaemia / vasospasm, conditional on aSAH admission | aSAH liability (aneurysm rupture; Bakker 2020) |
| **Population** | Adults admitted to ICU with aSAH (MIMIC-IV + eICU) | European-ancestry participants of the source GWAS |
| **Dominant bias** | Confounding by age (exposure is defined by age); ICU selection; survival to admission | Horizontal pleiotropy of menopause variants |
| **Direction of that bias** | Age raises vasospasm risk in the young → pushes the estimate *against* the protective hypothesis | Pleiotropy direction unknown; tested by Egger/MR-PRESSO/Steiger |
| **What it can identify** | Whether menopausal-state (mostly age) tracks DCI in admitted patients | Causal effect of menopausal timing on rupture, free of age confounding |
| **What it cannot identify** | An estrogen effect separable from age; incidence (only conditional-on-admission) | The DCI/ischaemia node (no DCI GWAS exists); a female-specific effect (outcome is sex-combined) |
| **Result** | Adjusted OR 0.86 (0.58–1.28), non-identifiable; spec curve significance all anti-protective | IVW OR 1.03/yr (0.97–1.09); excludes protection > OR 0.90/SD, not the null |
| **Validation** | Specification curve + age×sex difference-in-differences (null, 1.04) | Positive control ANM→breast cancer recovered (OR 1.055/yr, P=1.8×10⁻²⁵) |

Because the two arms estimate different causal nodes, their agreement is read as
convergent evidence against a protective estrogen effect, not as two estimates of
one quantity.
