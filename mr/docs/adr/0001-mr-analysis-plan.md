# ADR-0001: MR analysis plan (pre-specified)

- Status: accepted
- Date: 2026-07-01

## Context

Two published two-sample MRs already tested sex hormones vs aneurysm/aSAH
(Molenberg 2022, Stroke; Tan/Wu 2025, J Clin Neurosci) and **disagree on the SHBG
direction of effect**. Continuous age-at-natural-menopause (Ruth 2021) has never
been MR-tested against rupture-stratified intracranial-aneurysm outcomes. This is
the genuine, narrow gap this arm addresses. It is causal and, unlike the ICU
pilot, immune to age/ascertainment confounding.

## Question & estimand

Causal effect of genetically-proxied (a) age at natural menopause, (b) SHBG,
(c) bioavailable/total testosterone, (d) age at menarche, (e) estradiol on
intracranial aneurysm (IA) and, stratified, on aneurysmal SAH (rupture) vs
unruptured IA. Estimand: per-SD (or per-year for menopause) causal log-OR.

## Design (locked before fitting)

- **Two-sample MR.** Outcome = Bakker 2020 **UKB-excluded** files (avoids overlap
  with UKB-derived exposures); report cross-ancestry and rupture-stratified.
- **Instruments:** genome-wide significant (p<5e-8), LD-clumped r²<0.001/10 Mb;
  report F-statistic (exclude weak instruments F<10); Steiger filtering to drop
  SNPs explaining more outcome than exposure variance.
- **Primary estimator:** IVW (random effects). **Sensitivity:** MR-Egger
  (intercept = directional pleiotropy), weighted median, MR-PRESSO (outliers),
  leave-one-out. Concordance across methods is required before any causal claim.
- **Heterogeneity:** Cochran's Q, I².
- **Multivariable MR** for SHBG vs bioavailable testosterone to disentangle the
  androgen/SHBG axis (the source of the published conflict).
- **Multiplicity:** primary = age-at-menopause → aSAH; others are secondary with
  FDR across the exposure×outcome grid.

## Reproducibility

- Estimators implemented and unit-tested against synthetic instruments with a
  known causal effect (`simulate.py`, `tests/`).
- GWAS are NOT committed (large; some restricted redistribution). `sources.py`
  documents exact provenance; a fetch step downloads them locally into `data/gwas/`.
- Everything downstream of harmonized summary stats runs deterministically.

## Interpretation guard

A positive finding is credible only if IVW, weighted median, and MR-Egger agree in
direction, the Egger intercept is null, and it survives leave-one-out, and even
then it is a genetic-liability estimate, to be triangulated with the (future)
registry target-trial emulation, not equated with an HRT treatment effect.
