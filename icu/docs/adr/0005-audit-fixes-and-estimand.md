# ADR-0005: Audit fixes, estimand reframing, and the honest verdict

- Status: accepted
- Date: 2026-07-01

## Context

A 7-agent adversarial methodological audit (`docs/audit/adversarial_audit.md`)
reviewed the null primary result and the goal of finding a positive, literature-
consistent finding. Its verdict: the null is credible; the dominant problem is
**non-identifiability** (menopause is a deterministic function of age → no age
overlap across strata), not competing risk; and every "correction" that would
produce the hypothesised OR>1 does so mechanically or by opening a collider.

## Decisions

1. **Bug fix:** the primary outcome model no longer co-adjusts a smooth age term
   with the age-defined menopause step (`models.DEFAULT_COVARIATES` drops `age_c`).
   Co-adjusting age out of an age-defined exposure produced a false-precision CI.
   The outcome model now uses the same confounders as the propensity model.

2. **Estimand reframing:** the estrogen-specific controlled direct effect is NOT
   identifiable in this cohort. We do not claim it. The honest confirmatory
   question becomes the *female-specific deviation in DCI odds across the
   menopausal transition, net of the shared ageing gradient*, the **age×sex
   spline difference-in-differences** (men as ageing reference; `multiverse.age_sex_did`).

3. **Multiverse over fork-selection:** report the menopause OR across all defensible
   forks (`multiverse.spec_curve`) rather than selecting one. Prohibited: elevating
   the eICU-only estimate, sweeping the age cutoff to significance, survivor-
   restriction / mortality adjustment (opens Severity→Death←DCI collider), and
   presenting death-inclusive composites as anything but labelled sensitivities.

4. **Source-stratified reporting:** given the ~10× ascertainment gap and k=2, the
   pooled random-effects OR is demoted; sources are reported stratified.

## Consequences (empirical, real data)

- Specification curve: 71/72 forks converged; **41 significant, all in the
  anti-hypothesis direction (OR<1); zero significant OR>1**. No cherry-pickable
  positive exists.
- Age×sex spline DiD: OR-ratio **1.04 (0.62–1.76)**, no menopause inflection.
- **Verdict: no defensible, literature-consistent positive finding is recoverable
  in this ICU pilot.** The credible causal route is the MR / registry-TTE arms.
