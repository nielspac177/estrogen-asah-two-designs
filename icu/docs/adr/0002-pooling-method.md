# ADR-0002: Pooling MIMIC-IV and eICU

- Status: accepted
- Date: 2026-07-01

## Context

The two sources differ in era (MIMIC-IV mixes ICD-9/10; eICU is APACHE-coded),
coding of DCI, available severity (eICU has APACHE IV, MIMIC does not here), and
site structure (MIMIC ≈ one hospital; eICU = many). Naive row-pooling risks
Simpson-style confounding by source.

## Decision

Harmonize to one schema, then analyze with **`source` (and eICU `hospital_id`)
as a clustering / stratification variable**, not as an ignorable covariate:

- Primary model: multilevel logistic regression with a random intercept for
  `hospital_id` (MIMIC contributes a single cluster).
- Report **source-stratified estimates** alongside the pooled estimate and a
  fixed/random-effects meta-analytic combination; disagreement between sources is
  a finding, not noise.
- Patient IDs are namespaced (`mimic_iv:*`, `eicu:*`) so pooling cannot create
  collisions; eICU is deduplicated by `uniquepid` upstream.

## Consequences

- One combined frame feeds every downstream step, but source is always available
  for stratification and as a random effect.
- Severity is only partially comparable (APACHE in eICU, absent in MIMIC); severity
  adjustment is therefore source-aware (see analysis; limitation to report).
