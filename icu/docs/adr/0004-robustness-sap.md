# ADR-0004: Robustness / competing-risk analysis plan (audit response)

- Status: accepted
- Date: 2026-07-01

## Context

The primary logistic analysis returned a null (post-vs-pre OR 0.86, 0.58–1.28). An
adversarial methods audit flagged competing mortality (postmenopausal women die
more, suppressing DCI ascertainment), age/menopause collinearity, outcome
misclassification, and researcher degrees of freedom as threats. This ADR
pre-specifies the robustness sweep so any positive finding is interpretable and
not the product of p-hacking.

## Decision

Run a **fixed, fully-reported** set of alternative analyses (`scripts/run_60_robust.py`);
report every result regardless of significance:

1. **Survivor-restricted** (hospital survivors only), competing-mortality fix #1.
2. **Age-restricted 45–55**, isolates menopause from ageing (peri-menopausal band).
3. **Overlap-weighted (ATO)**, exact covariate balance, no extrapolation.
4. **Outcome variants**, vasospasm-only, rescue-procedure-only, delayed-infarction
   (ADR-0003 sensitivity).
5. **Competing-events multinomial**, {neither, DCI, death-without-DCI}; death as a
   separate competing outcome rather than a collider.

## Multiplicity & interpretation rules

- These are **sensitivity / exploratory** analyses, not confirmatory. The
  confirmatory estimand remains the pre-specified primary logistic.
- With ~7 contrasts, apply a Bonferroni threshold of 0.05/7 ≈ 0.007 before calling
  any single result "significant"; otherwise treat as hypothesis-generating.
- **Direction matters more than the p-value.** The mechanistic prior predicts
  postmenopausal (lower estrogen) → MORE DCI (OR > 1). A significant result in the
  *opposite* direction does not confirm the hypothesis and must not be reported as
  if it did.
- No result from this sweep is promoted to a headline without external replication.
