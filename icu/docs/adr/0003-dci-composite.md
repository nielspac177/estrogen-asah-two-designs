# ADR-0003: DCI / vasospasm composite outcome

- Status: accepted
- Date: 2026-07-01

## Context

Delayed cerebral ischemia has no single reliable code, especially in the ICD-9
era. Candidate signals: coded cerebral vasospasm, rescue endovascular therapy
(cerebral angioplasty / intra-arterial vasodilator), and delayed cerebral
infarction. Each alone is insensitive; nimodipine is standard of care for all
aSAH and non-discriminating.

## Decision

Primary outcome = **composite (any of):** coded vasospasm ∪ rescue cerebral
angioplasty/IA vasodilator ∪ delayed cerebral infarction. Rationale: maximize
sensitivity for clinically significant vasospasm/DCI in a modest cohort.

Pre-specified **sensitivity analyses vary the composite:**
1. vasospasm-coded only (most specific);
2. rescue-procedure only (captures treated, clinically significant vasospasm);
3. exclude delayed infarction (which can have non-DCI causes).

Nimodipine is explicitly excluded as a marker (see `codelists/dci_procedures`).

## Consequences

- One boolean `dci_composite` feature drives the primary analysis; components are
  retained so sensitivity variants are a filter, not a re-extraction.
- Composite specificity is imperfect (a known limitation to report); the
  sensitivity variants bound its effect on the estimate.
