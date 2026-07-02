# ADR-0001: Cohort, exposure, and outcome definitions

- Status: accepted
- Date: 2026-07-01

## Context

We are running a hypothesis-generating pilot on whether endogenous estrogen state
is associated with delayed cerebral ischemia (DCI) after aneurysmal subarachnoid
hemorrhage (aSAH), pooling MIMIC-IV v3.1 and eICU-CRD v2.0. Definitions must be
codelist-driven, versioned, and identical in spirit across two very different
schemas. Feasibility counts were established by direct queries (see
`estrogen_aneurysm_tte/ROADMAP.md`).

## Decision

**Cohort, aneurysmal SAH, adults (≥18):**
- MIMIC-IV: nontraumatic SAH (ICD-9 `430`, ICD-10 `I60*`) AND evidence of an
  aneurysm, a cerebral-aneurysm diagnosis (`codelists/aneurysm`) OR an
  aneurysm-securing procedure (clip/coil, `codelists/aneurysm_procedures`).
- eICU: APACHE/admission dx "Subarachnoid hemorrhage/intracranial aneurysm"
  (± "surgery for") OR a hemorrhagic-stroke SAH `diagnosisstring`
  (`codelists/eicu_asah_strings`), especially "from ruptured berry aneurysm".
- **Exclusions:** traumatic SAH, AVM-related SAH. One row per patient-admission;
  eICU deduplicated by `uniquepid` to the first qualifying unit stay.

**Primary exposure, menopausal state proxy (women):** premenopausal (age < 51)
vs postmenopausal (age ≥ 51). Age 51 ≈ median US natural menopause. Men retained
as a reference arm for a sex-difference secondary. eICU age ">89" → 90.

**Secondary exposure (exploratory, low power):** HRT/estrogen at admission, 
eICU `admissionDrug` home meds; MIMIC-IV `prescriptions` started within 48 h
(`codelists/hrt`).

**Primary outcome, DCI/vasospasm composite (any of):**
1. coded cerebral vasospasm (ICD-10 `I67848`/`I67841`; eICU diagnosis "…with
   vasospasm"), `codelists/vasospasm`;
2. rescue cerebral angioplasty / intra-arterial vasodilator
   (`codelists/dci_procedures`);
3. delayed cerebral infarction (where timing supports a >72 h window).

**Secondary outcomes:** in-hospital mortality; poor discharge disposition
(death/hospice); ICU length of stay.

## Consequences

- Definitions live in `config/codelists/*.yaml`, each version-controlled and unit
  tested; changing a definition is a reviewable diff.
- ICD-9 has no clean cerebral-vasospasm code, so DCI in the MIMIC ICD-9 era leans
  on procedures + delayed infarction, a known sensitivity limitation, to be
  reported.
- The composite trades specificity for power; sensitivity analyses will vary its
  components (ADR-0003).
