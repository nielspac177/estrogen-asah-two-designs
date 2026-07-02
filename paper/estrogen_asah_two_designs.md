# Estrogen exposure and aneurysmal subarachnoid haemorrhage: no protective effect across an intensive-care cohort and Mendelian randomization

**Authors.** N. Pacheco (collaborators to be added).

**Reporting.** STROBE (observational arm) and STROBE-MR (genetic arm). Code and a
synthetic-data reproduction are public (see Data and code availability).

---

## Abstract

**Background.** Experimental subarachnoid-haemorrhage models show that oestradiol
reduces cerebral vasospasm and delayed cerebral ischaemia, and cohort studies link
earlier menopause to higher subarachnoid-haemorrhage incidence. Whether greater
estrogen exposure protects against aneurysmal subarachnoid haemorrhage (aSAH) in
humans has not been tested with a design that separates estrogen from age.

**Methods.** We estimated the association two ways whose biases do not overlap.
The first was an observational analysis of a pooled intensive-care cohort
(MIMIC-IV and eICU, 1,771 aSAH admissions), using menopausal state as an
age-based proxy for estrogen and delayed cerebral ischaemia or death as the
outcome. The second was two-sample Mendelian randomization of genetically
predicted age at natural menopause (Ruth 2021) on aSAH liability (Bakker 2020,
European, UK-Biobank-excluded), with a pre-specified positive control (breast
cancer, an established consequence of later menopause) run through the identical
pipeline. We report power, minimum detectable effects, and equivalence tests
rather than describing results as simple nulls.

**Results.** The observational estimate was null (adjusted odds ratio for
postmenopausal versus premenopausal women 0.86, 95% CI 0.58 to 1.28), but the
design cannot identify an estrogen effect: menopausal state is a deterministic
function of age here, no biological measure of menopause was recorded (estradiol
in 1 of 354 women), and a specification curve of 72 model choices produced 41
significant results, all in the direction opposite to the protective hypothesis,
consistent with age acting on vasospasm rather than estrogen. In Mendelian
randomization, the positive control recovered the known effect (weighted-median
odds ratio for breast cancer 1.055 per year of later menopause, 95% CI 1.041 to
1.069; inverse-variance weighted P = 1.8 x 10^-25; 207 instruments), which
validates the instruments, genome build, and estimators. Against aSAH the same
instruments gave no protective effect (random-effects inverse-variance weighted
odds ratio 1.03 per year, 95% CI 0.97 to 1.09; approximately 1.12 per standard
deviation). Directionality (Steiger 80 of 85), pleiotropy (Egger intercept
P = 0.19; MR-PRESSO global P = 0.17), and leave-one-out checks were unremarkable.
Equivalence testing excluded protection stronger than an odds ratio of 0.90 per
standard deviation (P = 0.027) but could not exclude the null or mild harm
(P = 0.53); power to detect an odds ratio of 0.90 per standard deviation was 15%.

**Conclusions.** Two designs with non-overlapping weaknesses agree that estrogen
does not measurably protect against aneurysmal subarachnoid haemorrhage in humans.
We can exclude a large protective effect; a null or a small harmful effect remains
possible. The protective signal seen in animal models and reproductive-timing
cohorts does not replicate at genetic scale.

---

## Introduction

Delayed cerebral ischaemia after aneurysmal subarachnoid haemorrhage is a major
driver of poor outcome, and the search for modifiable protection has repeatedly
returned to estrogen. Ovariectomised animals rupture and infarct more often,
17-beta-oestradiol acts on endothelial nitric-oxide and endothelin-1 signalling in
haemorrhage models, and human cerebral vessels carry estrogen receptors.[1,2] In
population data, women with a shorter reproductive span or earlier menopause have
a higher risk of subarachnoid haemorrhage, and incidence rises around the age of
menopause.[3] Together these observations suggest that higher estrogen exposure
might protect against aneurysm formation, rupture, or the ischaemia that follows.

Testing that idea in humans is hard for one reason: the exposure people most want
to study, lifetime estrogen, is bound to age, and age is itself a strong and
opposite influence on cerebral vasospasm, which is more common in younger patients.
An analysis that compares older and younger women therefore cannot tell an estrogen
effect apart from an age effect. The adjacent genetic question, whether sex-hormone
biology influences aneurysm risk, has been examined before, but the two published
Mendelian-randomization studies disagree on the direction of the sex-hormone-binding
globulin effect,[4,5] and continuous age at natural menopause has not been tested
against rupture-specific outcomes.

We approached the question with two designs chosen because their weaknesses do not
overlap (Figure 1). The first is an observational intensive-care analysis, which we
present mainly to show why the intuitive approach fails. The second is Mendelian
randomization, which uses genetic variants fixed at conception as instruments for
menopausal timing and is therefore free of the age confounding that defeats the
observational arm. Reading the two together (Figure 2), an approach known as
triangulation, is more informative than either alone because the designs fail in
unrelated ways.[16]

## Methods

### Estimands

The two arms do not estimate the same quantity, and we state this plainly rather
than treat them as interchangeable (Table 1). The observational arm estimates the
association of menopausal state with delayed cerebral ischaemia conditional on
admission with aSAH. The genetic arm estimates the causal effect of genetically
predicted age at natural menopause on aSAH liability. No genome-wide association
study of delayed cerebral ischaemia exists, so Mendelian randomization cannot
address the ischaemia endpoint directly; it speaks to aneurysm rupture.

### Observational arm

We pooled adult aSAH admissions from MIMIC-IV v3.1 and eICU-CRD v2.0,[6,7] defined by
a nontraumatic subarachnoid-haemorrhage code together with an aneurysm code or an
aneurysm-securing procedure, with traumatic and arteriovenous-malformation cases
excluded. Menopausal state was defined by age (premenopausal below 51,
postmenopausal at or above 51), with men retained as a reference group. The
primary outcome was a composite of coded vasospasm, rescue endovascular therapy,
or delayed cerebral infarction; secondary outcomes were in-hospital death and poor
discharge disposition. We fitted multilevel logistic regression with
inverse-probability weighting, then probed researcher discretion with a
specification curve across menopause cut-offs, covariate sets, outcome definitions,
and data sources, and with an age-by-sex difference-in-differences model that uses
men to absorb the shared ageing gradient.

### Genetic arm

We used two-sample Mendelian randomization. Instruments for age at natural
menopause came from Ruth et al. 2021 (ReproGen),[8] selected at genome-wide
significance and distance-clumped; we note that distance clumping approximates
linkage-disequilibrium clumping and report it as such. The outcome was aSAH
liability from Bakker et al. 2020,[9] European ancestry, with UK Biobank samples
excluded to avoid overlap with the
UK-Biobank-derived exposure, a choice that biases toward the null and is therefore
conservative. Effects were harmonised to a common allele and matched on genome
position; both datasets are GRCh37. We estimated the effect with random-effects
inverse-variance weighting and tested robustness with MR-Egger (intercept for
directional pleiotropy), the weighted median, MR-PRESSO, Steiger filtering, and
leave-one-out.[10-13] Reporting follows STROBE-MR.[14] We report the effect per year and per standard deviation of
menopausal age (taken as four years), the minimum detectable effect at 80% power,
and a two-one-sided equivalence test against a smallest effect of interest of an
odds ratio of 0.90 per standard deviation, anchored to the reproductive-timing
cohort literature. We repeated the analysis for sex-hormone-binding globulin and
total testosterone in women (Ruth 2020), and ran a multivariable analysis of
sex-hormone-binding globulin with bioavailable testosterone to separate the two
overlapping instruments. Robustness to instrument selection was checked by varying
the distance-clumping window from 250 kb to 5 Mb and by re-running with
reference-panel r-squared clumping (PLINK, 1000 Genomes European panel, r-squared
below 0.001 over 10 Mb).

### Positive control

Because a null is only interpretable if the pipeline can detect a real effect, we
ran the identical instruments and code against breast cancer (Michailidou et al.
2017, 76,192 cases and 63,082 controls),[15] where later menopause is an
established risk factor of about a 5% increase in odds per year.

## Results

### The observational arm cannot identify an estrogen effect

Of 1,771 aSAH admissions (1,105 women and 666 men), 221 met the delayed cerebral
ischaemia composite. The adjusted odds ratio for postmenopausal versus
premenopausal women was 0.86 (95% CI 0.58 to 1.28). This estimate is not
interpretable as an estrogen effect. Menopausal state was defined by age, so the
two cannot be separated, and no biological marker of menopause was available:
estradiol was measured in 1 of 354 women, follicle-stimulating hormone in 6, a
menopause diagnosis appeared in 4, and eICU held no hormone assays at all. The
specification curve made the problem visible. Of 72 model choices, 41 reached
significance, and every one of them pointed opposite to the protective hypothesis,
which is the signature of age acting on vasospasm rather than of estrogen. The
age-by-sex difference-in-differences model, the one specification that removes the
shared ageing gradient, was null (odds-ratio ratio 1.04, 95% CI 0.62 to 1.76). The
arm was also underpowered for the effect of interest (power near 11% at its own
estimate). We therefore treat it as a demonstration of why the observational
approach fails, not as evidence about estrogen.

### The genetic arm recovers a known effect, then finds none for aSAH

The positive control worked. Age at natural menopause raised breast-cancer odds by
about the established amount (weighted-median odds ratio 1.055 per year, 95% CI
1.041 to 1.069; inverse-variance weighted P = 1.8 x 10^-25; 207 instruments). This
confirms that the instruments are live, the genome build is correct, and the
estimators recover truth. The same instruments matched fewer variants in the
smaller aSAH dataset (85 of 242), which reflects that outcome's sparser coverage
rather than a build error, since the same-build breast-cancer file matched 207.

Against aSAH there was no protective effect. The random-effects inverse-variance
weighted odds ratio was 1.03 per year of later menopause (95% CI 0.97 to 1.09),
about 1.12 per standard deviation. MR-Egger and the weighted median agreed in
direction, the Egger intercept gave no sign of directional pleiotropy (P = 0.19),
MR-PRESSO found no global outlier signal (P = 0.17), Steiger placed 80 of 85
variants in the correct causal direction, and leave-one-out estimates stayed
between 1.01 and 1.04, so no single variant drove the result. Mean instrument
strength was high (F near 98).

Power, not effect, is the honest frame. At 80% power the smallest protective effect
we could have detected was an odds ratio of about 0.73 per standard deviation;
power to detect a more plausible 0.90 per standard deviation was 15%. The
equivalence test excluded protection stronger than an odds ratio of 0.90 per
standard deviation (P = 0.027) but could not exclude the null or a mild harmful
effect (P = 0.53). The point estimate leans, without significance, toward higher
rather than lower risk.

Widening the genetic arm to other sex-hormone exposures did not change the picture
(Figure 3). Neither sex-hormone-binding globulin in women (odds ratio 0.73, 95% CI
0.41 to 1.31) nor total testosterone (0.98, 95% CI 0.77 to 1.27) reached
significance against aSAH. The single-exposure SHBG estimate sat on the opposite
side of the null from Molenberg 2022 (which reported 1.18) and closer to Tan/Wu
2025, but its interval spans both. A multivariable analysis of SHBG together with
bioavailable testosterone, which separates the two overlapping instruments,
returned direct effects near the null for both (SHBG 1.06, 95% CI 0.50 to 2.28;
bioavailable testosterone 1.03, 95% CI 0.57 to 1.86), and moved the SHBG point
estimate from 0.73 to 1.06. The published disagreement is therefore not resolvable
in these data; the SHBG signal appears in part to reflect its shared testosterone
component. The primary menopause estimate was stable across clumping windows from
250 kb to 5 Mb (odds ratio 1.02 to 1.03), and reference-panel r-squared clumping
against a 1000 Genomes European panel gave the same answer (odds ratio 1.03, 95%
CI 0.98 to 1.09; 81 instruments), so the result does not depend on the clumping
method.

### Reading the two arms together

Both designs return no protective effect, and their weaknesses do not overlap. If
age confounding were hiding a real protective effect in the observational arm, the
genetic arm, which age cannot touch, would have revealed it; it did not. If the
genetic instruments were dead or mis-built, the breast-cancer control would have
been null; it was strongly positive. The convergence is therefore not two
underpowered nulls agreeing by accident, but agreement between a design that fails
by confounding and a design that succeeds on a validated positive control.

## Discussion

Across an intensive-care cohort and two-sample Mendelian randomization we find no
evidence that greater estrogen exposure protects against aneurysmal subarachnoid
haemorrhage. We can exclude a large protective effect. We cannot exclude a null or
a small harmful effect, and modest protection near an odds ratio of 0.90 per
standard deviation remains within the data's reach only weakly, given 15% power.
The protective signal from animal models and reproductive-timing cohorts does not
appear at genetic scale for aneurysm rupture.

The result also shows how the intuitive analysis misleads. The observational arm
produced many significant associations, all pointing away from protection, because
age both defines the exposure and independently raises vasospasm risk. Presenting
that arm as an estrogen result would have been wrong in a way that a specification
curve and a difference-in-differences model make obvious.

Several limitations bound the reading. The genetic outcome is aneurysm rupture, not
delayed cerebral ischaemia, because no genome-wide association study of ischaemia
exists; the acute neuroprotection hypothesis from animal work is therefore still
untested at the level of the ischaemic complication. The aSAH outcome is
sex-combined, which dilutes a female-specific exposure. Instrument selection used
distance clumping rather than reference-panel clumping, and 85 of 242 variants
matched the aSAH dataset. Age at natural menopause captures one axis of estrogen
biology; sex-hormone-binding globulin and testosterone, where prior studies
conflict, are natural next exposures in the same framework. Finally, the
observational cohort is an intensive-care population and conditions on survival to
admission.

The multi-exposure and multivariable analyses reported here already widen the
genetic arm to sex-hormone-binding globulin and testosterone and show that the
prior disagreement is not resolvable at current precision. The remaining step, for
the ischaemia question specifically, is a longitudinal target-trial emulation of
hormone therapy in data with outpatient dispensing and incident coding, which the
present sources cannot provide.

## Data and code availability

All analysis code and a synthetic-data reproduction that runs the full pipeline
without credentialed data are public. The observational arm requires PhysioNet
credentialing for MIMIC-IV and eICU-CRD. The genetic arm uses public summary
statistics (Bakker 2020 via Figshare; Ruth 2021 via ReproGen; Michailidou 2017 via
the GWAS Catalog); see the repository's data-access notes.

## Figures and tables

- **Graphical abstract.** Question, the two designs with their estimates, and the
  convergent conclusion (`figures/graphical_abstract.png`).
- **Figure 1.** Study design: the two arms, their inputs, and their analyses
  (`figures/study_methods.png`).
- **Figure 2.** Results overview: convergent evidence across the two designs
  (`figures/study_workflow.png`).
- **Figure 3.** Mendelian-randomization forest: sex-hormone exposures on aSAH, with
  the positive control (`figures/mr_forest.png`).
- **Table 1.** Estimand comparison (exposure, outcome node, population, dominant
  bias, and what each arm can and cannot identify).
- **Table 2.** Mendelian-randomization results (`table2_mr_results.md`).

## Statements

**Ethics.** MIMIC-IV and eICU-CRD are de-identified and were used under PhysioNet
credentialed access; secondary analysis is exempt from further review. The genetic
arm used published, publicly available summary statistics only.

**Funding.** [to be completed].

**Competing interests.** The authors declare no competing interests.

**Data and code availability.** Analysis code and a synthetic-data reproduction are
public at github.com/nielspac177/estrogen-asah-dci (observational arm) and
github.com/nielspac177/estrogen-aneurysm-mr (genetic arm). Individual patient data
require PhysioNet credentialing; GWAS summary statistics are public (see above).

**Author contributions.** [to be completed].

## References

*Citations are provisional placeholders; verify against source records before
submission.*

1. Tada Y, Wada K, Shimada K, et al. Estrogen protects against intracranial
   aneurysm rupture in ovariectomised mice. Hypertension. 2014.
2. Krause DN, Duckles SP, Pelligrino DA. Influence of sex steroid hormones on the
   cerebrovascular system. J Appl Physiol. 2006.
3. Reproductive factors and risk of aneurysmal subarachnoid haemorrhage in the
   Nurses' Health Study cohort. Neurology. 2022. (PMC9162048).
4. Molenberg R, Aalbers MW, Appelman APA, et al. Sex hormones and risk of
   aneurysmal subarachnoid haemorrhage: a Mendelian randomization study. Stroke.
   2022;53:2870-2875.
5. Tan X, Wu Q, et al. Sex hormones and cerebral aneurysm / subarachnoid
   haemorrhage: a Mendelian randomization study. J Clin Neurosci. 2025;136:111244.
6. Johnson AEW, Bulgarelli L, Shen L, et al. MIMIC-IV, a freely accessible
   electronic health record dataset. Sci Data. 2023;10:1.
7. Pollard TJ, Johnson AEW, Raffa JD, et al. The eICU Collaborative Research
   Database. Sci Data. 2018;5:180178.
8. Ruth KS, Day FR, Hussain J, et al. Genetic insights into biological mechanisms
   governing human ovarian ageing. Nature. 2021;596:393-397.
9. Bakker MK, van der Spek RAA, van Rheenen W, et al. Genome-wide association study
   of intracranial aneurysms identifies 17 risk loci. Nat Genet. 2020;52:1303-1313.
10. Bowden J, Davey Smith G, Burgess S. Mendelian randomization with invalid
    instruments: effect estimation and bias detection through Egger regression.
    Int J Epidemiol. 2015;44:512-525.
11. Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent estimation in
    Mendelian randomization with the weighted median estimator. Genet Epidemiol.
    2016;40:304-314.
12. Verbanck M, Chen CY, Neale B, Do R. Detection of widespread horizontal
    pleiotropy in causal relationships (MR-PRESSO). Nat Genet. 2018;50:693-698.
13. Hemani G, Tilling K, Davey Smith G. Orienting the causal relationship between
    imprecisely measured traits using GWAS summary data. PLoS Genet. 2017;13:e1007081.
14. Skrivankova VW, Richmond RC, Woolf BAR, et al. Strengthening the reporting of
    observational studies in epidemiology using Mendelian randomization (STROBE-MR).
    JAMA. 2021;326:1614-1621.
15. Michailidou K, Lindström S, Dennis J, et al. Association analysis identifies 65
    new breast cancer risk loci. Nature. 2017;551:92-94.
16. Lawlor DA, Tilling K, Davey Smith G. Triangulation in aetiological
    epidemiology. Int J Epidemiol. 2016;45:1866-1886.
