# Adversarial audit synthesis, triangulation paper strategy

6-agent panel (all findings rated **blocker**). The synthesis agent failed on a
schema cap; this is the human-distilled synthesis from the recovered critiques.

## Verdict: achievable, but only reframed as a BOUNDED NEGATIVE

Not "a triangulated null." The honest, defensible thesis:

> Across an ICU cohort and two-sample Mendelian randomization we find no evidence
> that greater cumulative estrogen exposure protects against aneurysm rupture; we
> can **exclude a large protective effect** (per-SD OR < ~0.74 at 80% power), the
> genetic estimate is if anything mildly anti-protective, but **modest protection
> cannot be excluded**. The animal/cohort protective signal does not replicate.

## Mandatory before any inference (blockers)

1. **Positive control (single most decisive analysis).** Run the identical
   pipeline on ANM → a known estrogen-mediated outcome (**breast cancer**, expected
   OR ~1.05/yr later menopause). Recovering it validates build, matching, power;
   failing it means the aSAH null is an artifact. Report BEFORE the aSAH result.
2. **Power / MDE + equivalence, not "null."** Use the true outcome N_eff (~17,000;
   ~4,250 case-equiv). Report per-SD (state SD 3.5–4.5 yr): IVW ≈ OR/SD 1.12
   (0.91–1.38); 80% MDE ≈ OR/SD 0.74. TOST vs SESOI OR 0.90/SD: strong protection
   rejected, harm not, equivalence-to-null fails → bounding language only.
3. **MR QC**: Steiger filtering, MR-PRESSO (global + outliers), leave-one-out /
   single-SNP, random-effects (multiplicative) IVW, per-SNP F, winner's-curse note.
4. **Build/clumping honesty**: 85/85 allele-concordant matches ⇒ same build
   (GRCh37); 65% attrition is outcome SNP coverage, not a build error, state this.
   Distance-clump = documented sensitivity; proper r² clumping if a panel is obtainable.
5. **Sample overlap**: Ruth exposures include UKB → use UKB-excluded Bakker (done);
   note this biases conservative (toward null).

## Reframe the architecture

- **Arm 1 (ICU) demoted** to a design-failure / negative-control demonstration
  (non-identifiable; spec curve shows age-confounding manufactures spurious
  anti-protective significance). It addresses **DCI conditional on admission**.
- **Arm 2 (MR) is the identification arm**; it addresses **aneurysm incidence /
  liability**, NOT DCI (no DCI GWAS exists, state this).
- The arms interrogate **different causal nodes** → drop "triangulation of one
  effect"; use an **estimand table + DAG** (Figure 1) as the load-bearing artifact.
- Strengthen Arm 2 into a **multi-exposure** module (menopause + menarche + SHBG +
  testosterone), one forest; adjudicate the SHBG conflict via MVMR.
- Target: **JAHA / European Stroke Journal / JNNP**, framed as a rigorous bounded null.
