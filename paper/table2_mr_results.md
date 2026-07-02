# Table 2. Mendelian-randomization results

Outcome is aneurysmal SAH (Bakker 2020, European, UK-Biobank-excluded) unless
stated. Distance clumping; random-effects IVW. OR = odds ratio.

| Analysis | Exposure | Instruments | OR (95% CI) | P | Notes |
|---|---|---|---|---|---|
| Primary | Age at natural menopause, per year | 85 | 1.03 (0.97–1.09) | 0.32 | Steiger 80/85; Egger intercept P 0.19; MR-PRESSO P 0.17; leave-one-out 1.01–1.04; mean F 98 |
| Primary, per SD | Age at natural menopause, per SD | 85 | 1.12 (0.91–1.38) | — | 80% power MDE OR 0.73; TOST excludes protection > OR 0.90/SD (P 0.027), not the null/harm (P 0.53) |
| **Positive control** | Age at natural menopause → **breast cancer** | 207 | **1.055 (1.041–1.069)** | 1.8×10⁻²⁵ | recovers the established effect; validates build, matching, power |
| Single-exposure | SHBG, women, per SD | 82 | 0.73 (0.41–1.31) | 0.29 | opposite side of null from Molenberg 2022 (1.18), nearer Tan/Wu 2025 |
| Single-exposure | Total testosterone, per SD | 58 | 0.98 (0.77–1.27) | 0.91 | Egger intercept P 0.04 |
| Multivariable | SHBG, women, per SD (adj. testosterone) | 106 | 1.06 (0.50–2.28) | 0.88 | SHBG point moves 0.73 → 1.06 when testosterone is held fixed |
| Multivariable | Bioavailable testosterone, women, per SD (adj. SHBG) | 106 | 1.03 (0.57–1.86) | 0.93 | |
| Sensitivity | Age at menopause across clumping windows 250 kb–5 Mb | 61–107 | 1.02–1.03 | — | estimate stable to clumping stringency |
| Sensitivity | Age at menopause, r² LD clumping (PLINK, 1000G EUR, r²<0.001/10 Mb) | 81 | 1.03 (0.98–1.09) | 0.22 | gold-standard clumping confirms the distance-clumped primary |

No exposure shows a significant effect on aSAH. The published SHBG direction
conflict is not resolvable at this precision, and the univariable SHBG signal
reflects, in part, its shared testosterone component.
