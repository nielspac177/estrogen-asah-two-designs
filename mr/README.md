# Estrogen / sex-hormones and intracranial aneurysm, Mendelian randomization

Two-sample MR of sex-hormone and reproductive-timing exposures (age at natural
menopause, SHBG, bioavailable/total testosterone, age at menarche, estradiol) on
**intracranial aneurysm (IA)** and **aneurysmal SAH**, using public GWAS summary
statistics. Companion causal arm to the `estrogen-asah-dci` ICU pilot, immune to
the age/ascertainment confounding that limits observational ICU data.

**Gap addressed** (see `docs/adr/0001`): resolve the published SHBG direction-of-
effect conflict (Molenberg 2022 vs Tan/Wu 2025) and test *continuous* age at
natural menopause (Ruth 2021) against rupture-stratified Bakker 2020 outcomes
(UKB-excluded, to avoid sample overlap).

## Status

Primary analysis complete. The GWAS have been downloaded, instruments selected and
harmonized, and the full analysis plan run, including sensitivity, power, and a
passing positive control. What remains is optional extension work, not a blocker.

### Results (age at natural menopause -> aneurysmal SAH)

No protective effect. Random-effects IVW odds ratio 1.03 per year of later
menopause (95% CI 0.97 to 1.09), about 1.12 per standard deviation. 85 instruments,
mean F 98, Steiger 80/85 in the correct direction, MR-Egger intercept P 0.19,
MR-PRESSO global P 0.17, leave-one-out stable (1.01 to 1.04). Equivalence testing
excludes protection stronger than an odds ratio of 0.90 per standard deviation
(P 0.027) but cannot exclude the null or mild harm (P 0.53); power to detect an
odds ratio of 0.90 per standard deviation was 15%. Report the finding as a bound,
not as "no effect".

### Positive control (age at natural menopause -> breast cancer): PASSES

The identical pipeline recovers the established effect of later menopause on breast
cancer: weighted-median odds ratio 1.055 per year (95% CI 1.041 to 1.069), IVW
P 1.8e-25, 207 instruments (Michailidou 2017, GRCh37 build). This validates the
instruments, genome build, matching, and power, so the aSAH null is credible
evidence rather than a broken pipeline.

### Multi-exposure MR (done)

Same pipeline over three sex-hormone exposures against aSAH: none is significant.
Age at natural menopause 1.03 (0.97 to 1.09). SHBG in women 0.73 (0.41 to 1.31),
which leans opposite to Molenberg 2022 (OR 1.18) and toward Tan/Wu 2025, but the
interval is too wide to resolve their conflict. Total testosterone 0.98 (0.77 to
1.27), with a mild directional-pleiotropy flag (Egger intercept P 0.04). Script:
`run_mr_multiexposure.py`.

### Multivariable MR (done)

SHBG and bioavailable testosterone in women, jointly, on aSAH (106 SNPs): both
null (SHBG 1.06 [0.50 to 2.28], bioavailable testosterone 1.03 [0.57 to 1.86]). The
SHBG point estimate moves from 0.73 (univariable) to 1.06 once testosterone is held
fixed, so the univariable SHBG signal was in part its shared testosterone. The data
cannot resolve the Molenberg 2022 vs Tan/Wu 2025 conflict. Script: `run_mvmr.py`.

### r^2 LD clumping (done)

Real reference-panel clumping (PLINK 2, 1000 Genomes European, r^2 < 0.001 over
10 Mb, with a chr:pos to rsID map for the menopause instruments) gives the same
answer as distance clumping: age at menopause on aSAH, IVW-RE OR 1.03 (0.98 to
1.09), 81 instruments. The distance-window sensitivity (250 kb to 5 Mb) is also
stable (1.02 to 1.03). Scripts: `run_mr_ld_clump.py`, `run_mr_clumping_sensitivity.py`.

## Run

```bash
uv sync --extra dev
uv run pytest                       # estimator + harmonization + QC tests (no data)
python scripts/run_mr_full.py       # ANM -> aSAH, full QC + power + equivalence
python scripts/run_mr_poscontrol.py # positive control: ANM -> breast cancer
```

GWAS summary statistics are not committed (large / restricted); `sources.py` lists
exact provenance, `docs/gwas_access.md` explains how to fetch them, and `data/gwas/`
is git-ignored.

## License

Code: MIT. Uses only public, published GWAS summary statistics per their terms.
