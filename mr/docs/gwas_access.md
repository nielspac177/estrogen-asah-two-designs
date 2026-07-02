# How to access the GWAS summary statistics (for another agent)

Two-sample MR needs only **public summary statistics**, no individual data, no
PLINK, no special genomics libraries. Everything below was done with `curl` +
Python (pandas). All files here are **GRCh37/hg19**.

## Outcome, intracranial aneurysm / aSAH (Bakker et al. 2020, Nat Genet)

Hosted on **Figshare**, public, no auth. Use the Figshare API to get download URLs:

```bash
curl -s "https://api.figshare.com/v2/articles/11303372" \
  | python3 -c "import sys,json; [print(f['size'],f['name'],f['download_url']) for f in json.load(sys.stdin)['files']]"
# download_url looks like https://ndownloader.figshare.com/files/<id>
curl -sL "https://ndownloader.figshare.com/files/29146425" -o bakker_sah_euro_noUKBB.txt.gz
```

- Files (~75–95 MB gz each): Stage_1/2, SAH-only, uIA-only, EastAsian, and
  **`*_excludingUKBB`** variants (Stage_1/2, SAH-only European & cross-ancestry).
- Columns: `CHR BP SNP A_EFF A_NONEFF Freq_EFF BETA SE P MAF A_MINOR N`.
  `N` is **effective** N = 4·Ncase·Ncontrol/(Ncase+Ncontrol) (max ≈17k for aSAH-Euro).
- **Pick the UKB-excluded file** if your exposure GWAS used UK Biobank (see overlap).

## Exposures, reproductive timing (ReproGen)

Direct download from reprogen.org. Find the link by scraping the page:

```bash
curl -sL "https://www.reprogen.org/data_download.html" | grep -ioE 'href="[^"]*(gz|zip)"'
curl -sL "https://www.reprogen.org/reprogen_ANM_201K_170621.txt.gz" -o anm_ruth2021.txt.gz  # age at natural menopause, 248 MB
# Menarche (Day 2017): Menarche_1KG_NatGen2017_WebsiteUpload.zip
```

- ANM columns: `SNP CHR POS Effect_Allele Other_Allele EAF Effect SE Pval N`.
- Gotcha: ANM `SNP` IDs are `chr:pos[:type]` (not rsIDs), so **match to the outcome
  by CHR:POS**, not rsID. INDEL rows exist.

## Exposures, sex hormones (Ruth et al. 2020, Nat Med; SHBG, testosterone)

**GWAS Catalog** (EBI FTP), accessions GCST90012103 (bioT), GCST90012109/…111 (SHBG),
GCST90012113 (total T). Harmonised TSVs under
`http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90012001-GCST90013000/<ACC>/harmonised/`.

## Instrument selection & matching (summary-stats only)

1. Filter exposure to genome-wide significant SNPs (`Pval < 5e-8`), do it with
   `awk` before loading (files are 10M+ rows).
2. **Clumping:** no PLINK here → **distance clumping** (keep lowest-p SNP, drop
   others within ±1 Mb/chr) as an approximation. Proper r²-clumping needs a 1000G
   EUR panel (downloadable; LD computable in Python) or an API (OpenGWAS/ieugwasr
   now needs a token; LDlink is free).
3. **Match** exposure↔outcome on `CHR:POS` (both GRCh37).
4. **Harmonize alleles:** align outcome effect to the exposure effect allele (flip
   `beta_out` sign on allele swap; handle strand via complement); drop ambiguous
   palindromes (A/T, C/G) at intermediate frequency using EAF.

## Two big gotchas

- **Sample overlap:** UKB-derived exposures (Ruth ANM 2021, Day menarche, Ruth
  SHBG/T 2020 all include UKB) → **use the UKB-excluded outcome** to avoid
  overlap bias (which pulls two-sample MR toward the confounded estimate).
- **Genome build:** confirm both are GRCh37. If allele-concordant CHR:POS matches
  are near zero, suspect an hg19/hg38 mismatch (liftover with `pyliftover`).
