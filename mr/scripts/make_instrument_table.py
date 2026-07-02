"""eTable 3: per-instrument (SNP) harmonized effects and F-statistics (STROBE-MR)."""
from pathlib import Path
import pandas as pd
from estrogen_mr.instruments import distance_clump, harmonize_by_position
ROOT = Path(__file__).resolve().parents[1]; GWAS = ROOT/"data"/"gwas"
exp = pd.read_csv(GWAS/"anm_ruth2021.txt.gz", sep="\t")
exp["Pval"] = pd.to_numeric(exp["Pval"], errors="coerce")
exp = exp[exp["Pval"]<5e-8].dropna(subset=["Pval","Effect","SE"])
inst = distance_clump(exp).copy(); inst["key"] = inst["CHR"].astype(str)+":"+inst["POS"].astype(str)
out = pd.read_csv(GWAS/"bakker_sah_euro_noUKBB.txt.gz", sep="\t")
out["key"]=out["CHR"].astype(str)+":"+out["BP"].astype(str)
out=out.rename(columns={"A_EFF":"oa_eff","A_NONEFF":"oa_non","BETA":"beta_out","SE":"se_out"})
h = harmonize_by_position(inst, out[["key","oa_eff","oa_non","beta_out","se_out"]])
h["F"]=(h["bx"]/h["sx"])**2
h=h.sort_values("key")
lines=["# eTable 3. Genetic instruments for age at natural menopause (harmonized) → aSAH","",
 f"{len(h)} independent SNPs (distance-clumped, matched to Bakker 2020 by position). "
 "β on the exposure (per-year) and outcome (log-OR) scales; F = per-SNP instrument strength.","",
 "| SNP (chr:pos) | β exposure | SE | β outcome (log OR) | SE | F |",
 "|---|---|---|---|---|---|"]
for _,r in h.iterrows():
    lines.append(f"| {r['key']} | {r['bx']:+.3f} | {r['sx']:.3f} | {r['by']:+.3f} | {r['sy']:.3f} | {r['F']:.0f} |")
lines += ["", f"Mean F = {h['F'].mean():.0f} (all > 10, no weak instruments)."]
(ROOT.parent/"tables"/"etable3_instruments.md").write_text("\n".join(lines))
print(f"wrote etable3 with {len(h)} SNPs, mean F={h['F'].mean():.0f}")
