"""Export per-SNP harmonized arrays for the MR sensitivity figures.

Writes mr/outputs/figure_data.json with bx/sx/by/sy (and leave-one-out ORs) for the
primary ANM -> aSAH analysis and the ANM -> breast cancer positive control. These
are the small summary-level arrays the scatter / funnel / leave-one-out plots need.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from estrogen_mr import qc
from estrogen_mr.instruments import distance_clump, harmonize_by_position

ROOT = Path(__file__).resolve().parents[1]
GWAS = ROOT / "data" / "gwas"
OUT = ROOT / "outputs"


def _anm_instruments():
    exp = pd.read_csv(GWAS / "anm_ruth2021.txt.gz", sep="\t")
    exp["Pval"] = pd.to_numeric(exp["Pval"], errors="coerce")
    exp = exp[exp["Pval"] < 5e-8].dropna(subset=["Pval", "Effect", "SE"])
    inst = distance_clump(exp).copy()
    inst["key"] = inst["CHR"].astype(str) + ":" + inst["POS"].astype(str)
    return inst


def _harm_to(inst, out_path, ren, chrom_col, pos_col):
    out = pd.read_csv(out_path, sep="\t")
    out["key"] = out[chrom_col].astype(str) + ":" + out[pos_col].astype(str)
    out = out.rename(columns=ren)
    return harmonize_by_position(inst, out[["key", "oa_eff", "oa_non", "beta_out", "se_out"]])


def main() -> None:
    OUT.mkdir(exist_ok=True)
    inst = _anm_instruments()

    asah = _harm_to(inst, GWAS / "bakker_sah_euro_noUKBB.txt.gz",
                    {"A_EFF": "oa_eff", "A_NONEFF": "oa_non", "BETA": "beta_out", "SE": "se_out"},
                    "CHR", "BP")
    bc = _harm_to(inst, GWAS / "bc_build37.tsv.gz",
                  {"effect_allele": "oa_eff", "other_allele": "oa_non",
                   "beta": "beta_out", "standard_error": "se_out"},
                  "chromosome", "base_pair_location")

    def pack(h, with_loo=True):
        d = {c: h[c].astype(float).tolist() for c in ("bx", "sx", "by", "sy")}
        if with_loo:
            d["loo_or"] = [x["or"] for x in
                          qc.leave_one_out(h["bx"].values, h["sx"].values,
                                           h["by"].values, h["sy"].values)]
        return d

    payload = {"primary_asah": pack(asah), "poscontrol_bc": pack(bc, with_loo=False)}
    (OUT / "figure_data.json").write_text(json.dumps(payload))
    print(f"wrote {OUT/'figure_data.json'}: aSAH {len(asah)} SNPs, BC {len(bc)} SNPs")


if __name__ == "__main__":
    main()
