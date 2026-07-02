"""Gold-standard refinement 1: multivariable MR of SHBG + bioavailable testosterone
(both in women, Ruth 2020) on aSAH, to separate the two overlapping instruments.

Single-exposure MR cannot tell whether an SHBG signal is really SHBG or the
bioavailable-testosterone it partly determines; that ambiguity is the likely
source of the Molenberg 2022 vs Tan/Wu 2025 disagreement. MVMR estimates each
direct effect holding the other fixed.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from estrogen_mr.instruments import distance_clump
from estrogen_mr.mvmr import mvmr_ivw

ROOT = Path(__file__).resolve().parents[1]
GWAS = ROOT / "data" / "gwas"
OUT = ROOT / "outputs"
SHBG = GWAS / "shbg_female.tsv.gz"
BIOT = GWAS / "biot_female.tsv.gz"
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _awk_sig(path):
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf:
        subprocess.run(["bash", "-c", f"gzcat '{path}' | awk 'NR==1 || ($9!=\"NA\" && $9+0<5e-8)'"],
                       stdout=tf, check=True)
        return tf.name


def _load_catalog(path):
    df = pd.read_csv(path, sep="\t").rename(columns={
        "chromosome": "CHR", "base_pair_location": "POS", "effect_allele": "EA",
        "other_allele": "OA", "beta": "B", "standard_error": "SE", "p_value": "Pval"})
    df["key"] = df["CHR"].astype(str) + ":" + df["POS"].astype(str)
    return df


def _extract_by_keys(path, keys):
    """awk full GWAS-catalog file for rows whose CHR:POS is in the key set."""
    kf = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
    kf.write("\n".join(keys))
    kf.close()
    of = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
    subprocess.run(["bash", "-c",
        f"gzcat '{path}' | awk 'FNR==NR{{k[$1];next}} FNR==1{{print;next}} (($2\":\"$3) in k)' "
        f"'{kf.name}' -"], stdout=of, check=True)
    of.close()
    return _load_catalog(of.name)


def _align(ea_ref, oa_ref, ea, oa, beta):
    ea_ref, oa_ref, ea, oa = ea_ref.upper(), oa_ref.upper(), str(ea).upper(), str(oa).upper()
    if {ea, oa} == {ea_ref, oa_ref}:
        return -beta if (ea, oa) == (oa_ref, ea_ref) else beta
    if {_COMP.get(ea), _COMP.get(oa)} == {ea_ref, oa_ref}:
        return -beta if (_COMP[ea], _COMP[oa]) == (oa_ref, ea_ref) else beta
    return None  # incompatible


def main() -> None:
    OUT.mkdir(exist_ok=True)
    sig = pd.concat([_load_catalog(_awk_sig(SHBG)), _load_catalog(_awk_sig(BIOT))], ignore_index=True)
    sig = sig.dropna(subset=["B", "SE", "CHR", "POS"])
    inst = distance_clump(sig).copy()
    inst["key"] = inst["CHR"].astype(int).astype(str) + ":" + inst["POS"].astype(int).astype(str)
    keys = set(inst["key"])

    shbg = _extract_by_keys(SHBG, keys).set_index("key")
    biot = _extract_by_keys(BIOT, keys).set_index("key")
    out = pd.read_csv(GWAS / "bakker_sah_euro_noUKBB.txt.gz", sep="\t")
    out["key"] = out["CHR"].astype(str) + ":" + out["BP"].astype(str)
    out = out.set_index("key")

    rows = []
    for k in keys:
        if k not in shbg.index or k not in biot.index or k not in out.index:
            continue
        s, b, o = shbg.loc[k], biot.loc[k], out.loc[k]
        if isinstance(s, pd.DataFrame) or isinstance(b, pd.DataFrame) or isinstance(o, pd.DataFrame):
            continue  # multiallelic collisions
        bb = _align(s["EA"], s["OA"], b["EA"], b["OA"], b["B"])          # bioT -> SHBG allele
        bo = _align(s["EA"], s["OA"], o["A_EFF"], o["A_NONEFF"], o["BETA"])  # outcome -> SHBG allele
        if bb is None or bo is None:
            continue
        rows.append({"b_shbg": s["B"], "b_biot": bb, "b_out": bo, "se_out": o["SE"]})
    d = pd.DataFrame(rows).dropna()

    res = mvmr_ivw(d[["b_shbg", "b_biot"]].values, d["b_out"].values, d["se_out"].values,
                   ["SHBG_women", "bioavailable_testosterone_women"])
    payload = res.as_dict()
    payload["note"] = ("MVMR direct effects on aSAH; per SD of each hormone. Compares to "
                       "Molenberg 2022 (higher SHBG -> more aSAH) vs Tan/Wu 2025 (opposite).")
    (OUT / "mr_mvmr_shbg_biot.json").write_text(json.dumps(payload, indent=2))

    print(f"MVMR (SHBG + bioavailable testosterone, women) -> aSAH, n_snps={res.n_snps}")
    for nm, v in payload["exposures"].items():
        print(f"  {nm:34s} OR {v['or']:.2f} ({v['or_ci_low']:.2f}-{v['or_ci_high']:.2f})  p={v['p']:.3g}")


if __name__ == "__main__":
    main()
