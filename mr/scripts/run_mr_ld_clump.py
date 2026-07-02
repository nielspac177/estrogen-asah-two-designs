"""Gold-standard refinement 2: real r^2 LD clumping (PLINK + 1000G EUR).

The menopause instruments are chr:pos only, so we map them to rsIDs via the 1000G
bim, LD-clump at r^2<0.001 / 10 Mb with PLINK against the EUR reference, map the
clumped index SNPs back to chr:pos, then re-run ANM -> aSAH. Compares directly
with the distance-clumped primary.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from estrogen_mr import qc
from estrogen_mr.instruments import f_statistic, harmonize_by_position

ROOT = Path(__file__).resolve().parents[1]
GWAS = ROOT / "data" / "gwas"
LD = ROOT / "data" / "ldref"
OUT = ROOT / "outputs"
PLINK = LD / "plink2"


def bim_map() -> dict:
    """chr:pos -> rsID from the 1000G EUR bim (cols: chr, rsid, cM, pos, a1, a2)."""
    m = {}
    with open(LD / "EUR.bim") as fh:
        for line in fh:
            c = line.split()
            m[f"{c[0]}:{c[3]}"] = c[1]
    return m


def main() -> None:
    OUT.mkdir(exist_ok=True)
    # sig ANM instruments
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf:
        subprocess.run(["bash", "-c",
            f"gzcat '{GWAS/'anm_ruth2021.txt.gz'}' | awk 'NR==1 || ($9!=\"NA\" && $9+0<5e-8)'"],
            stdout=tf, check=True)
        sig = pd.read_csv(tf.name, sep="\t")
    sig["key"] = sig["CHR"].astype(str) + ":" + sig["POS"].astype(str)
    sig["Pval"] = pd.to_numeric(sig["Pval"], errors="coerce")

    m = bim_map()
    sig["rsid"] = sig["key"].map(m)
    have = sig.dropna(subset=["rsid"]).drop_duplicates("rsid")
    print(f"ANM sig SNPs: {len(sig)}; found in 1000G EUR: {len(have)}")

    clump_in = OUT / "_anm_clump_in.tsv"
    have[["rsid", "Pval"]].rename(columns={"rsid": "SNP", "Pval": "P"}).to_csv(
        clump_in, sep="\t", index=False)

    subprocess.run([str(PLINK), "--bfile", str(LD / "EUR"), "--clump", str(clump_in),
                    "--clump-p1", "5e-8", "--clump-r2", "0.001", "--clump-kb", "10000",
                    "--clump-id-field", "SNP", "--clump-p-field", "P",
                    "--out", str(OUT / "_anm")], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    clumps = pd.read_csv(OUT / "_anm.clumps", sep=r"\s+")
    kept_rsid = set(clumps["ID"] if "ID" in clumps else clumps["SNP"])
    inst = have[have["rsid"].isin(kept_rsid)].copy()
    inst["key"] = inst["CHR"].astype(int).astype(str) + ":" + inst["POS"].astype(int).astype(str)
    print(f"instruments after r^2 clumping: {len(inst)}")

    out = pd.read_csv(GWAS / "bakker_sah_euro_noUKBB.txt.gz", sep="\t")
    out["key"] = out["CHR"].astype(str) + ":" + out["BP"].astype(str)
    out = out.rename(columns={"A_EFF": "oa_eff", "A_NONEFF": "oa_non",
                              "BETA": "beta_out", "SE": "se_out"})
    h = harmonize_by_position(inst, out[["key", "oa_eff", "oa_non", "beta_out", "se_out"]])
    bx, sx, by, sy = (h[c].values for c in ("bx", "sx", "by", "sy"))
    re = qc.ivw_random(bx, sx, by, sy)

    res = {"method": "r2 LD clumping (PLINK, 1000G EUR, r2<0.001/10Mb)",
           "n_in_1000G": int(len(have)), "n_clumped": int(len(inst)), "n_harmonized": int(len(h)),
           "mean_F": round(f_statistic(bx, sx), 1), "IVW_random": re.as_dict()}
    (OUT / "mr_ld_clumped.json").write_text(json.dumps(res, indent=2))
    v = re.as_dict()
    print(f"ANM -> aSAH, r^2 LD-clumped: IVW-RE OR {v['or']:.3f} "
          f"({v['or_ci_low']:.3f}-{v['or_ci_high']:.3f}), p={v['p']:.3g}, n={len(h)}, F={res['mean_F']}")
    print("  (distance-clumped primary was OR 1.03 (0.97-1.09) with 85 instruments)")


if __name__ == "__main__":
    main()
