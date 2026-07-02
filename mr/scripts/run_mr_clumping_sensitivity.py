"""Extension 2: clumping sensitivity for the ANM and SHBG -> aSAH estimates.

Proper r^2 LD clumping needs a genotype reference (1000G) and rsIDs; the ReproGen
ANM file carries chr:pos identifiers only, so LD clumping would require an extra
position->rsID mapping step. As a feasible and directly informative alternative,
we vary the distance-clumping window from 250 kb to 5 Mb. A stable estimate across
windows indicates the result is not an artifact of clumping stringency, which is
the property r^2 clumping would otherwise be relied on to guarantee. The passing
positive control (menopause -> breast cancer) already shows distance clumping
yields valid instruments.
"""

from __future__ import annotations

import json
import math
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from estrogen_mr import qc
from estrogen_mr.instruments import distance_clump, harmonize_by_position

ROOT = Path(__file__).resolve().parents[1]
GWAS = ROOT / "data" / "gwas"
OUT = ROOT / "outputs"
WINDOWS = [250_000, 500_000, 1_000_000, 2_000_000, 5_000_000]

EXPOSURES = {
    "age_at_menopause": ("anm_ruth2021.txt.gz", False),
    "SHBG_women": ("shbg_female.tsv.gz", True),
}


def sig_rows(path: Path, catalog: bool) -> pd.DataFrame:
    if catalog:
        ren = {"chromosome": "CHR", "base_pair_location": "POS", "effect_allele": "Effect_Allele",
               "other_allele": "Other_Allele", "effect_allele_frequency": "EAF",
               "beta": "Effect", "standard_error": "SE"}
    else:
        ren = {}
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf:
        subprocess.run(["bash", "-c", f"gzcat '{path}' | awk 'NR==1 || ($9!=\"NA\" && $9+0<5e-8)'"],
                       stdout=tf, check=True)
        tmp = tf.name
    df = pd.read_csv(tmp, sep="\t").rename(columns=ren)
    Path(tmp).unlink(missing_ok=True)
    df["Pval"] = pd.to_numeric(df.get("Pval", df.get("p_value")), errors="coerce")
    return df.dropna(subset=["Effect", "SE", "CHR", "POS"])


def outcome():
    out = pd.read_csv(GWAS / "bakker_sah_euro_noUKBB.txt.gz", sep="\t")
    out["key"] = out["CHR"].astype(str) + ":" + out["BP"].astype(str)
    return out.rename(columns={"A_EFF": "oa_eff", "A_NONEFF": "oa_non",
                               "BETA": "beta_out", "SE": "se_out"})[
        ["key", "oa_eff", "oa_non", "beta_out", "se_out"]]


def main() -> None:
    OUT.mkdir(exist_ok=True)
    out = outcome()
    results = {}
    for name, (fname, catalog) in EXPOSURES.items():
        sig = sig_rows(GWAS / fname, catalog)
        rows = []
        for w in WINDOWS:
            inst = distance_clump(sig, window=w).copy()
            inst["key"] = inst["CHR"].astype(int).astype(str) + ":" + inst["POS"].astype(int).astype(str)
            h = harmonize_by_position(inst, out)
            bx, sx, by, sy = (h[c].values for c in ("bx", "sx", "by", "sy"))
            r = qc.ivw_random(bx, sx, by, sy)
            rows.append({"window_kb": w // 1000, "n": int(len(h)),
                         "or": r.odds_ratio, "ci_low": math.exp(r.ci_low),
                         "ci_high": math.exp(r.ci_high), "p": r.p})
        results[name] = rows

    (OUT / "mr_clumping_sensitivity.json").write_text(json.dumps(results, indent=2))
    for name, rows in results.items():
        print(f"\n{name}: IVW-RE OR across clumping windows")
        for r in rows:
            print(f"  {r['window_kb']:>5} kb  n={r['n']:>3d}  "
                  f"OR {r['or']:.3f} ({r['ci_low']:.3f}-{r['ci_high']:.3f})  p={r['p']:.3g}")


if __name__ == "__main__":
    main()
