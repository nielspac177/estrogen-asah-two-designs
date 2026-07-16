"""FinnGen replication of the ANM -> aSAH MR (independent Finnish cohort, no UKB
overlap). FinnGen is GRCh38, so we harmonize by rsID (build-independent) rather
than by chr:pos: the ANM instruments are mapped to rsIDs via the 1000G EUR bim and
r2-clumped (reusing outputs/_anm.clumps), and matched to FinnGen on rsID.

Two steps:
  1. `uv run python scripts/run_mr_finngen.py` with no extract present writes the
     clumped instrument table (outputs/anm_instruments_rsid.tsv) and the rsID list
     (outputs/_finngen_rsids.txt), then exits telling you to stream-extract.
  2. Stream FinnGen rows for those rsIDs into data/gwas/finngen_R11_I9_SAH_extract.tsv
     (see the printed curl command), then re-run to compute the MR.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from estrogen_mr import qc
from estrogen_mr.instruments import f_statistic, harmonize_by_position
from estrogen_mr.methods import ivw, mr_egger, weighted_median

ROOT = Path(__file__).resolve().parents[1]
GWAS = ROOT / "data" / "gwas"
LD = ROOT / "data" / "ldref"
OUT = ROOT / "outputs"
N_EXP = 201323
ANM_SD_YEARS = 4.0
EXTRACT = GWAS / "finngen_R11_I9_SAH_extract.tsv"
FINNGEN_URL = ("https://storage.googleapis.com/finngen-public-data-r11/"
               "summary_stats/finngen_R11_I9_SAH.gz")
# FinnGen R11 I9_SAH: 3,151 cases / 393,987 controls (effective N approx 4kk/(k+k))
FINNGEN_N_CASES, FINNGEN_N_CTRL = 3151, 393987


def bim_map() -> dict:
    m = {}
    with open(LD / "EUR.bim") as fh:
        for line in fh:
            c = line.split()
            m[f"{c[0]}:{c[3]}"] = c[1]
    return m


def build_instruments() -> pd.DataFrame:
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf:
        subprocess.run(["bash", "-c",
            f"gzcat '{GWAS/'anm_ruth2021.txt.gz'}' | awk 'NR==1 || ($9!=\"NA\" && $9+0<5e-8)'"],
            stdout=tf, check=True)
        sig = pd.read_csv(tf.name, sep="\t")
    sig["key"] = sig["CHR"].astype(str) + ":" + sig["POS"].astype(str)
    sig["rsid"] = sig["key"].map(bim_map())
    have = sig.dropna(subset=["rsid"]).drop_duplicates("rsid")
    clumps = pd.read_csv(OUT / "_anm.clumps", sep=r"\s+")
    kept = set(clumps["ID"] if "ID" in clumps else clumps["SNP"])
    inst = have[have["rsid"].isin(kept)].copy()
    return inst[["rsid", "CHR", "POS", "Effect_Allele", "Other_Allele", "EAF", "Effect", "SE"]]


def load_extract(inst_rsids: set) -> pd.DataFrame:
    # FinnGen columns: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt ...
    df = pd.read_csv(EXTRACT, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    df = df.rename(columns={"rsids": "key", "alt": "oa_eff", "ref": "oa_non",
                            "beta": "beta_out", "sebeta": "se_out"})
    return df[df["key"].isin(inst_rsids)][["key", "oa_eff", "oa_non", "beta_out", "se_out"]]


def main() -> None:
    OUT.mkdir(exist_ok=True)
    inst = build_instruments().copy()
    inst["key"] = inst["rsid"]
    inst.to_csv(OUT / "anm_instruments_rsid.tsv", sep="\t", index=False)
    (OUT / "_finngen_rsids.txt").write_text("\n".join(inst["rsid"]) + "\n")
    print(f"clumped ANM instruments with rsIDs: {len(inst)}")

    if not EXTRACT.exists():
        print("\nFinnGen extract not found. Stream-extract the instrument rows with:\n")
        print(f"  (gzcat header + matched rows)\n"
              f"  curl -sL '{FINNGEN_URL}' | gunzip \\\n"
              f"    | awk 'NR==1 || FNR==NR' ...  # see run note; simplest below")
        print(f"\n  H=$(curl -sL '{FINNGEN_URL}' | gunzip 2>/dev/null | head -1);\n"
              f"  {{ echo \"$H\"; curl -sL '{FINNGEN_URL}' | gunzip \\\n"
              f"      | grep -F -w -f '{OUT/'_finngen_rsids.txt'}'; }} > '{EXTRACT}'")
        return

    out = load_extract(set(inst["rsid"]))
    n_out = int(4 * FINNGEN_N_CASES * FINNGEN_N_CTRL / (FINNGEN_N_CASES + FINNGEN_N_CTRL))
    h = harmonize_by_position(inst, out)
    bx, sx, by, sy = (h[c].values for c in ("bx", "sx", "by", "sy"))

    ivw_re = qc.ivw_random(bx, sx, by, sy)
    egg, egg_i = mr_egger(bx, sx, by, sy)
    wm = weighted_median(bx, sx, by, sy, n_boot=1000)
    st = qc.steiger(bx, sx, by, sy, N_EXP, n_out)
    loo = [x["or"] for x in qc.leave_one_out(bx, sx, by, sy)]
    eq = qc.tost(ivw_re.estimate * ANM_SD_YEARS, ivw_re.se * ANM_SD_YEARS, sesoi_or=0.90)

    res = {
        "outcome_label": "FinnGen R11 nontraumatic SAH (I9_SAH), Finnish, independent",
        "exposure": "age_at_natural_menopause_Ruth2021",
        "harmonized_by": "rsID (FinnGen GRCh38; build-independent match)",
        "n_cases": FINNGEN_N_CASES, "n_controls": FINNGEN_N_CTRL, "outcome_Neff": n_out,
        "n_harmonized": int(len(h)), "mean_F": round(f_statistic(bx, sx), 1),
        "primary_IVW_random": ivw_re.as_dict(), "IVW_fixed": ivw(bx, sx, by, sy).as_dict(),
        "MR_Egger": egg.as_dict(), "egger_intercept": egg_i,
        "weighted_median": wm.as_dict(),
        "steiger": {k: v for k, v in st.items() if k != "keep_mask"},
        "leave_one_out_or_range": [round(min(loo), 3), round(max(loo), 3)],
        "equivalence_TOST": eq,
    }
    (OUT / "mr_finngen_SAH.json").write_text(json.dumps(res, indent=2))
    v = ivw_re.as_dict()
    print(f"\nFinnGen SAH replication: harmonized {len(h)}, F={res['mean_F']}, Neff={n_out}")
    print(f"  IVW-RE OR/yr {v['or']:.3f} ({v['or_ci_low']:.3f}-{v['or_ci_high']:.3f}) p={v['p']:.3g}")
    print(f"  Egger {res['MR_Egger']['or']:.3f} (int p={egg_i['intercept_p']:.3g})  "
          f"WM {res['weighted_median']['or']:.3f}  LOO {res['leave_one_out_or_range']}")
    print("saved outputs/mr_finngen_SAH.json")


if __name__ == "__main__":
    main()
