"""eFigure. Observational-arm diagnostics: specification curve and covariate balance.

Left: the menopause->DCI odds ratio across all 72 defensible model specifications;
the point of the panel is that significant fits all fall on the anti-protective
side, the signature of age confounding. Right: standardized mean differences before
and after inverse-probability weighting. Reads icu/outputs/{spec_curve,balance}.csv.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": ["Avenir Next", "Arial"], "font.size": 9.5})
BLUE = "#0B6FB8"; GREY = "#8a97a3"; RED = "#c0392b"; INK = "#1b2733"

ROOT = Path(__file__).resolve().parents[1]
sc = pd.read_csv(ROOT / "icu" / "outputs" / "spec_curve.csv")
bal = pd.read_csv(ROOT / "icu" / "outputs" / "balance.csv")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6),
                               gridspec_kw={"width_ratios": [1.5, 1]})

# ---- A. specification curve ----
sc = sc[sc["converged"]].sort_values("or").reset_index(drop=True)
x = np.arange(len(sc))
sig = sc["p"] < 0.05
ax1.scatter(x[~sig], sc["or"][~sig], s=16, color=GREY, label="not significant")
ax1.scatter(x[sig], sc["or"][sig], s=18, color=BLUE, label="P < .05")
ax1.axhline(1.0, color=RED, ls="--", lw=1.2)
prim = sc[(sc.cutoff == 51) & (sc.adjust == "clinical") & (sc.outcome == "dci_composite")
          & (sc.source == "all")]
if len(prim):
    ax1.scatter([sc.index[sc["or"] == prim["or"].iloc[0]][0]], [prim["or"].iloc[0]],
                s=70, facecolor="none", edgecolor=INK, lw=1.8, label="primary")
ax1.set_yscale("log")
ax1.set_yticks([0.3, 0.5, 0.75, 1.0, 1.5])
ax1.set_yticklabels(["0.3", "0.5", "0.75", "1.0", "1.5"])
ax1.set_xlabel("Model specification (ranked by odds ratio)")
ax1.set_ylabel("OR, postmenopausal vs premenopausal")
ax1.set_title("A. Specification curve (72 models)", loc="left", fontweight="bold", fontsize=10.5)
ax1.text(2, 1.32, f"{int(sig.sum())} of {len(sc)} significant,\nall below OR = 1",
         fontsize=8.5, color=INK)
ax1.legend(fontsize=8, frameon=False, loc="lower right")

# ---- B. covariate balance (love plot) ----
y = np.arange(len(bal))
ax2.scatter(bal["smd_unweighted"], y, s=42, color=GREY, label="unweighted")
ax2.scatter(bal["smd_weighted"], y, s=42, color=BLUE, label="IPW-weighted")
ax2.axvline(0.1, color=RED, ls="--", lw=1, label="|SMD| = 0.1")
ax2.set_yticks(y); ax2.set_yticklabels(bal["covariate"])
ax2.set_xlabel("|Standardized mean difference|")
ax2.set_title("B. Covariate balance", loc="left", fontweight="bold", fontsize=10.5)
ax2.legend(fontsize=8, frameon=False)
ax2.margins(y=0.2)

fig.tight_layout()
fig.savefig(ROOT / "figures" / "observational_diagnostics.png", dpi=180,
            bbox_inches="tight", facecolor="white")
print("saved")
