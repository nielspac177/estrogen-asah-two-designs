"""eFigure. Why the observational menopause contrast is not identifiable.

(A) Age distribution by group: premenopausal and postmenopausal women occupy
    disjoint age ranges (zero common support), while men span the whole range and
    supply the ageing reference for the difference-in-differences. A positivity
    violation, not a collinearity problem (VIF is only ~2.3).
(B) Delayed cerebral ischaemia rate by age in women and men: DCI falls with age in
    BOTH sexes with no menopause-specific step at 51, so the apparent postmenopausal
    "protection" (OR 0.52) is the shared age gradient; the age x sex difference-in-
    differences that removes it is null (OR-ratio 1.04).

Reads the aggregate exports icu/outputs/{age_overlap.csv, age_dci_by_sex.csv, vif.json}.
"""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": ["Helvetica Neue", "Arial"], "font.size": 9.5})
INK = "#141414"; PRE = "#4a4a4a"; POST = "#a8a8a8"; MEN = "#141414"; ACC = "#b22222"

ROOT = Path(__file__).resolve().parents[0].parent
OUT = ROOT / "icu" / "outputs"
ov = pd.read_csv(OUT / "age_overlap.csv")
dci = pd.read_csv(OUT / "age_dci_by_sex.csv")
vif = json.loads((OUT / "vif.json").read_text())

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.4, 4.7))

# ---- Panel A: age distribution by group ----
x = ov["age_lo"].values
w = 4.6
axA.bar(x, ov["premenopausal"], width=w, align="edge", color=PRE, label="Premenopausal women")
axA.bar(x, ov["postmenopausal"], width=w, align="edge", color=POST, label="Postmenopausal women")
axA.step(np.append(x, x[-1] + 5), np.append(ov["men"], ov["men"].iloc[-1]),
         where="post", color=MEN, lw=1.4, label="Men (ageing reference)")
axA.axvline(51, color=ACC, ls="--", lw=1.3)
axA.text(51.6, axA.get_ylim()[1] * 0.96, "cutoff = 51 y\n(= menopause definition)",
         color=ACC, fontsize=8, va="top")
axA.set_xlabel("Age, y"); axA.set_ylabel("Patients")
axA.set_title("A  Zero common support across strata", loc="left", fontsize=10.5, fontweight="bold")
axA.legend(frameon=False, fontsize=8, loc="upper left")
axA.text(0.02, 0.60,
         "Premenopausal 19-50 y; postmenopausal 51-91 y.\n"
         "No woman appears on both sides of the cutoff:\n"
         "a positivity violation, so age cannot be adjusted\n"
         f"for (VIF only {vif['model_with_age']['VIF']['post']} - not collinearity).",
         transform=axA.transAxes, fontsize=7.8, color=INK, va="top",
         bbox=dict(facecolor="white", alpha=0.82, edgecolor="none", pad=1.5))
for s in ("top", "right"):
    axA.spines[s].set_visible(False)

# ---- Panel B: DCI rate vs age, by sex ----
for lab, mk, col in [("women", "o", INK), ("men", "s", "#9a9a9a")]:
    d = dci[dci["sexlab"] == lab].sort_values("age_lo")
    axB.plot(d["age_lo"] + 2.5, d["rate"] * 100, marker=mk, color=col, lw=1.6,
             ms=4, label=lab.capitalize())
axB.axvline(51, color=ACC, ls="--", lw=1.3)
axB.set_xlabel("Age, y"); axB.set_ylabel("DCI rate, %")
axB.set_title("B  Age gradient in both sexes, no step at menopause", loc="left",
              fontsize=10.5, fontweight="bold")
axB.legend(frameon=False, fontsize=8.5, loc="upper right")
axB.text(0.02, 0.16,
         "DCI falls with age in both sexes with no extra drop\n"
         "at 51, so the postmenopausal 'protection' (OR 0.52)\n"
         "is the age gradient; the age x sex difference-in-\n"
         "differences that removes it is null (1.04, 0.62-1.76).",
         transform=axB.transAxes, fontsize=7.8, color=INK, va="bottom",
         bbox=dict(facecolor="white", alpha=0.82, edgecolor="none", pad=1.5))
for s in ("top", "right"):
    axB.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig(ROOT / "figures" / "overlap_nonidentifiability.png", dpi=200,
            bbox_inches="tight", facecolor="white")
print("saved figures/overlap_nonidentifiability.png")
