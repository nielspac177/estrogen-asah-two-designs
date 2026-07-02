"""Table 1. Baseline characteristics, rendered as a JAMA-style table figure.

Same visual language as Figure 3: horizontal rules only, monochrome, shaded section
bands, footnotes. Reads icu/outputs/table1.csv.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({"font.family": ["Helvetica Neue", "Arial"], "font.size": 9.5})
INK = "#141414"; RULEC = "#141414"; SHADE = "#f0f0f0"

ROOT = Path(__file__).resolve().parents[1]
t = pd.read_csv(ROOT / "icu" / "outputs" / "table1.csv", index_col=0)


def v(row, col, fmt="{:.1f}"):
    return fmt.format(t.loc[row, col])


# section, label, csv-row, format
SPEC = [
    ("H", "Demographics", None),
    ("R", "Age, mean, y", "age_mean"),
    ("H", "Comorbidities, %", None),
    ("R", "Hypertension", "htn_pct"),
    ("R", "Current or former smoking", "smoking_pct"),
    ("R", "Diabetes", "diabetes_pct"),
    ("H", "Illness severity", None),
    ("R", "APACHE score, mean$^{a}$", "apache_mean"),
    ("H", "Outcomes, %", None),
    ("R", "DCI or vasospasm", "dci_pct"),
    ("R", "In-hospital mortality", "mortality_pct"),
]

# columns: (csv-col, header lines, x)
COLS = [
    ("overall", ["Overall", "(n = 1771)"], 45),
    ("premenopausal", ["Premenopausal", "women (n = 313)"], 62),
    ("postmenopausal", ["Postmenopausal", "women (n = 792)"], 79),
    ("male", ["Men", "(n = 666)"], 95),
]
X_LBL = 1.5

n = len(SPEC)
fig, ax = plt.subplots(figsize=(9.4, 0.52 * n + 2.4))
ax.axis("off"); ax.set_xlim(0, 100)
top = n + 1.0
ax.set_ylim(-2.2, top + 1.8)

# top rule + headers + header rule
ax.plot([0, 100], [top + 1.4, top + 1.4], color=RULEC, lw=1.6)
ax.text(X_LBL, top + 0.5, "Characteristic", ha="left", va="center", fontsize=9.5,
        fontweight="bold", color=INK)
for _, lines, x in COLS:
    ax.text(x, top + 0.9, lines[0], ha="center", va="center", fontsize=9.3, fontweight="bold", color=INK)
    ax.text(x, top + 0.0, lines[1], ha="center", va="center", fontsize=8.3, color=INK)
ax.plot([0, 100], [top - 0.55, top - 0.55], color=RULEC, lw=1.0)

y = top - 1.5
for kind, label, row in SPEC:
    if kind == "H":
        ax.add_patch(Rectangle((0, y - 0.5), 100, 1.0, facecolor=SHADE, edgecolor="none", zorder=0))
        ax.text(X_LBL, y, label, ha="left", va="center", fontsize=9.3, fontweight="bold", color=INK)
    else:
        ax.text(X_LBL + 1.2, y, label, ha="left", va="center", fontsize=9.2, color=INK)
        for col, _, x in COLS:
            ax.text(x, y, v(row, col), ha="center", va="center", fontsize=9.0, color=INK)
    y -= 1.0

bottom = y + 0.5
ax.plot([0, 100], [bottom, bottom], color=RULEC, lw=1.6)

fns = [
    "Abbreviations: APACHE, Acute Physiology and Chronic Health Evaluation; DCI, delayed cerebral ischaemia.",
    "Comorbidities are ascertained from diagnosis codes. $^{a}$APACHE score is recorded in eICU only.",
    "Menopausal strata use an age cutoff of 51 years, the population median age at natural menopause.",
]
for i, s in enumerate(fns):
    ax.text(0, bottom - 0.9 - i * 0.62, s, ha="left", va="top", fontsize=7.2, color="#333")

ax.text(0, top + 2.5, "Table 1. Baseline characteristics of the aSAH cohort, by menopausal stratum",
        ha="left", fontsize=11, fontweight="bold", color=INK)
fig.savefig(ROOT / "figures" / "table1.png", dpi=200, bbox_inches="tight", facecolor="white")
print("saved")
