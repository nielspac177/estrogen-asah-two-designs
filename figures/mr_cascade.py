"""Formation-to-rupture cascade + instrument controls, JAMA table-with-forest style.

Later age at natural menopause instrumented against unruptured IA (formation), aSAH
(rupture), and combined intracranial aneurysm; independent FinnGen replication; and
the positive (breast cancer) and negative (appendicitis) instrument controls. Same
visual language as results_forest_table.py. Data-driven: reads mr/outputs JSONs so
the numbers cannot go stale.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({"font.family": ["Helvetica Neue", "Arial"], "font.size": 9.5})
INK = "#141414"; MARK = "#222222"; RULEC = "#141414"; GREY = "#8a8a8a"; SHADE = "#f0f0f0"

ROOT = Path(__file__).resolve().parents[1]
MO = ROOT / "mr" / "outputs"
strat = json.loads((MO / "mr_outcomes_stratified.json").read_text())
finn = json.loads((MO / "mr_finngen_SAH.json").read_text())
pos = json.loads((MO / "mr_poscontrol_breastcancer.json").read_text())
neg = json.loads((MO / "mr_negcontrol_appendicitis.json").read_text())


def row(label, block, key="primary_IVW_random", fn=""):
    r = block[key]
    return ("R", label, str(block["n_harmonized"]),
            r["or"], r["or_ci_low"], r["or_ci_high"], r["p"], fn)


def ptxt(p):
    if p < 0.001:
        return "<.001"
    return f"{p:.2f}".lstrip("0")


ROWS = [
    ("H", "Disease course (age at natural menopause)", "", None, None, None, None, ""),
    row("Unruptured IA (formation)", strat["uIA_euro"], fn="a"),
    row("Aneurysmal SAH (rupture)", strat["aSAH_euro_noUKBB"], fn="a"),
    row("Combined intracranial aneurysm", strat["IA_combined_euro_noUKBB"], fn="a"),
    ("H", "Independent replication", "", None, None, None, None, ""),
    row("Aneurysmal SAH, FinnGen", finn, fn="b"),
    ("H", "Instrument controls", "", None, None, None, None, ""),
    row("Breast cancer (positive control)", pos, key="IVW", fn="c"),
    row("Appendicitis (negative control)", neg, key="IVW", fn="d"),
]

X_LBL, X_N = 1.5, 45
FX_L, FX_R = 46, 70
OR_LO, OR_HI = 0.4, 2.6
X_CI, X_P = 72.5, 99
TICKS = [0.5, 0.75, 1.0, 1.5, 2.0]


def fx(orr):
    return FX_L + (np.log(orr) - np.log(OR_LO)) / (np.log(OR_HI) - np.log(OR_LO)) * (FX_R - FX_L)


n = len(ROWS)
fig, ax = plt.subplots(figsize=(9.4, 0.52 * n + 2.4))
ax.axis("off"); ax.set_xlim(0, 100)
top = n + 1.0
ax.set_ylim(-3.6, top + 1.4)

ax.plot([0, 100], [top + 0.9, top + 0.9], color=RULEC, lw=1.6)
for x, t, ha in [(X_LBL, "Analysis", "left"), (X_N, "No.", "right"),
                 (X_CI, "OR (95% CI)", "left"), (X_P, "P Value", "right")]:
    ax.text(x, top + 0.1, t, ha=ha, va="center", fontsize=9.5, fontweight="bold", color=INK)
ax.plot([0, 100], [top - 0.55, top - 0.55], color=RULEC, lw=1.0)

ax.plot([fx(1.0), fx(1.0)], [-0.4, top - 0.7], color=GREY, lw=0.9, zorder=0)

y = top - 1.5
for kind, label, ntxt, orr, lo, hi, p, fn in ROWS:
    if kind == "H":
        ax.add_patch(Rectangle((0, y - 0.5), 100, 1.0, facecolor=SHADE, edgecolor="none", zorder=0))
        ax.text(X_LBL, y, label, ha="left", va="center", fontsize=9.3, fontweight="bold", color=INK)
    else:
        lab = label + (f"$^{{{fn}}}$" if fn else "")
        ax.text(X_LBL + 1.2, y, lab, ha="left", va="center", fontsize=9.2, color=INK)
        ax.text(X_N, y, ntxt, ha="right", va="center", fontsize=9.0, color=INK)
        ax.plot([fx(lo), fx(hi)], [y, y], color=MARK, lw=1.2, zorder=3)
        for xe in (fx(lo), fx(hi)):
            ax.plot([xe, xe], [y - 0.16, y + 0.16], color=MARK, lw=1.2, zorder=3)
        ax.scatter([fx(orr)], [y], marker="s", s=26, color=MARK, zorder=4)
        ax.text(X_CI, y, f"{orr:.2f} ({lo:.2f}-{hi:.2f})", ha="left", va="center",
                fontsize=9.0, color=INK)
        ax.text(X_P, y, ptxt(p), ha="right", va="center", fontsize=9.0, color=INK)
    y -= 1.0

bottom = y + 0.5
ax.plot([0, 100], [bottom, bottom], color=RULEC, lw=1.6)

for orr in TICKS:
    ax.plot([fx(orr), fx(orr)], [bottom, bottom - 0.28], color=INK, lw=0.9)
    ax.text(fx(orr), bottom - 0.85, f"{orr:g}", ha="center", va="top", fontsize=8, color=INK)
ax.annotate("", xy=(fx(0.52), bottom - 1.7), xytext=(fx(0.85), bottom - 1.7),
            arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))
ax.annotate("", xy=(fx(1.9), bottom - 1.7), xytext=(fx(1.18), bottom - 1.7),
            arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))
ax.text(fx(1.0) - 1.4, bottom - 2.3, "Favors lower risk", ha="right", fontsize=7.2, color=INK)
ax.text(fx(1.0) + 1.4, bottom - 2.3, "Favors higher risk", ha="left", fontsize=7.2, color=INK)

fns = [
    "Abbreviations: IA, intracranial aneurysm; aSAH, aneurysmal subarachnoid haemorrhage; OR, odds ratio.",
    "All estimates are random-effects inverse-variance weighted, per year of later age at natural menopause.",
    "$^{a}$Bakker 2020, European. Ruptured and unruptured IA are near-genetically identical (rg approx 0.97),",
    "so these are a consistency check across the disease course, not independent stages.  $^{b}$FinnGen R11,",
    "matched by rsID.  $^{c}$Different outcome, validates the pipeline.  $^{d}$Outcome with no plausible",
    "hormonal link; a null bounds residual pleiotropy.",
]
fy = bottom - 3.0
for i, t in enumerate(fns):
    ax.text(0, fy - i * 0.62, t, ha="left", va="top", fontsize=7.2, color="#333")

ax.text(0, top + 2.0, "Mendelian randomization: later menopause and intracranial aneurysm across the disease course",
        ha="left", fontsize=10.5, fontweight="bold", color=INK)
for ext in ("png", "pdf", "svg"):
    fig.savefig(ROOT / "figures" / f"mr_cascade.{ext}", dpi=200, bbox_inches="tight", facecolor="white")
print("saved mr_cascade.png/pdf/svg")
