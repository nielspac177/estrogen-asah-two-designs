"""Figure 3. Results table with an embedded forest plot, JAMA table style.

Horizontal rules only (no vertical gridlines or box), monochrome markers, square
points with capped confidence intervals, grouped section headers with light
shading, and footnotes. Values are the final estimates from both arms.
"""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({"font.family": ["Helvetica Neue", "Arial"], "font.size": 9.5})
INK = "#141414"; MARK = "#222222"; RULEC = "#141414"; GREY = "#8a8a8a"; SHADE = "#f0f0f0"

# kind, label, No., OR, lo, hi, P, footnote-mark
ROWS = [
    ("H", "Observational (ICU cohort)", "", None, None, None, "", ""),
    ("R", "Postmenopausal vs premenopausal", "1105", 0.86, 0.58, 1.28, ".46", "a"),
    ("R", "Age × sex difference-in-differences", "1771", 1.04, 0.62, 1.76, "", "a"),
    ("H", "Genetic (Mendelian randomization), aSAH", "", None, None, None, "", ""),
    ("R", "Age at menopause, IVW", "85", 1.03, 0.97, 1.09, ".32", "b"),
    ("R", "Age at menopause, MR-Egger", "85", 1.10, 0.98, 1.24, ".10", "b"),
    ("R", "Age at menopause, weighted median", "85", 1.09, 1.00, 1.19, ".06", "b"),
    ("R", "SHBG, women", "82", 0.73, 0.41, 1.31, ".29", "b"),
    ("R", "Total testosterone", "58", 0.98, 0.77, 1.27, ".91", "b"),
    ("R", "SHBG, women, MVMR", "106", 1.06, 0.50, 2.28, ".88", "c"),
    ("R", "Bioavailable testosterone, MVMR", "106", 1.03, 0.57, 1.86, ".93", "c"),
    ("H", "Positive control", "", None, None, None, "", ""),
    ("R", "Age at menopause and breast cancer", "207", 1.055, 1.041, 1.069, "<.001", "d"),
]

X_LBL, X_N = 1.5, 43
FX_L, FX_R = 46, 70
OR_LO, OR_HI = 0.4, 2.6
X_CI, X_P = 72.5, 99
TICKS = [0.5, 0.75, 1.0, 1.5, 2.0]

def fx(orr):
    return FX_L + (np.log(orr) - np.log(OR_LO)) / (np.log(OR_HI) - np.log(OR_LO)) * (FX_R - FX_L)

n = len(ROWS)
fig, ax = plt.subplots(figsize=(9.4, 0.52 * n + 2.2))
ax.axis("off"); ax.set_xlim(0, 100)
top = n + 1.0
ax.set_ylim(-3.6, top + 1.4)

# top rule + column headers + header rule
ax.plot([0, 100], [top + 0.9, top + 0.9], color=RULEC, lw=1.6)
for x, t, ha in [(X_LBL, "Analysis", "left"), (X_N, "No.", "right"),
                 (X_CI, "OR (95% CI)", "left"), (X_P, "P Value", "right")]:
    ax.text(x, top + 0.1, t, ha=ha, va="center", fontsize=9.5, fontweight="bold", color=INK)
ax.plot([0, 100], [top - 0.55, top - 0.55], color=RULEC, lw=1.0)

y0 = top - 1.5
# reference line at OR=1 spanning the data rows
ax.plot([fx(1.0), fx(1.0)], [-0.4, top - 0.7], color=GREY, lw=0.9, zorder=0)

y = y0
for kind, label, ntxt, orr, lo, hi, ptxt, fn in ROWS:
    if kind == "H":
        ax.add_patch(Rectangle((0, y - 0.5), 100, 1.0, facecolor=SHADE, edgecolor="none", zorder=0))
        ax.text(X_LBL, y, label, ha="left", va="center", fontsize=9.3, fontweight="bold", color=INK)
    else:
        lab = label + (f"$^{{{fn}}}$" if fn else "")
        ax.text(X_LBL + 1.2, y, lab, ha="left", va="center", fontsize=9.2, color=INK)
        ax.text(X_N, y, ntxt, ha="right", va="center", fontsize=9.0, color=INK)
        # CI line with end caps + square marker
        ax.plot([fx(lo), fx(hi)], [y, y], color=MARK, lw=1.2, zorder=3)
        for xe in (fx(lo), fx(hi)):
            ax.plot([xe, xe], [y - 0.16, y + 0.16], color=MARK, lw=1.2, zorder=3)
        ax.scatter([fx(orr)], [y], marker="s", s=26, color=MARK, zorder=4)
        ax.text(X_CI, y, f"{orr:.2f} ({lo:.2f}-{hi:.2f})", ha="left", va="center",
                fontsize=9.0, color=INK)
        ax.text(X_P, y, ptxt, ha="right", va="center", fontsize=9.0, color=INK)
    y -= 1.0

bottom = y + 0.5
ax.plot([0, 100], [bottom, bottom], color=RULEC, lw=1.6)

# forest axis
for orr in TICKS:
    ax.plot([fx(orr), fx(orr)], [bottom, bottom - 0.28], color=INK, lw=0.9)
    ax.text(fx(orr), bottom - 0.85, f"{orr:g}", ha="center", va="top", fontsize=8, color=INK)
ax.annotate("", xy=(fx(0.52), bottom - 1.7), xytext=(fx(0.85), bottom - 1.7),
            arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))
ax.annotate("", xy=(fx(1.9), bottom - 1.7), xytext=(fx(1.18), bottom - 1.7),
            arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))
ax.text(fx(1.0) - 1.4, bottom - 2.3, "Favors lower risk", ha="right", fontsize=7.2, color=INK)
ax.text(fx(1.0) + 1.4, bottom - 2.3, "Favors higher risk", ha="left", fontsize=7.2, color=INK)

# footnotes
fns = [
    "Abbreviations: aSAH, aneurysmal subarachnoid haemorrhage; IVW, inverse-variance weighted;",
    "MVMR, multivariable Mendelian randomization; OR, odds ratio; SHBG, sex hormone-binding globulin.",
    "$^{a}$Adjusted for hypertension, smoking, and diabetes, with cluster-robust SEs by hospital; age is",
    "excluded because it defines the exposure.  $^{b}$Two-sample MR, random-effects IVW unless noted.",
    "$^{c}$Multivariable MR, mutually adjusted for the other hormone.  $^{d}$Different outcome (breast",
    "cancer), shown to validate the genetic pipeline.",
]
fy = bottom - 3.0
for i, t in enumerate(fns):
    ax.text(0, fy - i * 0.62, t, ha="left", va="top", fontsize=7.2, color="#333")

ax.text(0, top + 2.0, "Table. Estimates for estrogen and aneurysmal SAH across both study designs",
        ha="left", fontsize=11, fontweight="bold", color=INK)
fig.savefig("/Volumes/Niels 2/MIMIC/estrogen-asah-two-designs/figures/results_forest_table.png",
            dpi=200, bbox_inches="tight", facecolor="white")
print("saved")
