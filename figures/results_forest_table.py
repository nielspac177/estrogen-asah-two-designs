"""Figure 3. Results table with an embedded forest plot (JAMA/NEJM style).

A single figure laid out as a table: text columns on the left and right, a central
forest column on a shared log-odds-ratio axis with a reference line at OR = 1.
Rows are grouped by design. Values are the final estimates from both arms.
"""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

plt.rcParams.update({"font.family": ["Avenir Next", "Arial"], "font.size": 9.5})
INK = "#1b2733"; MUTE = "#5b6b7a"; RED = "#c0392b"; GRID = "#e7ebef"
BLUE = "#0B6FB8"; GREEN = "#0E8F6E"; ORANGE = "#C77A0A"

# section, label, N text, OR, lo, hi, P text, colour
ROWS = [
    ("H", "Observational (ICU cohort)", "", None, None, None, "", BLUE),
    ("R", "Postmenopausal vs premenopausal", "1105", 0.86, 0.58, 1.28, ".46", BLUE),
    ("R", "Age × sex difference-in-differences", "1771", 1.04, 0.62, 1.76, "", BLUE),
    ("H", "Genetic (Mendelian randomization), aSAH", "", None, None, None, "", GREEN),
    ("R", "Age at menopause, IVW", "85", 1.03, 0.97, 1.09, ".32", GREEN),
    ("R", "Age at menopause, MR-Egger", "85", 1.10, 0.98, 1.24, ".10", GREEN),
    ("R", "Age at menopause, weighted median", "85", 1.09, 1.00, 1.19, ".06", GREEN),
    ("R", "SHBG, women", "82", 0.73, 0.41, 1.31, ".29", GREEN),
    ("R", "Total testosterone", "58", 0.98, 0.77, 1.27, ".91", GREEN),
    ("R", "SHBG, women (MVMR)", "106", 1.06, 0.50, 2.28, ".88", GREEN),
    ("R", "Bioavailable testosterone (MVMR)", "106", 1.03, 0.57, 1.86, ".93", GREEN),
    ("H", "Positive control", "", None, None, None, "", ORANGE),
    ("R", "Age at menopause → breast cancer", "207", 1.055, 1.041, 1.069, "<.001", ORANGE),
]

# column x-positions (0-100) and forest axis
X_LBL, X_N = 1.5, 44
FX_L, FX_R = 47, 74          # forest column pixels
OR_LO, OR_HI = 0.4, 2.6
X_CI, X_P = 76.5, 99
TICKS = [0.5, 0.75, 1.0, 1.5, 2.0]

def fx(orr):
    return FX_L + (np.log(orr) - np.log(OR_LO)) / (np.log(OR_HI) - np.log(OR_LO)) * (FX_R - FX_L)

n = len(ROWS)
fig, ax = plt.subplots(figsize=(9.6, 0.5 * n + 1.6))
ax.axis("off"); ax.set_xlim(0, 100)
top = n + 1.2
ax.set_ylim(-1.4, top + 1.2)

# column headers
hy = top
for x, t, ha in [(X_LBL, "Analysis", "left"), (X_N, "No.", "right"),
                 (X_CI, "OR (95% CI)", "left"), (X_P, "P", "right")]:
    ax.text(x, hy, t, ha=ha, va="center", fontsize=9.5, fontweight="bold", color=INK)
ax.text((FX_L + FX_R) / 2, hy, "Odds ratio (95% CI)", ha="center", va="center",
        fontsize=9.5, fontweight="bold", color=INK)
ax.plot([0.5, 99.5], [hy - 0.55, hy - 0.55], color=INK, lw=1.1)

# forest gridlines + reference line spanning data rows
y_bottom = 0.4
for orr in TICKS:
    ax.plot([fx(orr), fx(orr)], [y_bottom, hy - 0.7], color=GRID, lw=0.9, zorder=0)
ax.plot([fx(1.0), fx(1.0)], [y_bottom, hy - 0.7], color=RED, ls="--", lw=1.2, zorder=1)

# rows
y = top - 1.1
for kind, label, ntxt, orr, lo, hi, ptxt, c in ROWS:
    if kind == "H":
        ax.text(X_LBL, y, label, ha="left", va="center", fontsize=9.4,
                fontweight="bold", color=c)
    else:
        ax.text(X_LBL + 1.5, y, label, ha="left", va="center", fontsize=9.2, color=INK)
        ax.text(X_N, y, ntxt, ha="right", va="center", fontsize=9.0, color=MUTE)
        ax.plot([fx(lo), fx(hi)], [y, y], color=c, lw=2.0, solid_capstyle="round", zorder=3)
        ax.add_patch(Circle((fx(orr), y), 0.16, facecolor=c, edgecolor="white", lw=0.6,
                            zorder=4, transform=ax.transData))
        ax.text(X_CI, y, f"{orr:.2f} ({lo:.2f}–{hi:.2f})", ha="left", va="center",
                fontsize=9.0, color=INK)
        ax.text(X_P, y, ptxt, ha="right", va="center", fontsize=9.0, color=INK)
    y -= 1.0

# forest axis ticks
for orr in TICKS:
    ax.text(fx(orr), y_bottom - 0.55, f"{orr:g}", ha="center", va="top", fontsize=8, color=MUTE)
ax.annotate("", xy=(fx(0.62), y_bottom - 1.15), xytext=(fx(0.9), y_bottom - 1.15),
            arrowprops=dict(arrowstyle="->", color=MUTE, lw=1))
ax.annotate("", xy=(fx(1.6), y_bottom - 1.15), xytext=(fx(1.11), y_bottom - 1.15),
            arrowprops=dict(arrowstyle="->", color=MUTE, lw=1))
ax.text(fx(0.62), y_bottom - 1.7, "protective", ha="center", fontsize=7.6, color=MUTE)
ax.text(fx(1.6), y_bottom - 1.7, "harmful", ha="center", fontsize=7.6, color=MUTE)

ax.set_title("Estrogen and aneurysmal SAH: estimates across both designs",
             fontsize=11.5, fontweight="bold", pad=10, loc="left", x=0.01)
fig.savefig("/Volumes/Niels 2/MIMIC/estrogen-asah-two-designs/figures/results_forest_table.png",
            dpi=190, bbox_inches="tight", facecolor="white")
print("saved")
