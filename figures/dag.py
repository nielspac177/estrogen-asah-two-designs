"""Figure 1B. Causal diagrams (DAGs) for the two arms.

Left: the observational arm, where age both defines the exposure and drives the
outcome (confounding), and ICU admission is a selection node. Right: the genetic
arm, where variants assigned at conception instrument menopausal timing and bypass
age; the dashed path is the pleiotropy assumption that the sensitivity analyses test.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch

plt.rcParams.update({"font.family": ["Avenir Next", "Arial"], "font.size": 10})
INK = "#1b2733"; BLUE = "#0B6FB8"; GREEN = "#0E8F6E"; RED = "#c0392b"; MUTE = "#5b6b7a"

fig, ax = plt.subplots(figsize=(11, 4.8)); ax.axis("off")
ax.set_xlim(0, 100); ax.set_ylim(0, 50)


def node(x, y, text, ec=INK, fc="white", w=15, h=7, fs=9.5, tc=INK):
    ax.add_patch(Ellipse((x, y), w, h, facecolor=fc, edgecolor=ec, lw=1.6, zorder=2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=tc, zorder=3)


def arrow(p1, p2, color=INK, style="-|>", ls="-", lw=1.7, rad=0.0):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle=style, mutation_scale=15,
                 color=color, lw=lw, linestyle=ls, shrinkA=10, shrinkB=10,
                 connectionstyle=f"arc3,rad={rad}", zorder=1))


# ---------------- Panel A: observational ----------------
ax.text(24, 47, "A. Observational (ICU cohort)", ha="center", fontsize=11.5,
        fontweight="bold", color=BLUE)
node(12, 36, "Age", ec=RED)
node(12, 20, "Menopausal\nstate", ec=BLUE)
node(38, 20, "DCI / aSAH", ec=INK)
node(25, 6, "ICU admission\n(selection)", ec=MUTE, fc="#f2f4f6", w=18, tc=MUTE, fs=8.6)

arrow((12, 32.5), (12, 23.5), color=RED)                       # Age -> Menopause (defines)
ax.text(6.5, 28, "defines", fontsize=8, color=RED, rotation=90, va="center")
arrow((16, 33), (34, 22.5), color=RED)                         # Age -> DCI (confounding)
ax.text(27, 31, "vasospasm", fontsize=8, color=RED)
arrow((19.5, 20), (30.5, 20), color=BLUE)                      # Menopause -> DCI (effect?)
ax.text(24.5, 22, "?", fontsize=13, color=BLUE, ha="center", fontweight="bold")
arrow((13, 16.5), (21, 9), color=MUTE, lw=1.2)
arrow((37, 16.5), (30, 9), color=MUTE, lw=1.2)

# ---------------- Panel B: genetic (MR) ----------------
ax.text(76, 47, "B. Genetic (Mendelian randomization)", ha="center", fontsize=11.5,
        fontweight="bold", color=GREEN)
node(62, 36, "Genetic\nvariants", ec=GREEN)
node(62, 20, "Menopausal\ntiming", ec=BLUE)
node(90, 20, "aSAH", ec=INK)

arrow((62, 32.5), (62, 23.5), color=GREEN)                     # G -> Menopause (instrument)
ax.text(55, 28, "instrument", fontsize=8, color=GREEN, rotation=90, va="center")
arrow((69.5, 20), (82.5, 20), color=BLUE)                      # Menopause -> aSAH (causal?)
ax.text(76, 22, "?", fontsize=13, color=BLUE, ha="center", fontweight="bold")
arrow((66, 33), (86, 23.5), color=MUTE, ls=(0, (4, 3)), lw=1.3, rad=-0.25)   # pleiotropy
ax.text(80, 33.5, "pleiotropy\n(assumed 0; tested)", fontsize=8, color=MUTE, ha="center")
ax.text(62, 43, "assigned at conception,\nindependent of age", fontsize=8,
        color=MUTE, ha="center", style="italic")

ax.plot([50, 50], [3, 45], color="#e0e5ea", lw=1)              # divider

fig.savefig("/Volumes/Niels 2/MIMIC/estrogen-asah-two-designs/figures/dag.png",
            dpi=190, bbox_inches="tight", facecolor="white")
print("saved")
