"""Graphical abstract (journal TOC): question -> two designs -> convergent null."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

plt.rcParams.update({"font.family": ["Avenir Next", "Arial"], "font.size": 11})
INK = "#1b2733"; MUTE = "#5b6b7a"
BLUE = "#0B6FB8"; GREEN = "#0E8F6E"; ORANGE = "#E69F00"; RED = "#c0392b"; PUR = "#6D4AAF"

fig, ax = plt.subplots(figsize=(12, 5.4)); ax.axis("off")
ax.set_xlim(0, 122); ax.set_ylim(0, 54)

ax.text(61, 51, "Estrogen and aneurysmal subarachnoid haemorrhage",
        ha="center", fontsize=15, fontweight="bold", color=INK)

# ---------- Zone 1: the question ----------
ax.text(18, 44, "The question", ha="center", fontsize=11, fontweight="bold", color=ORANGE)
ax.plot([6, 30], [37, 37], color="#b5651d", lw=7, solid_capstyle="round", zorder=1)
ax.add_patch(Circle((18, 40), 3.2, facecolor="#e8a06b", edgecolor="#b5651d", lw=2, zorder=2))
ax.text(18, 31, "Does estrogen protect against\naneurysm rupture (aSAH)?",
        ha="center", fontsize=11, color=INK)
ax.text(18, 24.5, "Animal models: yes\nHumans: untested",
        ha="center", fontsize=9.5, color=MUTE, style="italic")

# ---------- Zone 2: two designs, mini forest ----------
ax.text(60, 44, "Two designs, non-overlapping biases", ha="center", fontsize=11,
        fontweight="bold", color=PUR)

x0, x1 = 42, 74
def fx(orr):
    lo, hi = np.log(0.5), np.log(2.0)
    return x0 + (np.log(orr) - lo) / (hi - lo) * (x1 - x0)

for orr in (0.5, 0.75, 1.0, 1.5, 2.0):
    ax.plot([fx(orr), fx(orr)], [13, 34], color="#eef1f4", lw=1, zorder=0)
    ax.text(fx(orr), 10.8, f"{orr:g}", ha="center", fontsize=8, color=MUTE)
ax.plot([fx(1.0), fx(1.0)], [13, 36], color=RED, ls="--", lw=1.4, zorder=1)
ax.text((x0 + x1) / 2, 8, "odds ratio (log scale)  ·  1.0 = no effect",
        ha="center", fontsize=8.5, color=MUTE)

rows = [
    ("Observational, age-by-sex DiD  (n=1,771)",   1.04, 0.62, 1.76, BLUE, 31),
    ("Genetic  (Mendelian randomization)",         1.03, 0.97, 1.09, GREEN, 22),
]
for label, orr, lo, hi, c, yy in rows:
    ax.plot([fx(lo), fx(hi)], [yy, yy], color=c, lw=3, solid_capstyle="round", zorder=3)
    ax.add_patch(Circle((fx(orr), yy), 0.65, facecolor=c, edgecolor="white", lw=1, zorder=4))
    ax.text(x0, yy + 2.7, label, ha="left", fontsize=9.3, color=INK)
    ax.text(x1 + 1.5, yy, f"OR {orr:.2f}", ha="left", va="center", fontsize=9.4,
            fontweight="bold", color=c)

ax.text(60, 15.5, "Positive control (menopause to breast cancer): recovered",
        ha="center", fontsize=8.8, color=ORANGE)

# ---------- Zone 3: conclusion ----------
ax.add_patch(FancyBboxPatch((97, 19), 22, 17, boxstyle="round,pad=0.4,rounding_size=2",
             facecolor="#12303a", edgecolor="#12303a"))
ax.text(108, 31, "No protective\neffect", ha="center", va="center", fontsize=13.5,
        fontweight="bold", color="white")
ax.text(108, 23.5, "two methods agree", ha="center", va="center", fontsize=9, color="#cfe0e6")

# arrows between zones (clear of labels)
ax.add_patch(FancyArrowPatch((32, 27), (40, 27), arrowstyle="-|>", mutation_scale=17,
             color=MUTE, lw=2))
ax.add_patch(FancyArrowPatch((88, 27), (95.5, 27), arrowstyle="-|>", mutation_scale=17,
             color=MUTE, lw=2))

fig.savefig("/Volumes/Niels 2/MIMIC/estrogen-asah-two-designs/figures/graphical_abstract.png",
            dpi=200, bbox_inches="tight", facecolor="white")
print("saved")
