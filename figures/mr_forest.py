"""Genetic-arm forest plot: sex-hormone / reproductive exposures on aSAH, plus the
positive control. Values are the final MR results (see estrogen-aneurysm-mr)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch  # noqa: F401

plt.rcParams.update({"font.family": ["Avenir Next", "Arial"], "font.size": 10})
BLUE = "#0B6FB8"; GREEN = "#0E8F6E"; ORANGE = "#E69F00"; GREY = "#8a97a3"; INK = "#1b2733"

# (label, OR, lo, hi, colour)  top-to-bottom
rows = [
    ("Age at natural menopause  (per year)",        1.03, 0.97, 1.09, BLUE),
    ("SHBG, women  (per SD, univariable)",           0.73, 0.41, 1.31, BLUE),
    ("Total testosterone  (per SD)",                 0.98, 0.77, 1.27, BLUE),
    ("SHBG, women  (per SD, MVMR adj. testosterone)", 1.06, 0.50, 2.28, GREEN),
    ("Bioavailable testosterone  (per SD, MVMR)",    1.03, 0.57, 1.86, GREEN),
    ("__sep__", 0, 0, 0, GREY),
    ("Age at menopause -> breast cancer  (positive control)", 1.055, 1.041, 1.069, ORANGE),
]

fig, ax = plt.subplots(figsize=(8.6, 4.4))
ylabels, yticks = [], []
y = len(rows)
for label, orr, lo, hi, c in rows:
    if label == "__sep__":
        ax.axhline(y, color=GREY, lw=0.8, ls=":")
        y -= 1
        continue
    ax.plot([lo, hi], [y, y], color=c, lw=2.2, solid_capstyle="round", zorder=2)
    ax.scatter([orr], [y], color=c, s=46, zorder=3)
    ax.text(3.15, y, f"{orr:.2f} ({lo:.2f}-{hi:.2f})", va="center", ha="left", fontsize=8.6, color=INK)
    ylabels.append(label)
    yticks.append(y)
    y -= 1

ax.axvline(1.0, color="#c0392b", ls="--", lw=1)
ax.set_xscale("log")
ax.minorticks_off()
ax.set_xticks([0.5, 0.75, 1.0, 1.5, 2.0])
ax.set_xticklabels(["0.5", "0.75", "1.0", "1.5", "2.0"])
ax.set_xlim(0.36, 2.5)
ax.set_yticks(yticks)
ax.set_yticklabels(ylabels, fontsize=8.8)
ax.set_xlabel("Odds ratio (log scale)")
ax.set_title("Mendelian randomization: sex-hormone exposures and aneurysmal SAH", fontsize=11.5,
             fontweight="bold", pad=12)
ax.text(3.15, len(rows) + 0.4, "OR (95% CI)", ha="left", fontsize=8.6, fontweight="bold", color=INK)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.tick_params(left=False)
ax.margins(y=0.08)
fig.savefig("/Volumes/Niels 2/MIMIC/estrogen-asah-two-designs/figures/mr_forest.png",
            dpi=180, bbox_inches="tight", facecolor="white")
print("saved")
