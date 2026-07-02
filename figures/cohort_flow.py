"""Figure 2. Participant flow (STROBE) for the observational cohort."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

plt.rcParams.update({"font.family": ["Avenir Next", "Arial"], "font.size": 10})
INK = "#1b2733"; BLUE = "#0B6FB8"; MUTE = "#5b6b7a"; PANEL = "#F6FAFD"

fig, ax = plt.subplots(figsize=(8.4, 7.2)); ax.axis("off")
ax.set_xlim(0, 100); ax.set_ylim(0, 100)


def box(x, y, w, h, lines, fc=PANEL, ec=BLUE, fs=10, tc=INK, bold_first=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                 boxstyle="round,pad=0.3,rounding_size=2", fc=fc, ec=ec, lw=1.5))
    total = (len(lines) - 1) * 3.4
    yy = y + total / 2
    for i, ln in enumerate(lines):
        ax.text(x, yy, ln, ha="center", va="center", fontsize=fs + (0.6 if bold_first and i == 0 else 0),
                fontweight="bold" if bold_first and i == 0 else "normal", color=tc)
        yy -= 3.4


def down(x, y1, y2):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>", mutation_scale=15,
                 color=MUTE, lw=1.7))


CX = 38
box(CX, 92, 58, 10, ["ICU admissions, MIMIC-IV + eICU",
                     "with a nontraumatic subarachnoid-haemorrhage code"], bold_first=True)
down(CX, 87, 79.5)
# exclusion note
ax.add_patch(FancyBboxPatch((72, 66), 26, 12, boxstyle="round,pad=0.3,rounding_size=2",
             fc="#FBF4E7", ec="#B9812A", lw=1.2))
ax.text(85, 72, "Excluded\ntraumatic SAH\narteriovenous malformation\nno aneurysm evidence",
        ha="center", va="center", fontsize=8.4, color="#6b4a12")
ax.add_patch(FancyArrowPatch((CX + 29, 72), (72, 72), arrowstyle="-|>", mutation_scale=13,
             color=MUTE, lw=1.3))

box(CX, 72, 58, 11, ["Aneurysmal SAH cohort", "n = 1,771",
                     "(SAH + aneurysm diagnosis or securing procedure)"], bold_first=True)
# split women / men
box(20, 50, 34, 10, ["Women", "n = 1,105"], bold_first=True)
box(60, 50, 30, 10, ["Men (reference)", "n = 666"], ec=MUTE, tc=MUTE, bold_first=True)
ax.add_patch(FancyArrowPatch((CX, 66.5), (20, 55.2), arrowstyle="-|>", mutation_scale=14,
             color=MUTE, lw=1.6))
ax.add_patch(FancyArrowPatch((CX, 66.5), (60, 55.2), arrowstyle="-|>", mutation_scale=14,
             color=MUTE, lw=1.6))
down(20, 44.8, 36)
# strata
box(20, 30, 40, 11, ["Premenopausal (age < 51)", "n = 313", "DCI 16.3%"], bold_first=True)
down(20, 44.8, 41.5)
box(20, 14, 40, 11, ["Postmenopausal (age ≥ 51)", "n = 792", "DCI 10.4%"], ec=BLUE, bold_first=True)
ax.add_patch(FancyArrowPatch((20, 24.5), (20, 19.5), arrowstyle="-|>", mutation_scale=14,
             color=MUTE, lw=1.6))
# primary contrast bracket
ax.text(46, 22, "primary\ncontrast", ha="center", va="center", fontsize=9,
        color=BLUE, fontweight="bold")
ax.annotate("", xy=(40, 30), xytext=(46, 25), arrowprops=dict(arrowstyle="-", color=BLUE, lw=1))
ax.annotate("", xy=(40, 14), xytext=(46, 19), arrowprops=dict(arrowstyle="-", color=BLUE, lw=1))

fig.savefig("/Volumes/Niels 2/MIMIC/estrogen-asah-two-designs/figures/cohort_flow.png",
            dpi=190, bbox_inches="tight", facecolor="white")
print("saved")
