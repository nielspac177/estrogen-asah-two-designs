"""Figure 2. Participant flow (STROBE) for the observational cohort, with codes."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

plt.rcParams.update({"font.family": ["Avenir Next", "Arial"], "font.size": 10})
INK = "#1b2733"; BLUE = "#0B6FB8"; MUTE = "#5b6b7a"; PANEL = "#F6FAFD"

fig, ax = plt.subplots(figsize=(9.6, 7.4)); ax.axis("off")
ax.set_xlim(-6, 108); ax.set_ylim(0, 100)


def box(x, y, w, h, title, subs=(), fc=PANEL, ec=BLUE, tc=INK):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                 boxstyle="round,pad=0.3,rounding_size=2", fc=fc, ec=ec, lw=1.5))
    n = 1 + len(subs)
    yy = y + (n - 1) * 2.9 / 2
    ax.text(x, yy, title, ha="center", va="center", fontsize=10.3, fontweight="bold", color=tc)
    yy -= 2.9
    for s in subs:
        ax.text(x, yy, s, ha="center", va="center", fontsize=8.3, color=MUTE)
        yy -= 2.9


def down(x, y1, y2):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>", mutation_scale=15,
                 color=MUTE, lw=1.7))


CX = 38
box(CX, 92, 62, 12, "ICU admissions with a nontraumatic SAH code",
    ["MIMIC-IV + eICU", "ICD-9 430; ICD-10 I60.x"])
down(CX, 86, 80)

# exclusion note (clear of the main column, to the right)
ex_x, ex_y = 90, 74
ax.add_patch(FancyBboxPatch((ex_x - 17, ex_y - 8), 34, 16, boxstyle="round,pad=0.3,rounding_size=2",
             fc="#FBF4E7", ec="#B9812A", lw=1.2))
ax.text(ex_x, ex_y + 5.5, "Excluded", ha="center", va="center", fontsize=9,
        fontweight="bold", color="#6b4a12")
for i, s in enumerate(["Traumatic SAH (800-804; S06.x)", "AVM (747.81; Q28.2)",
                       "No aneurysm evidence"]):
    ax.text(ex_x, ex_y + 2 - i * 3, s, ha="center", va="center", fontsize=7.2, color="#6b4a12")
ax.add_patch(FancyArrowPatch((CX + 31, 74), (ex_x - 17, 74), arrowstyle="-|>", mutation_scale=13,
             color=MUTE, lw=1.3))

box(CX, 74, 62, 13, "Aneurysmal SAH cohort   n = 1,771",
    ["SAH + aneurysm (437.3; I67.1), or", "securing procedure clip/coil (39.51/39.52; O03)"])
# split women / men
box(22, 52, 34, 10, "Women", ["n = 1,105"])
box(60, 52, 32, 10, "Men (reference)", ["n = 666"], ec=MUTE, tc=MUTE)
ax.add_patch(FancyArrowPatch((CX, 67.5), (22, 57.2), arrowstyle="-|>", mutation_scale=14,
             color=MUTE, lw=1.6))
ax.add_patch(FancyArrowPatch((CX, 67.5), (60, 57.2), arrowstyle="-|>", mutation_scale=14,
             color=MUTE, lw=1.6))
down(22, 46.8, 38)
# strata
box(22, 31, 42, 11, "Premenopausal (age < 51)", ["n = 313    DCI 16.3%"])
down(22, 25.5, 20.5)
box(22, 15, 42, 11, "Postmenopausal (age ≥ 51)", ["n = 792    DCI 10.4%"], ec=BLUE)
# primary contrast bracket
ax.text(50, 23, "primary\ncontrast", ha="left", va="center", fontsize=9,
        color=BLUE, fontweight="bold")
ax.annotate("", xy=(43, 31), xytext=(49, 25.5), arrowprops=dict(arrowstyle="-", color=BLUE, lw=1))
ax.annotate("", xy=(43, 15), xytext=(49, 20.5), arrowprops=dict(arrowstyle="-", color=BLUE, lw=1))

fig.savefig("/Volumes/Niels 2/MIMIC/estrogen-asah-two-designs/figures/cohort_flow.png",
            dpi=190, bbox_inches="tight", facecolor="white")
print("saved")
