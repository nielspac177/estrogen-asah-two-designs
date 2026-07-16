"""One-page workflow diagram: triangulating the role of estrogen in aSAH.

Design: swimlanes per arm, solid header banners, soft drop shadows for depth,
symmetric top-fork / bottom-merge right-angle connectors, emphasized result
numbers, restrained palette, clean sans-serif, text centered in every box.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

plt.rcParams.update({"font.family": ["Avenir Next", "Arial"], "font.size": 10})

INK = "#1b2733"; MUTE = "#5b6b7a"
BLUE = "#0B6FB8"; BLUE_BG = "#EAF3FB"; BLUE_PANEL = "#F6FAFD"
GREEN = "#0E8F6E"; GREEN_BG = "#E8F6F1"; GREEN_PANEL = "#F5FCF9"
PUR = "#6D4AAF"; PUR_BG = "#F1ECFA"
AMBER = "#B9812A"; AMBER_BG = "#FBF1DC"
SHADOW = (0, 0, 0, 0.07)

fig, ax = plt.subplots(figsize=(10.6, 8.6)); ax.axis("off")
ax.set_xlim(0, 100); ax.set_ylim(0, 100)


def _round(x, y, w, h, fc, ec, lw, rs):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                 boxstyle=f"round,pad=0.2,rounding_size={rs}", fc=fc, ec=ec, lw=lw))


def card(x, y, w, h, lines, *, head=None, fc="white", ec=BLUE, tc=INK, head_c=None,
         fs=9.3, lh=3.3, rs=2.2, lw=1.4, shadow=True, emph_last=False):
    if shadow:
        _round(x + 0.5, y - 0.6, w, h, SHADOW, "none", 0, rs)
    _round(x, y, w, h, fc, ec, lw, rs)
    items = ([("h", head)] if head else []) + [("b", ln) for ln in lines]
    total = (len(items) - 1) * lh
    top = y + total / 2
    for i, (kind, txt) in enumerate(items):
        yy = top - i * lh
        if kind == "h":
            ax.text(x, yy, txt, ha="center", va="center", fontsize=fs + 1.0,
                    color=head_c or ec, fontweight="bold")
        else:
            last = (i == len(items) - 1)
            ax.text(x, yy, txt, ha="center", va="center",
                    fontsize=fs + (1.3 if emph_last and last else 0),
                    fontweight="bold" if emph_last and last else "normal", color=tc)


def varrow(x, y1, y2, color=MUTE, lw=1.8):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>",
                 mutation_scale=14, color=color, lw=lw, shrinkA=0, shrinkB=0))


# title
ax.text(50, 97.3, "Estrogen and aneurysmal subarachnoid haemorrhage: two study designs",
        ha="center", fontsize=13.5, fontweight="bold", color=INK)

# question
card(50, 90, 92, 7.6, [
    "Does estrogen protect the brain's vessels: fewer aneurysm ruptures (aSAH) and less",
    "delayed cerebral ischaemia?   Protective in animal models · never tested in humans."],
    fc="#FFFBF2", ec=AMBER, fs=9.2, lh=3.1)

# top fork (question -> both arm headers), mirrors the bottom merge
fork_y = 83.5
ax.add_line(Line2D([50, 50], [86.2, fork_y], color=MUTE, lw=1.8))
ax.add_line(Line2D([26.5, 73.5], [fork_y, fork_y], color=MUTE, lw=1.8))
varrow(26.5, fork_y, 80.3, color=BLUE)
varrow(73.5, fork_y, 80.3, color=GREEN)

# swimlane panels
_round(26.5, 54.9, 43, 50.5, BLUE_PANEL, BLUE, 1.1, 3)
_round(73.5, 54.9, 43, 50.5, GREEN_PANEL, GREEN, 1.1, 3)

# solid header banners
card(26.5, 76, 39, 6.2, [], head="ARM 1 · Observational (ICU)", fc=BLUE, ec=BLUE,
     head_c="white", fs=9.6, shadow=False)
card(73.5, 76, 39, 6.2, [], head="ARM 2 · Genetic (MR)", fc=GREEN, ec=GREEN,
     head_c="white", fs=9.6, shadow=False)

# setup
card(26.5, 64.5, 38, 11, [
    "MIMIC-IV + eICU", "1,771 aSAH patients", "Menopausal state → DCI / death"],
    ec=BLUE, fs=9.1, lh=3.2)
card(73.5, 64.5, 38, 11, [
    "Public GWAS (100,000s of people)", "Menopause / SHBG gene variants",
    "→ aneurysm GWAS (Bakker 2020)"], ec=GREEN, fs=9.1, lh=3.2)

# weakness
card(26.5, 51, 38, 9, ["confounded by age", "(menopause ≈ age; no hormone data)"],
     head="WEAKNESS", fc=AMBER_BG, ec=AMBER, tc="#6b4a12", fs=8.7, lh=3.0)
card(73.5, 51, 38, 9, ["genetic pleiotropy", "(but immune to age confounding)"],
     head="WEAKNESS", fc=AMBER_BG, ec=AMBER, tc="#6b4a12", fs=8.7, lh=3.0)

# result (emphasize the OR)
card(26.5, 38, 38, 9, ["No protective effect", "DiD 1.04 (0.62–1.76)"],
     head="RESULT", fc="white", ec=BLUE, fs=9.2, lh=3.1, emph_last=True)
card(73.5, 38, 38, 9, ["No protective effect", "OR 1.03 (0.97–1.09)"],
     head="RESULT", fc="white", ec=GREEN, fs=9.2, lh=3.1, emph_last=True)

# intra-lane arrows
for x, c in [(26.5, BLUE), (73.5, GREEN)]:
    varrow(x, 72.8, 70.1, color=c)
    varrow(x, 59.0, 55.6, color=c)
    varrow(x, 46.5, 42.6, color=c)

# bottom merge (right-angle)
rail = 24
for x, c in [(26.5, BLUE), (73.5, GREEN)]:
    ax.add_line(Line2D([x, x], [33.5, rail], color=c, lw=1.8))
ax.add_line(Line2D([26.5, 73.5], [rail, rail], color=MUTE, lw=1.8))
varrow(50, rail, 21.4, color=PUR)

# convergence
card(50, 15.8, 86, 9, [
    "Two designs with different weaknesses point the same way.",
    "Age confounding would bias the ICU arm, but the genetic arm, which age cannot touch, agrees."],
    head="CONVERGENT EVIDENCE", fc=PUR_BG, ec=PUR, fs=8.9, lh=3.0)
varrow(50, 11.3, 7.9, color=PUR)

# conclusion
card(50, 4.3, 92, 5.6, [
    "Estrogen does not measurably protect against aneurysmal SAH in humans."],
    fc="#12303a", ec="#12303a", tc="white", fs=10.2, rs=1.8, shadow=True)

fig.savefig("/Volumes/Niels 2/MIMIC/estrogen-asah-two-designs/figures/study_workflow.png",
            dpi=180, bbox_inches="tight", facecolor="white")
print("saved")
