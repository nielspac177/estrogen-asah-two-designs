"""Methods-only workflow: the two-arm triangulation DESIGN (no results).

Same visual language as triangulation_workflow.py, but shows the analytical
pipeline of each arm and ends at the triangulation step — no effect estimates,
no conclusion.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

plt.rcParams.update({"font.family": ["Avenir Next", "Arial"], "font.size": 10})

INK = "#1b2733"; MUTE = "#5b6b7a"
BLUE = "#0B6FB8"; BLUE_PANEL = "#F6FAFD"
GREEN = "#0E8F6E"; GREEN_PANEL = "#F5FCF9"
PUR = "#6D4AAF"; PUR_BG = "#F1ECFA"
AMBER = "#B9812A"; AMBER_BG = "#FFFBF2"
SHADOW = (0, 0, 0, 0.07)

fig, ax = plt.subplots(figsize=(10.6, 8.8)); ax.axis("off")
ax.set_xlim(0, 100); ax.set_ylim(0, 100)


def _round(x, y, w, h, fc, ec, lw, rs):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                 boxstyle=f"round,pad=0.2,rounding_size={rs}", fc=fc, ec=ec, lw=lw))


def card(x, y, w, h, lines, *, head=None, fc="white", ec=BLUE, tc=INK, head_c=None,
         fs=9.1, lh=3.1, rs=2.2, lw=1.4, shadow=True):
    if shadow:
        _round(x + 0.5, y - 0.6, w, h, SHADOW, "none", 0, rs)
    _round(x, y, w, h, fc, ec, lw, rs)
    items = ([("h", head)] if head else []) + [("b", ln) for ln in lines]
    total = (len(items) - 1) * lh
    top = y + total / 2
    for i, (kind, txt) in enumerate(items):
        yy = top - i * lh
        if kind == "h":
            ax.text(x, yy, txt, ha="center", va="center", fontsize=fs + 0.4,
                    color=head_c or ec, fontweight="bold")
        else:
            ax.text(x, yy, txt, ha="center", va="center", fontsize=fs, color=tc)


def varrow(x, y1, y2, color=MUTE, lw=1.8):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>",
                 mutation_scale=14, color=color, lw=lw, shrinkA=0, shrinkB=0))


# title
ax.text(50, 97.3, "Study design: estrogen and aneurysmal SAH — two complementary approaches",
        ha="center", fontsize=13.5, fontweight="bold", color=INK)

# aim
card(50, 90.5, 92, 7, [
    "Aim: does estrogen influence the risk of aneurysmal SAH / delayed cerebral ischaemia?",
    "Estimate it two ways whose weaknesses do not overlap, then read the answers together."],
    fc=AMBER_BG, ec=AMBER, fs=9.2, lh=3.1)

# top fork
fork_y = 84
ax.add_line(Line2D([50, 50], [86.8, fork_y], color=MUTE, lw=1.8))
ax.add_line(Line2D([26.5, 73.5], [fork_y, fork_y], color=MUTE, lw=1.8))
varrow(26.5, fork_y, 80.8, color=BLUE)
varrow(73.5, fork_y, 80.8, color=GREEN)

# panels
_round(26.5, 52.5, 43, 55, BLUE_PANEL, BLUE, 1.1, 3)
_round(73.5, 52.5, 43, 55, GREEN_PANEL, GREEN, 1.1, 3)

# headers
card(26.5, 76.5, 39, 6, [], head="ARM 1 · Observational (ICU)", fc=BLUE, ec=BLUE,
     head_c="white", fs=9.6, shadow=False)
card(73.5, 76.5, 39, 6, [], head="ARM 2 · Genetic (Mendelian randomization)", fc=GREEN,
     ec=GREEN, head_c="white", fs=8.9, shadow=False)

# --- Arm 1 steps ---
card(26.5, 66, 38, 10.5, [
    "MIMIC-IV + eICU", "aneurysmal SAH cohort (n = 1,771)"],
    head="DATA", fc="white", ec=BLUE, fs=9.0, lh=3.0)
card(26.5, 52.5, 38, 12.5, [
    "Exposure: menopausal state (age proxy)",
    "Outcome: DCI / vasospasm, mortality",
    "Covariates: hypertension, smoking, diabetes"],
    head="DESIGN", fc="white", ec=BLUE, fs=8.8, lh=3.0)
card(26.5, 37.5, 38, 12.5, [
    "Multilevel logistic + IPTW",
    "Specification curve (all forks)",
    "Age × sex difference-in-differences"],
    head="ANALYSIS", fc="white", ec=BLUE, fs=8.8, lh=3.0)

# --- Arm 2 steps ---
card(73.5, 66, 38, 10.5, [
    "Public GWAS summary statistics",
    "Exposures: menopause, SHBG, testosterone"],
    head="DATA", fc="white", ec=GREEN, fs=9.0, lh=3.0)
card(73.5, 52.5, 38, 12.5, [
    "Instruments: genome-wide SNPs (p<5e-8)",
    "Clump + allele harmonization",
    "Outcome: aneurysm GWAS (Bakker 2020)"],
    head="INSTRUMENTS", fc="white", ec=GREEN, fs=8.8, lh=3.0)
card(73.5, 37.5, 38, 12.5, [
    "Two-sample MR: IVW, MR-Egger,",
    "weighted median",
    "Sensitivity: Steiger, MR-PRESSO, leave-one-out"],
    head="ANALYSIS", fc="white", ec=GREEN, fs=8.5, lh=3.0)

# intra-lane arrows
for x, c in [(26.5, BLUE), (73.5, GREEN)]:
    varrow(x, 73.3, 71.4, color=c)   # header -> data
    varrow(x, 60.6, 58.9, color=c)   # data -> design
    varrow(x, 46.2, 43.9, color=c)   # design -> analysis

# merge
rail = 22.5
for x, c in [(26.5, BLUE), (73.5, GREEN)]:
    ax.add_line(Line2D([x, x], [31.2, rail], color=c, lw=1.8))
ax.add_line(Line2D([26.5, 73.5], [rail, rail], color=MUTE, lw=1.8))
varrow(50, rail, 17.3, color=PUR)

# synthesis (design statement, not a result)
card(50, 10.5, 86, 10.5, [
    "Read the estrogen–aSAH association from both designs together. Where they agree, the",
    "case is stronger; where they differ, the disagreement points to which bias is at work."],
    head="SYNTHESIS ACROSS DESIGNS", fc=PUR_BG, ec=PUR, fs=8.8, lh=3.1)

fig.savefig("/Volumes/Niels 2/MIMIC/estrogen-asah-two-designs/figures/study_methods.png",
            dpi=180, bbox_inches="tight", facecolor="white")
print("saved")
