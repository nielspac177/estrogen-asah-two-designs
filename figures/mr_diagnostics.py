"""eFigure. MR diagnostic plots (STROBE-MR): scatter, funnel, leave-one-out, and the
positive-control scatter. Reads per-SNP arrays from mr/outputs/figure_data.json
(produced by mr/scripts/export_figure_data.py)."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": ["Avenir Next", "Arial"], "font.size": 9.5})
BLUE = "#0B6FB8"; GREEN = "#0E8F6E"; ORANGE = "#C77A0A"; RED = "#c0392b"; MUTE = "#5b6b7a"

ROOT = Path(__file__).resolve().parents[1]
d = json.loads((ROOT / "mr" / "outputs" / "figure_data.json").read_text())


def arrs(key):
    p = d[key]
    return (np.array(p["bx"]), np.array(p["sx"]), np.array(p["by"]), np.array(p["sy"]))


def ivw_slope(bx, by, sy):
    w = 1 / sy**2
    return np.sum(w * bx * by) / np.sum(w * bx**2)


def egger(bx, by, sy):
    s = np.sign(bx); bxo, byo = bx * s, by * s
    w = 1 / sy**2
    X = np.column_stack([np.ones_like(bxo), bxo])
    W = np.diag(w)
    b = np.linalg.solve(X.T @ W @ X, X.T @ W @ byo)
    return b[0], b[1]           # intercept, slope


fig, axes = plt.subplots(2, 2, figsize=(11, 9))
fig.suptitle("Mendelian randomization diagnostics: age at natural menopause",
             fontsize=13, fontweight="bold")

# ---- A. scatter (primary), SNPs oriented to positive exposure effect ----
ax = axes[0, 0]
bx, sx, by, sy = arrs("primary_asah")
sl = ivw_slope(bx, by, sy); it, es = egger(bx, by, sy)
s = np.sign(bx); bxo, byo = bx * s, by * s
ax.errorbar(bxo, byo, xerr=sx, yerr=sy, fmt="o", ms=3, color=BLUE, ecolor="#b9cfe0",
            elinewidth=0.7, alpha=0.9, zorder=2)
xg = np.linspace(0, bxo.max() * 1.05, 50)
ax.plot(xg, sl * xg, color=RED, lw=1.6, label="IVW")
ax.plot(xg, it + es * xg, color="#7B4FBE", lw=1.4, ls="--", label="MR-Egger")
ax.axhline(0, color="#ddd", lw=0.8)
ax.set_xlabel("SNP effect on age at menopause"); ax.set_ylabel("SNP effect on aSAH")
ax.set_title("A. Scatter", loc="left", fontweight="bold", fontsize=10.5)
ax.legend(fontsize=8, frameon=False)

# ---- B. funnel (primary) ----
ax = axes[0, 1]
iv = by / bx                     # per-SNP Wald ratio
strength = np.abs(bx) / sy       # instrument strength proxy
ax.scatter(iv, strength, s=14, color=BLUE, alpha=0.8)
ax.axvline(sl, color=RED, lw=1.4, label="IVW estimate")
ax.set_xlabel("SNP causal estimate (log OR)"); ax.set_ylabel("Instrument strength |bx|/se")
ax.set_title("B. Funnel", loc="left", fontweight="bold", fontsize=10.5)
ax.set_xlim(np.percentile(iv, 2), np.percentile(iv, 98))
ax.legend(fontsize=8, frameon=False)

# ---- C. leave-one-out (primary) ----
ax = axes[1, 0]
loo = np.array(d["primary_asah"]["loo_or"])
y = np.arange(len(loo))
ax.scatter(loo, y, s=10, color=BLUE, alpha=0.7)
ax.axvline(np.exp(sl), color=RED, lw=1.4, label="all SNPs")
ax.axvline(1.0, color=MUTE, ls="--", lw=0.9)
ax.set_xlabel("IVW odds ratio, leaving out each SNP"); ax.set_ylabel("SNP (dropped)")
ax.set_title("C. Leave-one-out", loc="left", fontweight="bold", fontsize=10.5)
ax.legend(fontsize=8, frameon=False)

# ---- D. positive-control scatter (oriented) ----
ax = axes[1, 1]
bx2, sx2, by2, sy2 = arrs("poscontrol_bc")
sl2 = ivw_slope(bx2, by2, sy2)
s2 = np.sign(bx2); bxo2, byo2 = bx2 * s2, by2 * s2
ax.errorbar(bxo2, byo2, xerr=sx2, yerr=sy2, fmt="o", ms=3, color=ORANGE, ecolor="#e6c79a",
            elinewidth=0.7, alpha=0.9, zorder=2)
xg2 = np.linspace(0, bxo2.max() * 1.05, 50)
ax.plot(xg2, sl2 * xg2, color=RED, lw=1.6, label="IVW")
ax.axhline(0, color="#ddd", lw=0.8)
ax.set_xlabel("SNP effect on age at menopause"); ax.set_ylabel("SNP effect on breast cancer")
ax.set_title("D. Positive control (breast cancer)", loc="left", fontweight="bold", fontsize=10.5)
ax.legend(fontsize=8, frameon=False)

fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(ROOT / "figures" / "mr_diagnostics.png", dpi=170, bbox_inches="tight", facecolor="white")
print("saved")
