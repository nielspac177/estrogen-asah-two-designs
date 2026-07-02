"""Publication figures (P8). Matplotlib, Agg backend, colourblind-safe palette.

Every function takes plain data + an output directory and returns the saved path,
so figures are reproducible from the analysis outputs alone.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Okabe-Ito (colourblind-safe)
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#E69F00", "#009E73", "#999999"
plt.rcParams.update({"figure.dpi": 130, "font.size": 10, "svg.fonttype": "none"})


def _save(fig, out_dir: Path, name: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / name
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def forest_plot(rows: list[dict], out_dir: Path, name: str = "forest.png") -> Path:
    """rows: [{label, or, ci_low, ci_high}] — ORs on a log axis with a null line at 1."""
    rows = [r for r in rows if r.get("or") == r.get("or") and r.get("or")]  # drop NaN
    labels = [r["label"] for r in rows]
    y = range(len(rows))
    fig, ax = plt.subplots(figsize=(6.4, 0.6 * len(rows) + 1.2))
    for i, r in zip(y, rows, strict=True):
        ax.plot([r["ci_low"], r["ci_high"]], [i, i], color=GREY, lw=1.8, zorder=1)
        ax.scatter([r["or"]], [i], color=BLUE, s=42, zorder=2)
    ax.axvline(1.0, color=ORANGE, ls="--", lw=1)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio for DCI (log scale)")
    ax.invert_yaxis()
    ax.set_title("DCI/vasospasm — adjusted odds ratios")
    return _save(fig, out_dir, name)


def rates_plot(rates, out_dir: Path, name: str = "dci_rates.png") -> Path:
    """rates: DataFrame with menopausal_stratum, n, rate."""
    order = ["premenopausal", "postmenopausal", "male"]
    r = rates.set_index("menopausal_stratum").reindex(order).dropna(how="all")
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    bars = ax.bar(r.index, r["rate"] * 100, color=[BLUE, GREEN, GREY][: len(r)])
    for b, (_, row) in zip(bars, r.iterrows(), strict=True):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.5,
                f"{row['rate']*100:.0f}%\n(n={int(row['n'])})",
                ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("DCI / vasospasm (%)")
    ax.set_title("Crude DCI rate by menopausal stratum")
    ax.margins(y=0.15)
    return _save(fig, out_dir, name)


def love_plot(balance, out_dir: Path, name: str = "love.png") -> Path:
    """balance: DataFrame with covariate, smd_unweighted, smd_weighted."""
    fig, ax = plt.subplots(figsize=(5.2, 0.5 * len(balance) + 1.4))
    y = range(len(balance))
    ax.scatter(balance["smd_unweighted"], y, color=GREY, label="unweighted")
    ax.scatter(balance["smd_weighted"], y, color=BLUE, label="IPW-weighted")
    ax.axvline(0.1, color=ORANGE, ls="--", lw=1, label="|SMD|=0.1")
    ax.set_yticks(list(y))
    ax.set_yticklabels(balance["covariate"])
    ax.set_xlabel("|Standardized mean difference|")
    ax.set_title("Covariate balance (pre vs post, women)")
    ax.legend(fontsize=8)
    return _save(fig, out_dir, name)


def cohort_flow(steps: list[tuple[str, int]], out_dir: Path, name: str = "cohort_flow.png") -> Path:
    """steps: [(label, n)] top-to-bottom CONSORT-style boxes."""
    fig, ax = plt.subplots(figsize=(5.2, 0.9 * len(steps) + 0.6))
    ax.axis("off")
    n = len(steps)
    for i, (label, count) in enumerate(steps):
        yv = n - i
        ax.add_patch(plt.Rectangle((0.1, yv - 0.4), 0.8, 0.7, fill=True,
                                   facecolor="#EAF2F8", edgecolor=BLUE))
        ax.text(0.5, yv - 0.05, f"{label}\nn = {count:,}", ha="center", va="center", fontsize=9)
        if i < n - 1:
            ax.annotate("", xy=(0.5, yv - 0.45), xytext=(0.5, yv - 0.6),
                        arrowprops=dict(arrowstyle="->", color=GREY))
    ax.set_xlim(0, 1)
    ax.set_ylim(0.3, n + 0.6)
    ax.set_title("Cohort selection")
    return _save(fig, out_dir, name)
