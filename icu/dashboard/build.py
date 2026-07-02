"""Build a self-contained results dashboard (P9).

Reads the analysis outputs and emits a single ``index.html`` with no external
dependencies: inline CSS (dark/light via prefers-color-scheme), KPI cards, and
figures embedded as base64. Only aggregate results are included — never PHI.
"""

from __future__ import annotations

import base64
import json
from pathlib import Path

import pandas as pd

_CSS = """
:root { color-scheme: light dark; --bg:#ffffff; --fg:#1a1a1a; --muted:#666;
  --card:#f5f7fa; --border:#e2e8f0; --accent:#0072B2; }
@media (prefers-color-scheme: dark) { :root { --bg:#0f1419; --fg:#e6edf3;
  --muted:#9aa7b4; --card:#161b22; --border:#30363d; --accent:#58a6ff; } }
* { box-sizing:border-box; } body { margin:0; background:var(--bg); color:var(--fg);
  font:16px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif; }
main { max-width:1000px; margin:0 auto; padding:1.5rem; }
h1 { font-size:1.5rem; margin:.2rem 0; } h2 { font-size:1.1rem; margin:2rem 0 .6rem; }
.banner { padding:.6rem 1rem; border-radius:8px; font-weight:600; }
.banner.syn { background:#E69F0022; border:1px solid #E69F00; }
.banner.real { background:#009E7322; border:1px solid #009E73; }
.grid { display:grid; gap:1rem; grid-template-columns:repeat(auto-fit,minmax(180px,1fr)); }
.card { background:var(--card); border:1px solid var(--border); border-radius:10px; padding:1rem; }
.card .v { font-size:1.6rem; font-weight:700; color:var(--accent); }
.card .l { color:var(--muted); font-size:.85rem; }
figure { margin:0 0 1rem; } figure img { width:100%; height:auto; border:1px solid var(--border); border-radius:8px; }
figcaption { color:var(--muted); font-size:.85rem; margin-top:.3rem; }
table { border-collapse:collapse; width:100%; font-size:.9rem; }
th,td { border:1px solid var(--border); padding:.35rem .6rem; text-align:right; }
th:first-child,td:first-child { text-align:left; }
.note { color:var(--muted); font-size:.85rem; border-left:3px solid var(--accent); padding:.4rem .8rem; }
"""


def _b64(path: Path) -> str:
    return base64.b64encode(path.read_bytes()).decode()


def _kpis(results: dict, meta: dict, rates: pd.DataFrame) -> list[tuple[str, str]]:
    prim = results.get("primary_post_vs_pre", {})
    kp = [("Cohort (n)", f"{meta.get('n', 0):,}")]
    if prim.get("converged"):
        kp.append(("Primary OR (post vs pre)",
                   f"{prim['or']:.2f} ({prim['ci_low']:.2f}–{prim['ci_high']:.2f})"))
    ev = results.get("evalue")
    if ev:
        kp.append(("E-value (point)", f"{ev['point']:.2f}"))
    if len(rates):
        overall = rates["events"].sum() / rates["n"].sum()
        kp.append(("DCI / vasospasm", f"{overall*100:.0f}%"))
    return kp


def build_dashboard(outputs: Path, dist: Path) -> Path:
    dist.mkdir(parents=True, exist_ok=True)
    results = json.loads((outputs / "results.json").read_text())
    meta = json.loads((outputs / "cohort_meta.json").read_text())
    rates = pd.read_csv(outputs / "crude_rates.csv")
    table1 = pd.read_csv(outputs / "table1.csv", index_col=0)

    synthetic = meta.get("mode") != "real"
    banner_cls = "syn" if synthetic else "real"
    banner_txt = ("SYNTHETIC DATA — illustrative only, no PhysioNet data"
                  if synthetic else "Real data (MIMIC-IV + eICU)")

    cards = "".join(
        f'<div class="card"><div class="v">{v}</div><div class="l">{lbl}</div></div>'
        for lbl, v in _kpis(results, meta, rates)
    )

    fig_dir = outputs / "figures"
    figs = [
        ("cohort_flow.png", "CONSORT-style cohort selection flow"),
        ("dci_rates.png", "Crude DCI/vasospasm rate by menopausal stratum"),
        ("forest.png", "Adjusted odds ratios for DCI (forest plot)"),
        ("love.png", "Covariate balance before and after IPW weighting"),
    ]
    fig_html = ""
    for fname, alt in figs:
        fp = fig_dir / fname
        if fp.exists():
            fig_html += (f'<figure><img alt="{alt}" src="data:image/png;base64,{_b64(fp)}">'
                         f'<figcaption>{alt}</figcaption></figure>')

    t1_html = table1.to_html(border=0, classes="t1", justify="left")

    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Estrogen &amp; DCI after aSAH — results</title>
<style>{_CSS}</style></head>
<body><main>
<h1>Estrogen state &amp; delayed cerebral ischemia after aSAH</h1>
<p class="banner {banner_cls}" role="status">{banner_txt}</p>
<h2>Key results</h2>
<section class="grid" aria-label="Key metrics">{cards}</section>
<h2>Figures</h2>
{fig_html}
<h2>Table 1 — baseline characteristics by stratum</h2>
{t1_html}
<h2>Interpretation</h2>
<p class="note">Hypothesis-generating pilot. A protective point estimate is consistent
with the estrogen-neuroprotection prior but is not confirmatory: the cohort is modest,
the exposure is a menopausal-state proxy confounded by age, and residual confounding is
plausible (see the E-value). Interpret directionally, not causally.</p>
</main></body></html>"""

    out = dist / "index.html"
    out.write_text(html, encoding="utf-8")
    return out
