"""Render every display to figures/final/ as numbered PNG + PDF + SVG.

Runs each figure script in-process and redirects its savefig to the numbered final
name in all three formats (vector PDF/SVG come straight from matplotlib). The
ordering prefix (1., 2., ...) keeps the files sorted as they appear in the paper.
"""
import runpy
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.figure as mfig
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT / "figures"
FINAL = FIGDIR / "final"
FINAL.mkdir(exist_ok=True)

# render order: main figures, table, supplementary figures, visual summary
MANIFEST = [
    ("study_methods.py", "1.Figure_1_study_design"),
    ("cohort_flow.py", "2.Figure_2_participant_flow"),
    ("results_forest_table.py", "3.Figure_3_results_forest_table"),
    ("mr_cascade.py", "4.Figure_4_cascade_forest"),
    ("table1_figure.py", "5.Table_1_baseline_characteristics"),
    ("dag.py", "6.eFigure_1_causal_diagrams"),
    ("overlap_nonidentifiability.py", "7.eFigure_2_nonidentifiability"),
    ("mr_diagnostics.py", "8.eFigure_3_mr_diagnostics"),
    ("observational_diagnostics.py", "9.eFigure_4_observational_diagnostics"),
    ("graphical_abstract.py", "10.Visual_graphical_abstract"),
    ("study_workflow.py", "11.Visual_workflow_overview"),
]

_orig = mfig.Figure.savefig
_state = {"name": None, "count": 0}


def _patched(self, fname=None, **kw):
    # Each figure script renders a single figure; ignore repeat savefig calls
    # (e.g. a script's own png/pdf/svg loop) and emit all three formats once.
    kw.pop("format", None)
    _state["count"] += 1
    if _state["count"] > 1:
        return
    for ext in ("png", "pdf", "svg"):
        _orig(self, str(FINAL / f"{_state['name']}.{ext}"), **kw)


mfig.Figure.savefig = _patched

for script, name in MANIFEST:
    _state["name"], _state["count"] = name, 0
    try:
        runpy.run_path(str(FIGDIR / script), run_name="__main__")
        print(f"  ok  {name:44s} <- {script}")
    except Exception as e:  # noqa: BLE001
        print(f"  FAIL {name:43s} <- {script}: {e}")
    plt.close("all")

mfig.Figure.savefig = _orig
n = len(list(FINAL.glob("*.png")))
print(f"\nfinal/ now holds {n} PNG (+ matching PDF/SVG). Path: {FINAL}")
