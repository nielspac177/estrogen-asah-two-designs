"""P9 — dashboard builds, is self-contained, and leaks no patient-level data."""

import json
import sys
from pathlib import Path

from estrogen_asah_dci import synthetic as syn
from estrogen_asah_dci.analysis.evalue import evalue_or
from estrogen_asah_dci.analysis.models import crude_dci_rates, primary_logistic
from estrogen_asah_dci.features.derive import add_features
from estrogen_asah_dci.report.tables import table_one
from estrogen_asah_dci.viz import figures

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "dashboard"))
from build import build_dashboard  # noqa: E402


def _make_outputs(out: Path):
    cohort = add_features(syn.simulate_cohort(n=600, seed=5))
    out.mkdir(parents=True, exist_ok=True)
    rates = crude_dci_rates(cohort)
    rates.to_csv(out / "crude_rates.csv", index=False)
    table_one(cohort).to_csv(out / "table1.csv")
    prim = primary_logistic(cohort)
    base = float(cohort["dci_composite"].astype("boolean").astype("float").mean())
    ev = evalue_or(prim.odds_ratio, base, prim.ci_low, prim.ci_high)
    (out / "results.json").write_text(json.dumps({
        "primary_post_vs_pre": prim.as_dict(),
        "evalue": {"point": ev.point, "ci_bound": ev.ci_bound},
    }))
    (out / "cohort_meta.json").write_text(json.dumps({"mode": "synthetic", "n": len(cohort)}))
    figures.rates_plot(rates, out / "figures")
    return cohort


def test_dashboard_builds_and_is_self_contained(tmp_path):
    out = tmp_path / "outputs"
    _make_outputs(out)
    idx = build_dashboard(out, tmp_path / "dist")
    html = idx.read_text()
    assert idx.stat().st_size > 10_000
    assert "Key results" in html and "SYNTHETIC" in html
    assert 'src="http' not in html          # no external assets (offline/reproducible)
    assert "data:image/png;base64," in html  # figure embedded


def test_dashboard_has_no_patient_ids(tmp_path):
    out = tmp_path / "outputs"
    _make_outputs(out)
    html = build_dashboard(out, tmp_path / "dist").read_text()
    # only aggregates should appear — never namespaced patient identifiers
    assert "mimic_iv:" not in html
    assert "eicu:uid" not in html
