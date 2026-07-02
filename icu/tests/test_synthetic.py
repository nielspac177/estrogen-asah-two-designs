"""P2 — synthetic fixtures load, are internally consistent, cover every branch."""

import duckdb
import pandas as pd

from estrogen_asah_dci import synthetic as syn
from estrogen_asah_dci.harmonize.common_schema import (
    COHORT_COLUMNS,
    coerce_cohort,
    empty_cohort,
    validate_cohort,
)


def test_ground_truth_counts():
    assert len(syn.SPECS) == 12
    assert sorted(syn.COHORT_IDS) == [1, 2, 3, 4, 5, 6, 7, 8, 12]  # 9 aSAH
    assert sorted(syn.DCI_IDS) == [1, 2, 4, 5, 8]                   # 5 DCI+


def test_branch_coverage():
    cohort = [s for s in syn.SPECS if s["in_cohort"]]
    women = [s for s in cohort if s["sex"] == "F"]
    assert any(s["age"] < 51 for s in women)            # premenopausal
    assert any(s["age"] >= 51 for s in women)           # postmenopausal
    assert any(s["sex"] == "M" for s in cohort)         # male reference arm
    assert any(s["vasospasm"] for s in cohort)          # DCI via each component
    assert any(s["dci_proc"] for s in cohort)
    assert any(s["infarct"] for s in cohort)
    assert any(s["hrt"] for s in cohort) and any(not s["hrt"] for s in cohort)
    assert any(s["coding"] == "i9" for s in cohort)     # ICD-9 and ICD-10
    assert any(s["coding"] == "i10" for s in cohort)
    # non-cohort patients present (unruptured, traumatic, unrelated)
    assert {9, 10, 11}.issubset({s["id"] for s in syn.SPECS if not s["in_cohort"]})


def test_mimic_tables_shape():
    t = syn.mimic_tables()
    assert set(t) == {
        "patients", "admissions", "diagnoses_icd",
        "procedures_icd", "prescriptions", "d_icd_diagnoses",
    }
    assert len(t["patients"]) == 12
    assert (t["prescriptions"]["drug"].str.contains("Estradiol")).all()


def test_mimic_loads_into_duckdb():
    con = syn.load_mimic_duckdb()
    n = con.execute("SELECT count(*) FROM mimiciv_hosp.patients").fetchone()[0]
    assert n == 12
    # SAH codes present
    m = con.execute(
        "SELECT count(*) FROM mimiciv_hosp.diagnoses_icd WHERE icd_code IN ('I609','430')"
    ).fetchone()[0]
    assert m >= 1


def test_eicu_tables_and_gz_roundtrip(tmp_path):
    d = syn.write_eicu(tmp_path / "eicu")
    con = duckdb.connect()
    n = con.execute(
        f"SELECT count(*) FROM read_csv_auto('{d}/patient.csv.gz')"
    ).fetchone()[0]
    assert n == 12
    # eICU 'with vasospasm' diagnosis string is present
    v = con.execute(
        f"SELECT count(*) FROM read_csv_auto('{d}/diagnosis.csv.gz') "
        "WHERE diagnosisstring ILIKE '%with vasospasm%'"
    ).fetchone()[0]
    assert v >= 1


def test_common_schema_helpers():
    e = empty_cohort()
    assert list(e.columns) == COHORT_COLUMNS
    assert len(e) == 0
    # a minimal valid row
    row = {c: pd.NA for c in COHORT_COLUMNS}
    row.update(source="mimic_iv", patient_id="mimic_iv:1", hospital_id="mimic_iv",
               age=44.0, sex="F")
    df = coerce_cohort(pd.DataFrame([row]))
    validate_cohort(df)
    assert df["sex"].iloc[0] == "F"
