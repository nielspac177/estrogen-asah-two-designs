"""P3 — MIMIC-IV extractor reproduces SPECS ground truth from encoded tables."""

import pytest

from estrogen_asah_dci import synthetic as syn
from estrogen_asah_dci.extract import mimic_iv


@pytest.fixture(scope="module")
def con():
    return syn.load_mimic_duckdb()


def _ids(df):
    return sorted(int(p.split(":")[1]) for p in df["patient_id"])


def test_asah_cohort_matches_ground_truth(con):
    cohort = mimic_iv.asah_cohort(con)
    got = sorted(int(x) for x in cohort["subject_id"])
    assert got == syn.COHORT_IDS  # excludes unruptured(9), traumatic(10), unrelated(11)


def test_build_cohort_row_count(con):
    df = mimic_iv.build_cohort(con)
    assert len(df) == len(syn.COHORT_IDS)
    assert _ids(df) == syn.COHORT_IDS


def _indexed(con):
    df = mimic_iv.build_cohort(con)
    return df.assign(_id=[int(p.split(":")[1]) for p in df["patient_id"]]).set_index("_id")


def test_flag_columns(con):
    df = _indexed(con)

    def ids_true(col):
        return sorted(df.index[df[col].fillna(False)])

    assert ids_true("vasospasm_dx") == [1, 5, 8]
    assert ids_true("dci_procedure") == [2, 8]
    assert ids_true("delayed_infarction") == [4]
    assert ids_true("hrt_exposure") == [2, 8, 12]
    assert ids_true("died") == [3]
    assert ids_true("poor_disposition") == [3, 4]      # died OR hospice
    assert ids_true("htn") == [1, 2, 3, 5, 6, 8]
    assert ids_true("smoking") == [1, 4, 5, 8]
    assert ids_true("diabetes") == [2, 6]


def test_demographics_and_modality(con):
    df = _indexed(con)
    assert df.loc[1, "sex"] == "F" and df.loc[1, "age"] == 44
    assert df.loc[5, "sex"] == "M"
    assert df.loc[1, "treatment_modality"] == "clip"
    assert df.loc[2, "treatment_modality"] == "coil"
    assert df.loc[3, "aneurysm_secured"] == False  # noqa: E712  (in via aneurysm dx, not secured)


def test_schema_valid(con):
    from estrogen_asah_dci.harmonize.common_schema import validate_cohort
    validate_cohort(mimic_iv.build_cohort(con))
