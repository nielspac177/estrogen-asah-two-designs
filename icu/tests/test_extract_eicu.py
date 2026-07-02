"""P4 — eICU extractor reproduces SPECS ground truth from encoded CSVs."""

import pytest

from estrogen_asah_dci import synthetic as syn
from estrogen_asah_dci.extract import eicu
from estrogen_asah_dci.harmonize.common_schema import validate_cohort


@pytest.fixture(scope="module")
def eicu_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("eicu_data")
    return str(syn.write_eicu(d))


def _indexed(eicu_dir):
    df = eicu.build_cohort(eicu_dir)
    return df.assign(_id=[int(p.split("-")[1]) for p in df["patient_id"]]).set_index("_id")


def test_cohort_matches_ground_truth(eicu_dir):
    df = _indexed(eicu_dir)
    assert sorted(df.index) == syn.COHORT_IDS


def test_flag_columns(eicu_dir):
    df = _indexed(eicu_dir)

    def ids_true(col):
        return sorted(df.index[df[col].fillna(False)])

    assert ids_true("vasospasm_dx") == [1, 5, 8]
    assert ids_true("dci_procedure") == [2, 8]
    assert ids_true("delayed_infarction") == [4]
    assert ids_true("hrt_exposure") == [2, 8, 12]
    assert ids_true("died") == [3]
    assert ids_true("poor_disposition") == [3, 4]
    assert ids_true("htn") == [1, 2, 3, 5, 6, 8]
    assert ids_true("smoking") == [1, 4, 5, 8]
    assert ids_true("diabetes") == [2, 6]


def test_demographics_modality_severity(eicu_dir):
    df = _indexed(eicu_dir)
    assert df.loc[7, "age"] == 51 and df.loc[7, "sex"] == "F"
    assert df.loc[5, "sex"] == "M"
    assert df.loc[1, "treatment_modality"] == "clip"
    assert df.loc[2, "treatment_modality"] == "coil"
    assert df["apache_score"].notna().all()          # eICU carries APACHE severity


def test_dci_ids_consistent_with_specs(eicu_dir):
    df = _indexed(eicu_dir)
    dci = df["vasospasm_dx"] | df["dci_procedure"] | df["delayed_infarction"]
    assert sorted(df.index[dci]) == syn.DCI_IDS


def test_schema_valid(eicu_dir):
    validate_cohort(eicu.build_cohort(eicu_dir))
