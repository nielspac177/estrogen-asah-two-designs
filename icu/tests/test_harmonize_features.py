"""P5–P6 — pooling + derived features on the synthetic cohort."""

from estrogen_asah_dci import synthetic as syn
from estrogen_asah_dci.extract import eicu, mimic_iv
from estrogen_asah_dci.features.derive import add_features, dci_composite
from estrogen_asah_dci.harmonize.combine import combine


def _pooled(tmp_path):
    con = syn.load_mimic_duckdb()
    m = mimic_iv.build_cohort(con)
    d = eicu.build_cohort(str(syn.write_eicu(tmp_path / "eicu")))
    return combine(m, d)


def test_combine_row_count_and_uniqueness(tmp_path):
    pooled = _pooled(tmp_path)
    # 9 aSAH per source, namespaced ids -> 18 unique rows
    assert len(pooled) == 2 * len(syn.COHORT_IDS)
    assert pooled["patient_id"].is_unique
    assert set(pooled["source"]) == {"mimic_iv", "eicu"}


def test_menopausal_stratum_boundaries():
    con = syn.load_mimic_duckdb()
    df = add_features(mimic_iv.build_cohort(con))
    df = df.assign(_id=[int(p.split(":")[1]) for p in df["patient_id"]]).set_index("_id")
    assert df.loc[8, "menopausal_stratum"] == "premenopausal"   # age 50
    assert df.loc[7, "menopausal_stratum"] == "postmenopausal"  # age 51 (boundary)
    assert df.loc[5, "menopausal_stratum"] == "male"


def test_dci_composite_matches_specs():
    con = syn.load_mimic_duckdb()
    df = add_features(mimic_iv.build_cohort(con))
    df = df.assign(_id=[int(p.split(":")[1]) for p in df["patient_id"]]).set_index("_id")
    assert sorted(df.index[df["dci_composite"].fillna(False)]) == syn.DCI_IDS


def test_features_present_and_pooled_dci(tmp_path):
    pooled = add_features(_pooled(tmp_path))
    for col in ["menopausal_stratum", "dci_composite", "female", "age_c"]:
        assert col in pooled.columns
    # each source contributes the same 5 DCI+ cases -> 10 pooled
    assert int(pooled["dci_composite"].fillna(False).sum()) == 2 * len(syn.DCI_IDS)


def test_dci_composite_is_or_of_components():
    con = syn.load_mimic_duckdb()
    df = mimic_iv.build_cohort(con)
    expected = (
        df["vasospasm_dx"].fillna(False)
        | df["dci_procedure"].fillna(False)
        | df["delayed_infarction"].fillna(False)
    )
    assert (dci_composite(df).fillna(False) == expected).all()
