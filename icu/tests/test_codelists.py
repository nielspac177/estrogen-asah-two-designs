"""P1 — every codelist loads, carries expected codes, and matchers behave."""

import pytest

from estrogen_asah_dci.codelists import (
    any_pattern,
    available_codelists,
    code_matches,
    load_all,
    load_codelist,
)

EXPECTED = {
    "sah",
    "aneurysm",
    "aneurysm_procedures",
    "vasospasm",
    "dci_procedures",
    "hrt",
    "eicu_asah_strings",
}


def test_all_expected_codelists_present():
    assert EXPECTED.issubset(set(available_codelists()))


def test_all_codelists_load():
    lists = load_all()
    for name, cl in lists.items():
        assert cl.name == name
        assert cl.description


@pytest.mark.parametrize(
    "name,field,expected",
    [
        ("sah", "icd9", "430"),
        ("sah", "icd10", "I60"),
        ("aneurysm", "icd9", "4373"),
        ("vasospasm", "icd10", "I67848"),
        ("dci_procedures", "eicu_treatment_patterns", "cerebral angioplasty"),
    ],
)
def test_expected_codes_present(name, field, expected):
    assert expected in load_codelist(name).field(field)


def test_code_matches_prefix_and_dot_insensitive():
    sah = load_codelist("sah")
    assert code_matches("430", sah.field("icd9"))
    assert code_matches("I60.7", sah.field("icd10"))   # dotted form
    assert code_matches("I609", sah.field("icd10"))
    assert not code_matches("431", sah.field("icd9"))   # ICH, not SAH


def test_any_pattern_matches_eicu_strings():
    vaso = load_codelist("vasospasm")
    pats = vaso.field("eicu_diagnosis_patterns")
    s = "hemorrhagic stroke|subarachnoid hemorrhage|with vasospasm"
    assert any_pattern(s, pats)
    assert not any_pattern("neurologic|stroke|ischemic", pats)


def test_nimodipine_is_not_a_dci_marker():
    # standard of care for all aSAH — must not be in the DCI rescue list
    dci = load_codelist("dci_procedures")
    assert not any_pattern(
        "neurologic|ICH/ cerebral infarct|nimodipine",
        dci.field("eicu_treatment_patterns"),
    )
