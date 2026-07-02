"""MIMIC-IV cohort extraction (DuckDB, read-only).

Builds the harmonized aSAH cohort from the ``mimiciv_hosp`` schema. SQL is kept
minimal — it selects the small SAH-related slices — and phenotype logic reuses
the shared, tested codelist matchers so it is identical to the eICU extractor.

Every function accepts a DuckDB connection so it can run against either the real
``mimic4.db`` or the in-memory synthetic database (`synthetic.load_mimic_duckdb`).
"""

from __future__ import annotations

import duckdb
import pandas as pd

from ..codelists import load_codelist
from ..harmonize.common_schema import coerce_cohort

SCHEMA = "mimiciv_hosp"


def connect(db_path: str, *, read_only: bool = True) -> duckdb.DuckDBPyConnection:
    return duckdb.connect(db_path, read_only=read_only)


# ---- SQL condition builders (from codelists) ----

def _prefix_or(col: str, prefixes: list[str]) -> str:
    if not prefixes:
        return "FALSE"
    return "(" + " OR ".join(f"{col} LIKE '{p}%'" for p in prefixes) + ")"


def _icd_cond(code_col: str, ver_col: str, icd9: list[str], icd10: list[str]) -> str:
    parts = []
    if icd9:
        parts.append(f"({ver_col}=9 AND {_prefix_or(code_col, icd9)})")
    if icd10:
        parts.append(f"({ver_col}=10 AND {_prefix_or(code_col, icd10)})")
    return "(" + " OR ".join(parts) + ")" if parts else "FALSE"


# ---- cohort building blocks ----

def sah_admissions(con, schema: str = SCHEMA) -> pd.DataFrame:
    """Subject/admission pairs with a nontraumatic SAH diagnosis (3.2)."""
    sah = load_codelist("sah")
    cond = _icd_cond("icd_code", "icd_version", sah.field("icd9"), sah.field("icd10"))
    return con.execute(
        f"SELECT DISTINCT subject_id, hadm_id FROM {schema}.diagnoses_icd WHERE {cond}"
    ).df()


def _subjects_with_dx(con, schema, icd9, icd10) -> set[int]:
    cond = _icd_cond("icd_code", "icd_version", icd9, icd10)
    return set(
        con.execute(
            f"SELECT DISTINCT subject_id FROM {schema}.diagnoses_icd WHERE {cond}"
        ).df()["subject_id"]
    )


def _secured_modality(con, schema) -> pd.DataFrame:
    """Subjects with an aneurysm-securing procedure + modality (clip/coil) (3.3)."""
    proc = load_codelist("aneurysm_procedures")
    clip9 = [c for c in proc.field("icd9_proc") if c in ("3951", "3952")]
    coil9 = [c for c in proc.field("icd9_proc") if c in ("3972", "3975", "3976")]
    clip10 = [c for c in proc.field("icd10_proc") if c.startswith("03V")]
    coil10 = [c for c in proc.field("icd10_proc") if c.startswith("03L")]
    clip = _icd_cond("icd_code", "icd_version", clip9, clip10)
    coil = _icd_cond("icd_code", "icd_version", coil9, coil10)
    return con.execute(
        f"""
        SELECT subject_id,
               CASE WHEN {clip} THEN 'clip' WHEN {coil} THEN 'coil' END AS modality
        FROM {schema}.procedures_icd
        WHERE {clip} OR {coil}
        """
    ).df().dropna(subset=["modality"]).drop_duplicates("subject_id")


def asah_cohort(con, schema: str = SCHEMA) -> pd.DataFrame:
    """aSAH = SAH AND (aneurysm dx OR securing procedure). One row per subject (3.3)."""
    sah = sah_admissions(con, schema)
    aneurysm = load_codelist("aneurysm")
    aneu_subjects = _subjects_with_dx(con, schema, aneurysm.field("icd9"), aneurysm.field("icd10"))
    secured = _secured_modality(con, schema)
    secured_subjects = set(secured["subject_id"])

    keep = sah[sah["subject_id"].isin(aneu_subjects | secured_subjects)].copy()
    # first qualifying admission per subject
    keep = keep.sort_values(["subject_id", "hadm_id"]).drop_duplicates("subject_id")
    keep = keep.merge(secured, on="subject_id", how="left")
    keep["aneurysm_secured"] = keep["subject_id"].isin(secured_subjects)
    return keep.rename(columns={"modality": "treatment_modality"})


def build_cohort(con, schema: str = SCHEMA) -> pd.DataFrame:
    """Full harmonized MIMIC-IV cohort (3.10)."""
    cohort = asah_cohort(con, schema)
    if cohort.empty:
        from ..harmonize.common_schema import empty_cohort
        return empty_cohort()

    subjects = tuple(int(x) for x in cohort["subject_id"])
    in_list = f"({','.join(str(s) for s in subjects)})" if subjects else "(NULL)"

    patients = con.execute(
        f"SELECT subject_id, gender, anchor_age FROM {schema}.patients "
        f"WHERE subject_id IN {in_list}"
    ).df()
    admissions = con.execute(
        f"SELECT subject_id, hadm_id, admittime, dischtime, hospital_expire_flag, "
        f"discharge_location FROM {schema}.admissions WHERE subject_id IN {in_list}"
    ).df()
    dx = con.execute(
        f"SELECT subject_id, icd_code, icd_version FROM {schema}.diagnoses_icd "
        f"WHERE subject_id IN {in_list}"
    ).df()
    proc = con.execute(
        f"SELECT subject_id, icd_code, icd_version FROM {schema}.procedures_icd "
        f"WHERE subject_id IN {in_list}"
    ).df()
    rx = con.execute(
        f"SELECT subject_id, drug FROM {schema}.prescriptions WHERE subject_id IN {in_list}"
    ).df()

    vaso = load_codelist("vasospasm")
    dci = load_codelist("dci_procedures")
    hrt = load_codelist("hrt")
    com = load_codelist("comorbidities")

    def dx_subjects(icd9, icd10) -> set:
        if dx.empty:
            return set()
        cond9 = dx["icd_version"].eq(9) & _match(dx["icd_code"], icd9)
        cond10 = dx["icd_version"].eq(10) & _match(dx["icd_code"], icd10)
        return set(dx.loc[cond9 | cond10, "subject_id"])

    def proc_subjects(icd9, icd10) -> set:
        if proc.empty:
            return set()
        cond9 = proc["icd_version"].eq(9) & _match(proc["icd_code"], icd9)
        cond10 = proc["icd_version"].eq(10) & _match(proc["icd_code"], icd10)
        return set(proc.loc[cond9 | cond10, "subject_id"])

    vaso_s = dx_subjects([], vaso.field("icd10"))
    dciproc_s = proc_subjects(dci.field("icd9_proc"), dci.field("icd10_proc"))
    infarct_s = dx_subjects(
        dci.field("delayed_infarction_icd9"), dci.field("delayed_infarction_icd10")
    )
    htn_s = dx_subjects(com.field("htn_icd9"), com.field("htn_icd10"))
    smk_s = dx_subjects(com.field("smoking_icd9"), com.field("smoking_icd10"))
    dm_s = dx_subjects(com.field("diabetes_icd9"), com.field("diabetes_icd10"))

    hrt_pat = "|".join(hrt.field("drug_patterns"))
    if rx.empty:
        hrt_s = set()
    else:
        hrt_s = set(rx.loc[rx["drug"].str.contains(hrt_pat, case=False, na=False), "subject_id"])

    adm = cohort[["subject_id", "hadm_id", "treatment_modality", "aneurysm_secured"]].merge(
        admissions, on=["subject_id", "hadm_id"], how="left"
    )
    adm = adm.merge(patients, on="subject_id", how="left")

    delta = pd.to_datetime(adm["dischtime"]) - pd.to_datetime(adm["admittime"])
    los = delta.dt.total_seconds() / 86400.0
    disloc = adm["discharge_location"].fillna("").str.upper()

    out = pd.DataFrame({
        "source": "mimic_iv",
        "patient_id": "mimic_iv:" + adm["subject_id"].astype(str),
        "hospital_id": "mimic_iv",
        "age": adm["anchor_age"].astype("Float64"),
        "sex": adm["gender"],
        "vasospasm_dx": adm["subject_id"].isin(vaso_s),
        "dci_procedure": adm["subject_id"].isin(dciproc_s),
        "delayed_infarction": adm["subject_id"].isin(infarct_s),
        "hrt_exposure": adm["subject_id"].isin(hrt_s),
        "aneurysm_secured": adm["aneurysm_secured"].fillna(False),
        "treatment_modality": adm["treatment_modality"],
        "htn": adm["subject_id"].isin(htn_s),
        "smoking": adm["subject_id"].isin(smk_s),
        "diabetes": adm["subject_id"].isin(dm_s),
        # GCS from chartevents is a real-data refinement (see ADR-0001); MIMIC has no APACHE IV
        "severity_gcs": pd.NA,
        "apache_score": pd.NA,
        "died": adm["hospital_expire_flag"].fillna(0).astype(int).astype(bool),
        "poor_disposition": adm["hospital_expire_flag"].fillna(0).astype(int).astype(bool)
                            | disloc.str.contains("HOSPICE"),
        "icu_los_days": los.astype("Float64"),
    })
    return coerce_cohort(out)


def _match(codes: pd.Series, prefixes: list[str]) -> pd.Series:
    """Vectorized prefix match (dot-insensitive) for a code column."""
    if not prefixes:
        return pd.Series(False, index=codes.index)
    norm = codes.astype(str).str.replace(".", "", regex=False).str.upper()
    pat = "^(?:" + "|".join(p.replace(".", "").upper() for p in prefixes) + ")"
    return norm.str.contains(pat, regex=True, na=False)
