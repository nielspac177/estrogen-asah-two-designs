"""eICU-CRD cohort extraction (gzipped CSVs via DuckDB).

Mirrors the MIMIC-IV extractor and produces the same harmonized schema. Files
are read by explicit name (never globbed) so macOS AppleDouble ``._*`` sidecars
on exFAT volumes are ignored. Phenotype logic reuses the shared codelist
matchers, so definitions are identical across sources.
"""

from __future__ import annotations

import duckdb
import pandas as pd

from ..codelists import any_pattern, load_codelist
from ..harmonize.common_schema import coerce_cohort, empty_cohort

_TABLES = {
    "patient": "patient.csv.gz",
    "diagnosis": "diagnosis.csv.gz",
    "admissionDx": "admissionDx.csv.gz",
    "treatment": "treatment.csv.gz",
    "admissionDrug": "admissionDrug.csv.gz",
    "apachePatientResult": "apachePatientResult.csv.gz",
}


def _read(con, eicu_dir, table: str, cols: str = "*") -> pd.DataFrame:
    path = f"{eicu_dir}/{_TABLES[table]}"
    return con.execute(f"SELECT {cols} FROM read_csv_auto('{path}', all_varchar=true)").df()


def _parse_age(val) -> float | None:
    if val is None or str(val).strip() == "":
        return None
    s = str(val).strip()
    return 90.0 if s.startswith(">") else float(s)


def asah_stays(diagnosis: pd.DataFrame, admitdx: pd.DataFrame) -> pd.DataFrame:
    """Identify aSAH unit stays via admission dx / diagnosis strings (4.1)."""
    cl = load_codelist("eicu_asah_strings")
    inc_admit = cl.field("admitdx_include")
    inc_diag = cl.field("diagnosis_include_patterns")
    exc = cl.field("exclude_patterns")

    stays: set = set()
    if not admitdx.empty:
        m = admitdx["admitdxname"].apply(
            lambda s: any(inc.lower() in str(s).lower() for inc in inc_admit)
            and not any_pattern(s, exc)
        )
        stays |= set(admitdx.loc[m, "patientunitstayid"].astype(int))
    if not diagnosis.empty:
        m = diagnosis["diagnosisstring"].apply(
            lambda s: any_pattern(s, inc_diag) and not any_pattern(s, exc)
        )
        stays |= set(diagnosis.loc[m, "patientunitstayid"].astype(int))
    return pd.DataFrame({"patientunitstayid": sorted(stays)}, dtype="int64")


def build_cohort(eicu_dir: str, con=None) -> pd.DataFrame:
    """Full harmonized eICU cohort (4.9)."""
    if con is None:
        con = duckdb.connect()

    patient = _read(con, eicu_dir, "patient")
    diagnosis = _read(con, eicu_dir, "diagnosis")
    admitdx = _read(con, eicu_dir, "admissionDx")
    treatment = _read(con, eicu_dir, "treatment")
    admdrug = _read(con, eicu_dir, "admissionDrug")
    apache = _read(con, eicu_dir, "apachePatientResult")

    stays = asah_stays(diagnosis, admitdx)
    if stays.empty:
        return empty_cohort()

    patient["patientunitstayid"] = patient["patientunitstayid"].astype(int)
    cohort = stays.merge(patient, on="patientunitstayid", how="left")

    # dedup to one stay per person (first qualifying = lowest stay id) (4.2)
    cohort = cohort.sort_values("patientunitstayid").drop_duplicates("uniquepid")

    def stays_matching_diag(patterns) -> set:
        if diagnosis.empty:
            return set()
        m = diagnosis["diagnosisstring"].apply(lambda s: any_pattern(s, patterns))
        return set(diagnosis.loc[m, "patientunitstayid"].astype(int))

    def stays_matching_treat(patterns) -> set:
        if treatment.empty:
            return set()
        m = treatment["treatmentstring"].apply(lambda s: any_pattern(s, patterns))
        return set(treatment.loc[m, "patientunitstayid"].astype(int))

    vaso = load_codelist("vasospasm").field("eicu_diagnosis_patterns")
    dci_treat = load_codelist("dci_procedures").field("eicu_treatment_patterns")
    com = load_codelist("comorbidities")
    hrt_pat = load_codelist("hrt").field("drug_patterns")

    vaso_s = stays_matching_diag(vaso)
    dci_s = stays_matching_treat(dci_treat)
    infarct_s = stays_matching_diag(["cerebral infarction"])
    htn_s = stays_matching_diag(com.field("htn_eicu"))
    dm_s = stays_matching_diag(com.field("diabetes_eicu"))
    smk_s = stays_matching_diag(com.field("smoking_eicu"))
    clip_s = stays_matching_treat(["clip"])
    coil_s = stays_matching_treat(["coil"])

    hrt_s: set = set()
    if not admdrug.empty:
        m = admdrug["drugname"].apply(lambda s: any_pattern(s, hrt_pat))
        hrt_s = set(admdrug.loc[m, "patientunitstayid"].astype(int))

    if not apache.empty:
        apache["patientunitstayid"] = apache["patientunitstayid"].astype(int)
        apache["apachescore"] = pd.to_numeric(apache["apachescore"], errors="coerce")
        apache_map = apache.groupby("patientunitstayid")["apachescore"].max()
    else:
        apache_map = pd.Series(dtype="float64")

    def _col(name: str) -> pd.Series:
        return cohort.get(name, pd.Series("", index=cohort.index)).fillna("").astype(str)

    sid = cohort["patientunitstayid"].astype(int)
    disloc = _col("hospitaldischargelocation").str.lower()
    died = _col("hospitaldischargestatus").str.lower().eq("expired")
    offset = pd.to_numeric(cohort.get("unitdischargeoffset"), errors="coerce")

    def modality(row_id):
        if row_id in clip_s:
            return "clip"
        if row_id in coil_s:
            return "coil"
        return pd.NA

    out = pd.DataFrame({
        "source": "eicu",
        "patient_id": "eicu:" + cohort["uniquepid"].astype(str),
        "hospital_id": "eicu:" + _col("hospitalid"),
        "age": [_parse_age(a) for a in cohort["age"]],
        "sex": cohort["gender"].map({"Female": "F", "Male": "M"}),
        "vasospasm_dx": sid.isin(vaso_s).values,
        "dci_procedure": sid.isin(dci_s).values,
        "delayed_infarction": sid.isin(infarct_s).values,
        "hrt_exposure": sid.isin(hrt_s).values,
        "aneurysm_secured": sid.isin(clip_s | coil_s).values,
        "treatment_modality": [modality(i) for i in sid],
        "htn": sid.isin(htn_s).values,
        "smoking": sid.isin(smk_s).values,
        "diabetes": sid.isin(dm_s).values,
        "severity_gcs": pd.NA,
        "apache_score": sid.map(apache_map).astype("Float64").values,
        "died": died.values,
        "poor_disposition": (died | disloc.str.contains("hospice")).values,
        "icu_los_days": (offset / 1440.0).values,
    })
    return coerce_cohort(out)
