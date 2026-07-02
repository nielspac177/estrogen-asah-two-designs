"""Deterministic synthetic data — the reproducibility keystone.

No PhysioNet data ships with this repo. Instead a small, fully specified cohort
(`SPECS`) is rendered into MIMIC-IV-shaped and eICU-shaped source tables so the
extractors and the whole pipeline can run — in tests and in CI — with zero real
data. `SPECS` is the single source of ground truth: extraction tests assert that
the extractors reproduce these flags from the encoded tables.

Every branch is covered: premenopausal / postmenopausal women, men, DCI via each
composite component, HRT +/-, ICD-9 and ICD-10 encodings, and non-cohort
patients (unruptured aneurysm, traumatic SAH, unrelated admission).
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

# --- ground truth: one dict per patient-admission ---------------------------
# in_cohort := is_sah AND (has_aneurysm OR secured); coding column picks ICD-9/10.
SPECS: list[dict] = [
    # id  sex age  sah   aneu  secured  coding  vaso   dciproc infarct hrt  htn  smk  dm   died  hospice los  in_cohort
    dict(id=1,  sex="F", age=44, is_sah=1, has_aneurysm=1, secured="clip", coding="i10", vasospasm=1, dci_proc=0, infarct=0, hrt=0, htn=1, smoking=1, diabetes=0, died=0, hospice=0, los=12.0, in_cohort=1),
    dict(id=2,  sex="F", age=58, is_sah=1, has_aneurysm=0, secured="coil", coding="i10", vasospasm=0, dci_proc=1, infarct=0, hrt=1, htn=1, smoking=0, diabetes=1, died=0, hospice=0, los=18.0, in_cohort=1),
    dict(id=3,  sex="F", age=62, is_sah=1, has_aneurysm=1, secured=None,   coding="i9",  vasospasm=0, dci_proc=0, infarct=0, hrt=0, htn=1, smoking=0, diabetes=0, died=1, hospice=0, los=4.0,  in_cohort=1),
    dict(id=4,  sex="F", age=49, is_sah=1, has_aneurysm=1, secured="coil", coding="i10", vasospasm=0, dci_proc=0, infarct=1, hrt=0, htn=0, smoking=1, diabetes=0, died=0, hospice=1, los=22.0, in_cohort=1),
    dict(id=5,  sex="M", age=55, is_sah=1, has_aneurysm=1, secured="clip", coding="i10", vasospasm=1, dci_proc=0, infarct=0, hrt=0, htn=1, smoking=1, diabetes=0, died=0, hospice=0, los=15.0, in_cohort=1),
    dict(id=6,  sex="M", age=70, is_sah=1, has_aneurysm=0, secured="coil", coding="i9",  vasospasm=0, dci_proc=0, infarct=0, hrt=0, htn=1, smoking=0, diabetes=1, died=0, hospice=0, los=9.0,  in_cohort=1),
    dict(id=7,  sex="F", age=51, is_sah=1, has_aneurysm=1, secured=None,   coding="i10", vasospasm=0, dci_proc=0, infarct=0, hrt=0, htn=0, smoking=0, diabetes=0, died=0, hospice=0, los=7.0,  in_cohort=1),
    dict(id=8,  sex="F", age=50, is_sah=1, has_aneurysm=1, secured="clip", coding="i10", vasospasm=1, dci_proc=1, infarct=0, hrt=1, htn=1, smoking=1, diabetes=0, died=0, hospice=0, los=20.0, in_cohort=1),
    dict(id=9,  sex="F", age=60, is_sah=0, has_aneurysm=1, secured=None,   coding="i10", vasospasm=0, dci_proc=0, infarct=0, hrt=0, htn=1, smoking=0, diabetes=0, died=0, hospice=0, los=3.0,  in_cohort=0),  # unruptured aneurysm, no SAH
    dict(id=10, sex="M", age=65, is_sah=1, has_aneurysm=0, secured=None,   coding="trauma", vasospasm=0, dci_proc=0, infarct=0, hrt=0, htn=0, smoking=1, diabetes=0, died=0, hospice=0, los=5.0, in_cohort=0),  # traumatic SAH
    dict(id=11, sex="F", age=72, is_sah=0, has_aneurysm=0, secured=None,   coding="other", vasospasm=0, dci_proc=0, infarct=0, hrt=1, htn=1, smoking=0, diabetes=1, died=0, hospice=0, los=6.0, in_cohort=0),  # unrelated (pneumonia)
    dict(id=12, sex="F", age=40, is_sah=1, has_aneurysm=1, secured="coil", coding="i10", vasospasm=0, dci_proc=0, infarct=0, hrt=1, htn=0, smoking=0, diabetes=0, died=0, hospice=0, los=11.0, in_cohort=1),
]

# derived ground-truth sets (single source of truth for tests)
COHORT_IDS = [s["id"] for s in SPECS if s["in_cohort"]]
DCI_IDS = [s["id"] for s in SPECS if s["in_cohort"] and (s["vasospasm"] or s["dci_proc"] or s["infarct"])]


# --- ICD code choices per phenotype -----------------------------------------
_SAH = {"i10": ("I609", 10), "i9": ("430", 9)}
_ANEURYSM = {"i10": ("I671", 10), "i9": ("4373", 9)}
_VASOSPASM = ("I67848", 10)
_INFARCT = {"i10": ("I639", 10), "i9": ("434", 9)}
_CLIP = {"i10": ("03VG0CZ", 10), "i9": ("3951", 9)}
_COIL = {"i10": ("03LG3DZ", 10), "i9": ("3975", 9)}
_DCI_PROC = {"i10": ("037G3ZZ", 10), "i9": ("0061", 9)}
_HTN = ("I10", 10)
_SMOKING = ("F17210", 10)
_DIABETES = ("E119", 10)


def _title(code: str) -> str:
    m = {
        "I609": "Nontraumatic subarachnoid hemorrhage, unspecified",
        "430": "Subarachnoid hemorrhage",
        "I671": "Cerebral aneurysm, nonruptured",
        "4373": "Cerebral aneurysm, nonruptured",
        "I67848": "Other cerebrovascular vasospasm and vasoconstriction",
        "I639": "Cerebral infarction, unspecified",
        "434": "Occlusion of cerebral arteries",
        "I10": "Essential (primary) hypertension",
        "F17210": "Nicotine dependence, cigarettes, uncomplicated",
        "E119": "Type 2 diabetes mellitus without complications",
    }
    return m.get(code, code)


# ---------------------------------------------------------------------------
# MIMIC-IV-shaped tables
# ---------------------------------------------------------------------------
def mimic_tables() -> dict[str, pd.DataFrame]:
    patients, admissions, dx, proc, rx = [], [], [], [], []
    d_icd: dict[tuple, str] = {}
    base = pd.Timestamp("2150-01-01")
    for s in SPECS:
        sid = s["id"]
        hadm = 100000 + sid
        patients.append(dict(subject_id=sid, gender=s["sex"], anchor_age=s["age"]))
        admit = base + pd.Timedelta(days=sid)
        disch = admit + pd.Timedelta(days=s["los"])
        admissions.append(dict(
            subject_id=sid, hadm_id=hadm, admittime=admit, dischtime=disch,
            hospital_expire_flag=int(s["died"]),
            discharge_location=("HOSPICE" if s["hospice"] else ("DIED" if s["died"] else "HOME")),
        ))

        def add_dx(code_ver, sid=sid, hadm=hadm):
            code, ver = code_ver
            dx.append(dict(subject_id=sid, hadm_id=hadm, icd_code=code, icd_version=ver))
            d_icd[(code, ver)] = _title(code)

        def add_proc(code_ver, sid=sid, hadm=hadm):
            code, ver = code_ver
            proc.append(dict(subject_id=sid, hadm_id=hadm, icd_code=code, icd_version=ver))

        coding = s["coding"]
        if s["is_sah"]:
            if coding == "trauma":
                # traumatic SAH: encode as head-injury SAH to be excluded upstream
                dx.append(dict(subject_id=sid, hadm_id=hadm, icd_code="S066X9A", icd_version=10))
                d_icd[("S066X9A", 10)] = "Traumatic subarachnoid hemorrhage"
            else:
                add_dx(_SAH[coding])
        if s["has_aneurysm"]:
            add_dx(_ANEURYSM[coding if coding in ("i9", "i10") else "i10"])
        if s["secured"] == "clip":
            add_proc(_CLIP[coding if coding in ("i9", "i10") else "i10"])
        elif s["secured"] == "coil":
            add_proc(_COIL[coding if coding in ("i9", "i10") else "i10"])
        if s["vasospasm"]:
            add_dx(_VASOSPASM)
        if s["dci_proc"]:
            add_proc(_DCI_PROC[coding if coding in ("i9", "i10") else "i10"])
        if s["infarct"]:
            add_dx(_INFARCT[coding if coding in ("i9", "i10") else "i10"])
        if s["htn"]:
            add_dx(_HTN)
        if s["smoking"]:
            add_dx(_SMOKING)
        if s["diabetes"]:
            add_dx(_DIABETES)
        if coding == "other":
            add_dx(("J189", 10))  # pneumonia, unrelated
            d_icd[("J189", 10)] = "Pneumonia, unspecified organism"
        if s["hrt"]:
            rx.append(dict(subject_id=sid, hadm_id=hadm, drug="Estradiol",
                           starttime=admit + pd.Timedelta(hours=6)))

    d_icd_rows = [dict(icd_code=c, icd_version=v, long_title=t) for (c, v), t in d_icd.items()]
    return {
        "patients": pd.DataFrame(patients),
        "admissions": pd.DataFrame(admissions),
        "diagnoses_icd": pd.DataFrame(dx),
        "procedures_icd": pd.DataFrame(proc),
        "prescriptions": pd.DataFrame(rx),
        "d_icd_diagnoses": pd.DataFrame(d_icd_rows),
    }


def load_mimic_duckdb(con=None):
    """Create schema ``mimiciv_hosp`` with synthetic tables in a DuckDB connection."""
    import duckdb

    if con is None:
        con = duckdb.connect()
    con.execute("CREATE SCHEMA IF NOT EXISTS mimiciv_hosp")
    for name, df in mimic_tables().items():
        con.register(f"_syn_{name}", df)
        con.execute(f"CREATE OR REPLACE TABLE mimiciv_hosp.{name} AS SELECT * FROM _syn_{name}")
        con.unregister(f"_syn_{name}")
    return con


# ---------------------------------------------------------------------------
# eICU-shaped tables
# ---------------------------------------------------------------------------
def eicu_tables() -> dict[str, pd.DataFrame]:
    patient, diagnosis, admitdx, treatment, admdrug, apache = [], [], [], [], [], []
    for s in SPECS:
        uid = f"uid-{s['id']}"
        stay = 2000 + s["id"]
        age = ">89" if s["age"] > 89 else str(s["age"])
        patient.append(dict(
            patientunitstayid=stay, uniquepid=uid, gender=("Female" if s["sex"] == "F" else "Male"),
            age=age, hospitalid=10 + (s["id"] % 3),
            unitdischargestatus=("Expired" if s["died"] else "Alive"),
            hospitaldischargestatus=("Expired" if s["died"] else "Alive"),
            hospitaldischargelocation=("Hospice" if s["hospice"] else ("Death" if s["died"] else "Home")),
            unitdischargeoffset=int(s["los"] * 24 * 60),
        ))
        apache.append(dict(patientunitstayid=stay, apachescore=40 + s["id"]))
        # cohort encoding
        if s["is_sah"]:
            if s["coding"] == "trauma":
                diagnosis.append(dict(patientunitstayid=stay,
                    diagnosisstring="neurologic|trauma - CNS|intracranial injury|with subarachnoid hemorrhage",
                    icd9code=""))
            else:
                base = "neurologic|disorders of vasculature|stroke|hemorrhagic stroke|subarachnoid hemorrhage"
                if s["has_aneurysm"]:
                    base += "|from ruptured berry aneurysm"
                diagnosis.append(dict(patientunitstayid=stay, diagnosisstring=base, icd9code="430"))
                name = "Subarachnoid hemorrhage/intracranial aneurysm"
                if s["secured"]:
                    name += ", surgery for"
                admitdx.append(dict(patientunitstayid=stay, admitdxname=name))
        elif s["coding"] == "other":
            diagnosis.append(dict(patientunitstayid=stay,
                diagnosisstring="pulmonary|disorders of the airways|pneumonia", icd9code="486"))
        elif s["has_aneurysm"]:
            diagnosis.append(dict(patientunitstayid=stay,
                diagnosisstring="neurologic|disorders of vasculature|cerebral aneurysm|unruptured", icd9code="4373"))
        # DCI components
        if s["vasospasm"]:
            diagnosis.append(dict(patientunitstayid=stay,
                diagnosisstring="neurologic|disorders of vasculature|stroke|hemorrhagic stroke|subarachnoid hemorrhage|with vasospasm",
                icd9code=""))
        if s["dci_proc"]:
            treatment.append(dict(patientunitstayid=stay,
                treatmentstring="neurologic|ICH/ cerebral infarct|angiogram|with cerebral angioplasty"))
        if s["infarct"]:
            diagnosis.append(dict(patientunitstayid=stay,
                diagnosisstring="neurologic|disorders of vasculature|stroke|cerebral infarction", icd9code="434"))
        # covariates
        for present, dstr, icd in [
            (s["htn"], "cardiovascular|hypertension|essential hypertension", "401.9"),
            (s["smoking"], "social history|smoking|current smoker", ""),
            (s["diabetes"], "endocrine|glucose metabolism|diabetes mellitus", "250"),
        ]:
            if present:
                diagnosis.append(dict(patientunitstayid=stay, diagnosisstring=dstr, icd9code=icd))
        if s["secured"] == "clip":
            treatment.append(dict(patientunitstayid=stay,
                treatmentstring="surgery|neurosurgery|aneurysm clipping"))
        elif s["secured"] == "coil":
            treatment.append(dict(patientunitstayid=stay,
                treatmentstring="neurologic|procedures|aneurysm coiling"))
        if s["hrt"]:
            admdrug.append(dict(patientunitstayid=stay, drugname="ESTRADIOL 1 MG PO TABS"))

    return {
        "patient": pd.DataFrame(patient),
        "diagnosis": pd.DataFrame(diagnosis),
        "admissionDx": pd.DataFrame(admitdx),
        "treatment": pd.DataFrame(treatment),
        "admissionDrug": pd.DataFrame(admdrug),
        "apachePatientResult": pd.DataFrame(apache),
    }


def simulate_cohort(n: int = 2000, seed: int = 0, log_or_post: float | None = None) -> pd.DataFrame:
    """A large *harmonized* cohort with a known injected effect, for estimator validation.

    Unlike the SPECS fixtures (which test extraction), this produces a ready-to-analyze
    cohort where postmenopausal vs premenopausal DCI odds follow ``exp(log_or_post)``,
    so analysis code can be checked for recovering a known truth.
    """
    import numpy as np

    from .harmonize.common_schema import coerce_cohort

    if log_or_post is None:
        log_or_post = float(np.log(0.6))  # protective, matching the hypothesis direction
    rng = np.random.default_rng(seed)
    source = rng.choice(["mimic_iv", "eicu"], n)
    hosp = np.where(source == "mimic_iv", "mimic_iv",
                    np.array(["eicu:" + str(h) for h in rng.integers(1, 6, n)]))
    sex = rng.choice(["F", "M"], n, p=[0.6, 0.4])
    age = np.clip(rng.normal(58, 12, n), 18, 95)
    htn = rng.random(n) < 0.5
    smoking = rng.random(n) < 0.3
    diabetes = rng.random(n) < 0.2
    post = (sex == "F") & (age >= 51)
    male = sex == "M"
    # No independent age term: menopause is age-defined, and the real analysis does
    # not (cannot) co-adjust age, so the injected OR must be recoverable unadjusted.
    lp = -1.0 + log_or_post * post + 0.3 * male + 0.4 * htn + 0.5 * smoking
    p = 1.0 / (1.0 + np.exp(-lp))
    dci = rng.random(n) < p
    died = rng.random(n) < 0.15
    out = pd.DataFrame({
        "source": source,
        "patient_id": [f"{s}:{i}" for i, s in enumerate(source)],
        "hospital_id": hosp,
        "age": age,
        "sex": sex,
        "vasospasm_dx": dci,          # route all DCI through one component -> composite==dci
        "dci_procedure": False,
        "delayed_infarction": False,
        "hrt_exposure": rng.random(n) < 0.05,
        "aneurysm_secured": rng.random(n) < 0.8,
        "treatment_modality": rng.choice(["clip", "coil"], n),
        "htn": htn,
        "smoking": smoking,
        "diabetes": diabetes,
        "severity_gcs": pd.NA,
        "apache_score": np.where(source == "eicu", rng.normal(60, 15, n), np.nan),
        "died": died,
        "poor_disposition": died | (rng.random(n) < 0.05),
        "icu_los_days": np.clip(rng.normal(12, 6, n), 0, None),
    })
    return coerce_cohort(out)


def write_eicu(dirpath: str | Path) -> Path:
    """Write synthetic eICU tables as gzipped CSVs (mirrors the real v2.0 layout)."""
    dirpath = Path(dirpath)
    dirpath.mkdir(parents=True, exist_ok=True)
    for name, df in eicu_tables().items():
        with gzip.open(dirpath / f"{name}.csv.gz", "wt", newline="", encoding="utf-8") as fh:
            df.to_csv(fh, index=False)
    return dirpath
