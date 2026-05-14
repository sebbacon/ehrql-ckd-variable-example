"""
Generate realistic dummy tables for the CKD Stage 3+ phenotype (PH42) study.

Run: uv run python scripts/generate_dummy_tables.py

Outputs to dummy-tables/:
  patients.csv, practice_registrations.csv, clinical_events.csv, medications.csv

Design notes
------------
Each patient is assigned a profile before any events are generated, so that the
mapping between input fractions and output classification is transparent.

CKD patients are split across four profiles that exercise all four identification
strategies in the dataset definition:

  "creatinine_and_egfr"  — qualifying creatinine + qualifying recorded eGFR
                           → all four flags True
  "creatinine_only"      — qualifying creatinine, no recorded eGFR events
                           → ckd_by_recorded_egfr False
  "recorded_egfr_only"   — recorded eGFR ≤ 60, no creatinine events
                           → only ckd_by_recorded_egfr True
  "borderline_scr"       — creatinine in the zone where CKD-EPI and MDRD
                           disagree (~60 boundary); no recorded eGFR
                           → some creatinine flags True, others False

Creatinine values are computed from a target eGFR using the inverse 2021
CKD-EPI formula so values are physiologically plausible for the patient's
age and sex.

Events are generated across the patient's full history back to 2005, not just
the measurement window, so that time-series patterns are realistic.
"""

import csv
import os
from datetime import date, timedelta

import numpy as np

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
RNG = np.random.default_rng(42)

# ---------------------------------------------------------------------------
# Population parameters
# ---------------------------------------------------------------------------
N_PATIENTS = 5_000
STUDY_END = date(2020, 3, 31)
HISTORY_START = date(2005, 1, 1)  # How far back to generate events

# ---------------------------------------------------------------------------
# Prevalence
# ---------------------------------------------------------------------------
# CKD Stage 3+ in a UK GP-registered adult population (QOF register rates)
CKD_PREVALENCE = 0.06  # 6% → ~300 patients

# Of non-CKD patients, proportion with at least one creatinine measurement
# (routine blood tests, diabetes monitoring, hypertension follow-up, etc.)
FRAC_NON_CKD_WITH_SCR = 0.35

# ---------------------------------------------------------------------------
# CKD profile weights (must sum to 1.0)
# ---------------------------------------------------------------------------
# "creatinine_and_egfr": qualifying creatinine readings + recorded eGFR
# "creatinine_only":     qualifying creatinine, no recorded eGFR code
# "recorded_egfr_only":  recorded eGFR ≤ 60, no creatinine code
# "borderline_scr":      creatinine near threshold (CKD-EPI vs MDRD may differ)
CKD_PROFILE_WEIGHTS = {
    "creatinine_and_egfr": 0.50,
    "creatinine_only":     0.30,
    "recorded_egfr_only":  0.10,
    "borderline_scr":      0.10,
}

# ---------------------------------------------------------------------------
# Measurement frequency (events per year)
# ---------------------------------------------------------------------------
CKD_SCR_PER_YEAR    = 2.0   # CKD patients monitored ~twice yearly
CKD_EGFR_PER_YEAR   = 1.5   # Recorded eGFR often accompanies creatinine result
NON_CKD_SCR_PER_YEAR = 0.25  # Occasional creatinine for other reasons
NON_CKD_EGFR_PER_YEAR = 0.10

# ---------------------------------------------------------------------------
# CKD eGFR target ranges (ml/min/1.73m²) and creatinine noise
# ---------------------------------------------------------------------------
CKD_EGFR_RANGE      = (15, 59)    # Stage 3 and below
BORDERLINE_EGFR_RANGE = (55, 65)  # Straddles the 60 threshold
NORMAL_EGFR_RANGE   = (65, 110)   # Non-CKD patients
SCR_NOISE_FRAC      = 0.08        # ±8% biological variability on creatinine

# ---------------------------------------------------------------------------
# Normal creatinine ranges for non-CKD patients (µmol/L)
# Used when we don't back-compute from eGFR
# ---------------------------------------------------------------------------
NORMAL_SCR_FEMALE_RANGE = (50, 90)
NORMAL_SCR_MALE_RANGE   = (60, 110)

# ---------------------------------------------------------------------------
# SNOMED codes
# ---------------------------------------------------------------------------
SCR_CODE  = "1000731000000107"  # Serum creatinine level
EGFR_CODE = "1011481000000105"  # eGFR using creatinine CKD-EPI per 1.73m²

# ---------------------------------------------------------------------------
# Age distribution (UK GP-registered adult population, weighted toward older)
# CKD Stage 3+ prevalence rises steeply with age so we oversample older patients
# ---------------------------------------------------------------------------
AGE_MEAN_CKD     = 72   # CKD patients are predominantly elderly
AGE_SD_CKD       = 10
AGE_MEAN_NON_CKD = 50
AGE_SD_NON_CKD   = 18
AGE_MIN          = 18
AGE_MAX          = 95

# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "dummy-tables")

UMOL_TO_MGDL = 88.42


# ---------------------------------------------------------------------------
# Formulae
# ---------------------------------------------------------------------------

def ckd_epi_2021(scr_mgdl: float, sex: str, age: float) -> float | None:
    """2021 CKD-EPI (race-free). Returns None for unknown sex."""
    if sex == "female":
        kappa, alpha = 0.7, (-0.241 if scr_mgdl <= 0.7 else -1.200)
        return 142 * (scr_mgdl / 0.7) ** alpha * 0.9938 ** age
    if sex == "male":
        kappa, alpha = 0.9, (-0.302 if scr_mgdl <= 0.9 else -1.200)
        return 142 * (scr_mgdl / 0.9) ** alpha * 0.9938 ** age
    return None


def scr_from_egfr_ckd_epi(egfr: float, sex: str, age: float) -> float:
    """
    Inverse 2021 CKD-EPI (>kappa branch, appropriate for CKD patients).
    Returns creatinine in mg/dL.
    """
    kappa = 0.7 if sex == "female" else 0.9
    age_factor = 0.9938 ** age
    return kappa * (egfr / (142 * age_factor)) ** (-1 / 1.2)


def random_date(start: date, end: date) -> date:
    delta = (end - start).days
    return start + timedelta(days=int(RNG.integers(0, delta + 1)))


def dates_across_years(start: date, end: date, rate_per_year: float) -> list[date]:
    """Generate Poisson-distributed event dates between start and end."""
    n_years = max((end - start).days / 365.25, 0)
    n_events = RNG.poisson(rate_per_year * n_years)
    return sorted(random_date(start, end) for _ in range(n_events))


def clamp_age(age: float) -> int:
    return max(AGE_MIN, min(AGE_MAX, int(round(age))))


# ---------------------------------------------------------------------------
# Build patient profiles
# ---------------------------------------------------------------------------

def build_patients() -> list[dict]:
    patients = []
    n_ckd = int(round(N_PATIENTS * CKD_PREVALENCE))
    n_non_ckd = N_PATIENTS - n_ckd

    sexes = ["female", "male"]  # intersex excluded from formula; rare in population
    ckd_profiles = list(CKD_PROFILE_WEIGHTS.keys())
    ckd_profile_probs = list(CKD_PROFILE_WEIGHTS.values())

    pid = 1

    # ---- CKD patients ----
    for _ in range(n_ckd):
        sex = RNG.choice(sexes)
        age_at_study_end = clamp_age(RNG.normal(AGE_MEAN_CKD, AGE_SD_CKD))
        dob = STUDY_END - timedelta(days=int(age_at_study_end * 365.25))
        dob = dob.replace(day=1)  # dates of birth rounded to first of month

        profile = RNG.choice(ckd_profiles, p=ckd_profile_probs)

        if profile == "borderline_scr":
            target_egfr = float(RNG.uniform(*BORDERLINE_EGFR_RANGE))
        else:
            target_egfr = float(RNG.uniform(*CKD_EGFR_RANGE))

        patients.append({
            "patient_id": pid,
            "sex": sex,
            "date_of_birth": dob,
            "age_at_study_end": age_at_study_end,
            "ckd": True,
            "profile": profile,
            "target_egfr": target_egfr,
        })
        pid += 1

    # ---- Non-CKD patients ----
    for _ in range(n_non_ckd):
        sex = RNG.choice(sexes)
        age_at_study_end = clamp_age(RNG.normal(AGE_MEAN_NON_CKD, AGE_SD_NON_CKD))
        dob = STUDY_END - timedelta(days=int(age_at_study_end * 365.25))
        dob = dob.replace(day=1)

        patients.append({
            "patient_id": pid,
            "sex": sex,
            "date_of_birth": dob,
            "age_at_study_end": age_at_study_end,
            "ckd": False,
            "profile": "non_ckd",
            "target_egfr": float(RNG.uniform(*NORMAL_EGFR_RANGE)),
        })
        pid += 1

    return patients


# ---------------------------------------------------------------------------
# Generate events for a patient
# ---------------------------------------------------------------------------

def make_scr_event(pid: int, event_date: date, scr_umol: float) -> dict:
    return {
        "patient_id": pid,
        "date": event_date.isoformat(),
        "snomedct_code": SCR_CODE,
        "numeric_value": round(scr_umol, 1),
    }


def make_egfr_event(pid: int, event_date: date, egfr: float) -> dict:
    return {
        "patient_id": pid,
        "date": event_date.isoformat(),
        "snomedct_code": EGFR_CODE,
        "numeric_value": round(egfr, 1),
    }


def generate_events(p: dict) -> list[dict]:
    pid = p["patient_id"]
    sex = p["sex"]
    dob = p["date_of_birth"]
    profile = p["profile"]
    target_egfr = p["target_egfr"]
    events = []

    if profile == "non_ckd":
        if RNG.random() > FRAC_NON_CKD_WITH_SCR:
            return []
        for d in dates_across_years(HISTORY_START, STUDY_END, NON_CKD_SCR_PER_YEAR):
            age = (d - dob).days / 365.25
            # Non-CKD: random normal creatinine, not back-computed
            if sex == "female":
                scr_umol = float(RNG.uniform(*NORMAL_SCR_FEMALE_RANGE))
            else:
                scr_umol = float(RNG.uniform(*NORMAL_SCR_MALE_RANGE))
            events.append(make_scr_event(pid, d, scr_umol))
        return events

    # ---- CKD patient ----
    # Determine whether this patient has creatinine events and/or recorded eGFR
    has_scr  = profile in ("creatinine_and_egfr", "creatinine_only",  "borderline_scr")
    has_egfr = profile in ("creatinine_and_egfr", "recorded_egfr_only")

    if has_scr:
        for d in dates_across_years(HISTORY_START, STUDY_END, CKD_SCR_PER_YEAR):
            age = max((d - dob).days / 365.25, 1.0)
            # Compute creatinine from target eGFR, add biological noise
            scr_mgdl = scr_from_egfr_ckd_epi(target_egfr, sex, age)
            noise = float(RNG.normal(1.0, SCR_NOISE_FRAC))
            scr_umol = scr_mgdl * noise * UMOL_TO_MGDL
            events.append(make_scr_event(pid, d, max(scr_umol, 10.0)))

    if has_egfr:
        for d in dates_across_years(HISTORY_START, STUDY_END, CKD_EGFR_PER_YEAR):
            age = max((d - dob).days / 365.25, 1.0)
            # Recorded eGFR: target eGFR with small noise
            noise = float(RNG.normal(1.0, 0.05))
            egfr = target_egfr * noise
            events.append(make_egfr_event(pid, d, max(egfr, 1.0)))

    return events


# ---------------------------------------------------------------------------
# Write tables
# ---------------------------------------------------------------------------

def write_patients(patients: list[dict]) -> None:
    path = os.path.join(OUT_DIR, "patients.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["patient_id", "date_of_birth", "sex", "date_of_death"])
        w.writeheader()
        for p in patients:
            w.writerow({
                "patient_id": p["patient_id"],
                "date_of_birth": p["date_of_birth"].isoformat(),
                "sex": p["sex"],
                "date_of_death": "",
            })


def write_registrations(patients: list[dict]) -> None:
    path = os.path.join(OUT_DIR, "practice_registrations.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["patient_id", "start_date", "end_date", "practice_pseudo_id"])
        w.writeheader()
        for p in patients:
            w.writerow({
                "patient_id": p["patient_id"],
                "start_date": "2000-01-01",
                "end_date": "",
                "practice_pseudo_id": 1,
            })


def write_clinical_events(all_events: list[dict]) -> None:
    path = os.path.join(OUT_DIR, "clinical_events.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["patient_id", "date", "snomedct_code", "numeric_value"])
        w.writeheader()
        for e in sorted(all_events, key=lambda x: (x["patient_id"], x["date"])):
            w.writerow(e)


def write_medications(patients: list[dict]) -> None:
    # No medications relevant to this phenotype; write empty table with headers only
    path = os.path.join(OUT_DIR, "medications.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["patient_id", "date", "dmd_code"])
        w.writeheader()


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(patients: list[dict], all_events: list[dict]) -> None:
    ckd = [p for p in patients if p["ckd"]]
    profiles = {}
    for p in ckd:
        profiles[p["profile"]] = profiles.get(p["profile"], 0) + 1

    scr_pids  = {e["patient_id"] for e in all_events if e["snomedct_code"] == SCR_CODE}
    egfr_pids = {e["patient_id"] for e in all_events if e["snomedct_code"] == EGFR_CODE}

    print(f"\n{'─' * 50}")
    print(f"  Dummy data summary")
    print(f"{'─' * 50}")
    print(f"  Total patients:            {len(patients):>6,}")
    print(f"  CKD Stage 3+ patients:     {len(ckd):>6,}  ({100 * len(ckd) / len(patients):.1f}%)")
    print(f"  Non-CKD patients:          {len(patients) - len(ckd):>6,}")
    print()
    print(f"  CKD profiles:")
    for profile, count in sorted(profiles.items()):
        expected_flag = {
            "creatinine_and_egfr": "all four flags True",
            "creatinine_only":     "ckd_by_recorded_egfr False",
            "recorded_egfr_only":  "only ckd_by_recorded_egfr True",
            "borderline_scr":      "some creatinine flags may differ",
        }.get(profile, "")
        print(f"    {profile:<25} {count:>5}  ({expected_flag})")
    print()
    print(f"  Patients with creatinine events:    {len(scr_pids):>5,}")
    print(f"  Patients with recorded eGFR events: {len(egfr_pids):>5,}")
    print(f"  Total clinical events:              {len(all_events):>5,}")
    print(f"{'─' * 50}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Building patient profiles...")
    patients = build_patients()

    print("Generating clinical events...")
    all_events = []
    for p in patients:
        all_events.extend(generate_events(p))

    print("Writing dummy tables...")
    write_patients(patients)
    write_registrations(patients)
    write_clinical_events(all_events)
    write_medications(patients)

    print_summary(patients, all_events)
    print("Done. Dummy tables written to dummy-tables/")
