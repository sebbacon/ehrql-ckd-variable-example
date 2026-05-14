"""
Assurance tests for the CKD Stage 3+ phenotype (PH42).

Algorithm rules being tested:
  R1 — Most recent qualifying eGFR ≤ 60 ml/min/1.73m² before study end date
  R2 — At least one further qualifying eGFR ≤ 60 recorded STRICTLY MORE THAN
       90 days before the R1 measurement

Four strategies are tested:
  ckd_by_ckd_epi_2021    — eGFR computed from creatinine using 2021 CKD-EPI
  ckd_by_mdrd            — eGFR computed from creatinine using MDRD
  ckd_by_era_appropriate — MDRD pre-2014, CKD-EPI from 2014
  ckd_by_recorded_egfr   — eGFR values recorded directly in the EHR

SNOMED codes used:
  1000731000000107  Serum creatinine level
  1011481000000105  eGFR using creatinine CKD-EPI per 1.73m² (recorded eGFR)

Creatinine values (female, DOB 1970-01-01, measurements 2019-2020 → age ≈ 49-50):
  200 µmol/L → eGFR ≈ 25 ml/min under all formulas  (clearly below 60)
   70 µmol/L → eGFR ≈ 89 ml/min under all formulas  (clearly above 60)
"""

from datetime import date

from dataset_definition import dataset

SCR_CODE = "1000731000000107"   # Serum creatinine
EGFR_CODE = "1011481000000105"  # Recorded eGFR (CKD-EPI creatinine)

LOW_SCR_UMOL = 200.0   # eGFR ≈ 25 ml/min — below threshold
HIGH_SCR_UMOL = 70.0   # eGFR ≈ 89 ml/min — above threshold
LOW_EGFR = 35.0        # Recorded eGFR clearly below 60
HIGH_EGFR = 85.0       # Recorded eGFR clearly above 60

test_data = {
    # ----------------------------------------------------------------
    # Patient 1 — IN POPULATION (via creatinine strategies)
    # R1: most recent low eGFR from creatinine 200 µmol/L on 2020-01-01
    # R2: confirming creatinine 200 µmol/L on 2019-01-01 (365 days before)
    # No recorded eGFR → ckd_by_recorded_egfr is False
    # ----------------------------------------------------------------
    1: {
        "patients": {"date_of_birth": date(1970, 1, 1), "sex": "female"},
        "practice_registrations": [
            {"start_date": date(2000, 1, 1), "end_date": None}
        ],
        "clinical_events": [
            {"date": date(2019, 1, 1), "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
            {"date": date(2020, 1, 1), "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
        ],
        "expected_in_population": True,
        "expected_columns": {
            "sex": "female",
            "ckd_by_ckd_epi_2021": True,
            "ckd_by_mdrd": True,
            "ckd_by_era_appropriate": True,
            "ckd_by_recorded_egfr": False,
            "index_creatinine_date": date(2020, 1, 1),
            "index_creatinine_umol": LOW_SCR_UMOL,
            "recorded_egfr_index_date": None,
            "recorded_egfr_index_value": None,
        },
    },

    # ----------------------------------------------------------------
    # Patient 2 — EXCLUDED (R2 fails for all strategies: only one low reading)
    # ----------------------------------------------------------------
    2: {
        "patients": {"date_of_birth": date(1970, 1, 1), "sex": "female"},
        "practice_registrations": [
            {"start_date": date(2000, 1, 1), "end_date": None}
        ],
        "clinical_events": [
            {"date": date(2020, 1, 1), "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
        ],
        "expected_in_population": False,
    },

    # ----------------------------------------------------------------
    # Patient 3 — EXCLUDED (confirming exactly 90 days before — not strictly >90)
    # ----------------------------------------------------------------
    3: {
        "patients": {"date_of_birth": date(1970, 1, 1), "sex": "male"},
        "practice_registrations": [
            {"start_date": date(2000, 1, 1), "end_date": None}
        ],
        "clinical_events": [
            {"date": date(2019, 10, 3), "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
            {"date": date(2020, 1, 1),  "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
        ],
        "expected_in_population": False,
    },

    # ----------------------------------------------------------------
    # Patient 4 — IN POPULATION (boundary: confirming exactly 91 days before)
    # ----------------------------------------------------------------
    4: {
        "patients": {"date_of_birth": date(1970, 1, 1), "sex": "female"},
        "practice_registrations": [
            {"start_date": date(2000, 1, 1), "end_date": None}
        ],
        "clinical_events": [
            {"date": date(2019, 10, 2), "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
            {"date": date(2020, 1, 1),  "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
        ],
        "expected_in_population": True,
        "expected_columns": {
            "ckd_by_ckd_epi_2021": True,
            "ckd_by_mdrd": True,
            "ckd_by_era_appropriate": True,
            "ckd_by_recorded_egfr": False,
            "index_creatinine_date": date(2020, 1, 1),
            "index_creatinine_umol": LOW_SCR_UMOL,
        },
    },

    # ----------------------------------------------------------------
    # Patient 5 — EXCLUDED (latest creatinine gives eGFR > 60, older low reading ignored)
    # ----------------------------------------------------------------
    5: {
        "patients": {"date_of_birth": date(1970, 1, 1), "sex": "female"},
        "practice_registrations": [
            {"start_date": date(2000, 1, 1), "end_date": None}
        ],
        "clinical_events": [
            {"date": date(2018, 1, 1), "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
            {"date": date(2020, 1, 1), "snomedct_code": SCR_CODE, "numeric_value": HIGH_SCR_UMOL},
        ],
        "expected_in_population": False,
    },

    # ----------------------------------------------------------------
    # Patient 6 — EXCLUDED (no creatinine or eGFR recorded at all)
    # ----------------------------------------------------------------
    6: {
        "patients": {"date_of_birth": date(1970, 1, 1), "sex": "female"},
        "practice_registrations": [
            {"start_date": date(2000, 1, 1), "end_date": None}
        ],
        "clinical_events": [],
        "expected_in_population": False,
    },

    # ----------------------------------------------------------------
    # Patient 7 — EXCLUDED (not registered at study end date)
    # ----------------------------------------------------------------
    7: {
        "patients": {"date_of_birth": date(1970, 1, 1), "sex": "male"},
        "practice_registrations": [
            {"start_date": date(2000, 1, 1), "end_date": date(2019, 6, 1)}
        ],
        "clinical_events": [
            {"date": date(2019, 1, 1), "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
            {"date": date(2019, 5, 1), "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
        ],
        "expected_in_population": False,
    },

    # ----------------------------------------------------------------
    # Patient 8 — EXCLUDED (all readings after study end date)
    # ----------------------------------------------------------------
    8: {
        "patients": {"date_of_birth": date(1970, 1, 1), "sex": "female"},
        "practice_registrations": [
            {"start_date": date(2000, 1, 1), "end_date": None}
        ],
        "clinical_events": [
            {"date": date(2020, 4, 1), "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
            {"date": date(2020, 6, 1), "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
        ],
        "expected_in_population": False,
    },

    # ----------------------------------------------------------------
    # Patient 9 — EXCLUDED via creatinine strategies (sex is intersex — formula
    # returns None). No recorded eGFR either, so excluded entirely.
    # ----------------------------------------------------------------
    9: {
        "patients": {"date_of_birth": date(1970, 1, 1), "sex": "intersex"},
        "practice_registrations": [
            {"start_date": date(2000, 1, 1), "end_date": None}
        ],
        "clinical_events": [
            {"date": date(2019, 1, 1), "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
            {"date": date(2020, 1, 1), "snomedct_code": SCR_CODE, "numeric_value": LOW_SCR_UMOL},
        ],
        "expected_in_population": False,
    },

    # ----------------------------------------------------------------
    # Patient 10 — IN POPULATION via recorded eGFR only
    # No creatinine recorded, so all three creatinine strategies give False.
    # Two recorded eGFR values ≤ 60 more than 90 days apart → qualifies via
    # ckd_by_recorded_egfr. Demonstrates coverage of the creatinine gap.
    # ----------------------------------------------------------------
    10: {
        "patients": {"date_of_birth": date(1970, 1, 1), "sex": "intersex"},
        "practice_registrations": [
            {"start_date": date(2000, 1, 1), "end_date": None}
        ],
        "clinical_events": [
            {"date": date(2019, 1, 1), "snomedct_code": EGFR_CODE, "numeric_value": LOW_EGFR},
            {"date": date(2020, 1, 1), "snomedct_code": EGFR_CODE, "numeric_value": LOW_EGFR},
        ],
        "expected_in_population": True,
        "expected_columns": {
            "ckd_by_ckd_epi_2021": False,
            "ckd_by_mdrd": False,
            "ckd_by_era_appropriate": False,
            "ckd_by_recorded_egfr": True,
            "index_creatinine_date": None,
            "index_creatinine_umol": None,
            "recorded_egfr_index_date": date(2020, 1, 1),
            "recorded_egfr_index_value": LOW_EGFR,
        },
    },

    # ----------------------------------------------------------------
    # Patient 11 — IN POPULATION via both creatinine and recorded eGFR
    # Both strategies qualify — all four flags are True.
    # ----------------------------------------------------------------
    11: {
        "patients": {"date_of_birth": date(1970, 1, 1), "sex": "female"},
        "practice_registrations": [
            {"start_date": date(2000, 1, 1), "end_date": None}
        ],
        "clinical_events": [
            # Creatinine readings qualifying under all three formula strategies
            {"date": date(2019, 1, 1), "snomedct_code": SCR_CODE,  "numeric_value": LOW_SCR_UMOL},
            {"date": date(2020, 1, 1), "snomedct_code": SCR_CODE,  "numeric_value": LOW_SCR_UMOL},
            # Recorded eGFR readings also qualifying (>90 days apart)
            {"date": date(2019, 3, 1), "snomedct_code": EGFR_CODE, "numeric_value": LOW_EGFR},
            {"date": date(2020, 2, 1), "snomedct_code": EGFR_CODE, "numeric_value": LOW_EGFR},
        ],
        "expected_in_population": True,
        "expected_columns": {
            "ckd_by_ckd_epi_2021": True,
            "ckd_by_mdrd": True,
            "ckd_by_era_appropriate": True,
            "ckd_by_recorded_egfr": True,
        },
    },
}
