"""
CKD Stage 3+ phenotype (PH42) — modified CALIBER/eGFR-CKD-EPI algorithm.
Source: https://phenotypes.healthdatagateway.org/phenotypes/PH42/version/2933/detail/

Four CKD identification strategies are provided so researchers can choose:

  ckd_by_ckd_epi_2021    — eGFR computed from creatinine using 2021 CKD-EPI
                           (race-free) applied throughout
  ckd_by_mdrd            — eGFR computed from creatinine using 4-variable MDRD
                           applied throughout
  ckd_by_era_appropriate — creatinine-derived eGFR using MDRD before 2014-01-01
                           and CKD-EPI 2021 from 2014-01-01 onward; closest to
                           the original PH42 algorithm intent
  ckd_by_recorded_egfr   — eGFR values recorded directly in the EHR (no formula
                           applied); reflects what the clinician actually saw and
                           covers patients with a recorded eGFR but no creatinine

A patient qualifies under a given strategy if:
  1. Their most recent qualifying eGFR ≤ 60 ml/min/1.73m² before the study end
     date (the "index measurement").
  2. They also have at least one qualifying eGFR ≤ 60 recorded strictly more
     than 90 days before that index measurement.

The population includes anyone who qualifies under at least one strategy.

Creatinine formula notes:
  - Creatinine is recorded in µmol/L in UK EHRs; formulas use mg/dL (÷ 88.42).
  - Patients with sex other than male/female are excluded by the creatinine
    strategies (both formulas are sex-specific); they may still qualify via
    ckd_by_recorded_egfr if eGFR values were recorded directly.
  - MDRD uses 175 as the leading coefficient (IDMS-standardised creatinine,
    standard in UK labs since ~2006) with no race adjustment.
  - The 2014 formula cutoff matches NICE CG182, which recommended UK labs switch
    from MDRD to CKD-EPI from that year.
"""

from ehrql import case, codelist_from_csv, create_dataset, days, when
from ehrql.tables.tpp import clinical_events, patients, practice_registrations

# ============================================================
# WARNING: All codelists below must be carefully reviewed.
# The wrong codelist may have been selected, or a codelist
# may be out of date. Do not use this definition in
# production without expert clinical review of every codelist.
# ============================================================

# CONFIDENCE: HIGH — UK Renal Registry, 2022, purpose-built for TPP (subset
#   where all codes carry attached numeric values). Covers serum and plasma
#   creatinine variants.
# Other candidates: user/richard-croker/creatinine-tests/3af470de (broader,
#   not confirmed to have numeric values in TPP for all codes)
creatinine_codelist = codelist_from_csv(
    "codelists/ukrr-creatinine-tests-level.csv", column="code"
)

# CONFIDENCE: HIGH — UK Renal Registry, 2022, purpose-built for TPP (subset
#   where all codes carry attached numeric values). Covers CKD-EPI creatinine,
#   CKD-EPI cystatin C, MDRD, and generic GFR codes.
# Other candidates: pincer/egfr/v1.8 (PRIMIS, 2021, 8 codes; adds two codes not
#   confirmed to carry numeric values in TPP)
egfr_codelist = codelist_from_csv(
    "codelists/ukrr-egfr-tests-level.csv", column="code"
)

dataset = create_dataset()

study_end_date = "2020-03-31"
CKD_EGFR_THRESHOLD = 60
PERSISTENCE_GAP = days(90)

# Creatinine is stored in µmol/L; both formulas require mg/dL
UMOL_TO_MGDL = 88.42

# Formula cutoff: MDRD before this date, CKD-EPI 2021 from this date onward
FORMULA_CUTOFF_DATE = "2014-01-01"

# Sex-specific creatinine kappa values for CKD-EPI 2021
FEMALE_KAPPA = 0.7
MALE_KAPPA = 0.9


def ckd_epi_2021(scr_mgdl, sex, age_years):
    """
    2021 CKD-EPI equation (race-free).
    Returns eGFR in ml/min/1.73m², or None for non-male/female sex.
    Reference: Inker et al., NEJM 2021, doi:10.1056/NEJMoa2102953
    """
    return case(
        when((sex == "female") & (scr_mgdl <= FEMALE_KAPPA)).then(
            142 * (scr_mgdl / FEMALE_KAPPA) ** -0.241 * 0.9938 ** age_years
        ),
        when((sex == "female") & (scr_mgdl > FEMALE_KAPPA)).then(
            142 * (scr_mgdl / FEMALE_KAPPA) ** -1.200 * 0.9938 ** age_years
        ),
        when((sex == "male") & (scr_mgdl <= MALE_KAPPA)).then(
            142 * (scr_mgdl / MALE_KAPPA) ** -0.302 * 0.9938 ** age_years
        ),
        when((sex == "male") & (scr_mgdl > MALE_KAPPA)).then(
            142 * (scr_mgdl / MALE_KAPPA) ** -1.200 * 0.9938 ** age_years
        ),
        otherwise=None,
    )


def mdrd_4var(scr_mgdl, sex, age_years):
    """
    4-variable MDRD equation, IDMS-standardised creatinine, no race adjustment.
    Reference: Levey et al., Ann Intern Med 2006, doi:10.7326/0003-4819-145-4-200608150-00004
    """
    return case(
        when(sex == "female").then(
            175 * scr_mgdl ** -1.154 * age_years ** -0.203 * 0.742
        ),
        when(sex == "male").then(
            175 * scr_mgdl ** -1.154 * age_years ** -0.203
        ),
        otherwise=None,
    )


def apply_ckd_criteria(event_frame, egfr_series):
    """
    Apply the two-reading CKD Stage 3+ criteria to any measurement event frame.

    event_frame  — EventFrame of measurements (creatinine or recorded eGFR)
    egfr_series  — EventSeries giving the eGFR value (ml/min/1.73m²) for each
                   row; may be computed (creatinine strategies) or the raw
                   numeric_value (recorded-eGFR strategy)

    Returns (has_ckd, index_date, index_numeric_value) where index_numeric_value
    is the raw numeric_value of the index event (creatinine in µmol/L, or
    recorded eGFR in ml/min/1.73m², depending on the frame used).
    """
    low_events = event_frame.where(egfr_series <= CKD_EGFR_THRESHOLD)
    latest = low_events.sort_by(clinical_events.date).last_for_patient()
    has_confirming = low_events.where(
        clinical_events.date.is_before(latest.date - PERSISTENCE_GAP)
    ).exists_for_patient()
    return latest.exists_for_patient() & has_confirming, latest.date, latest.numeric_value


# --- Creatinine-derived eGFR strategies ---
creatinine_events = (
    clinical_events.where(clinical_events.snomedct_code.is_in(creatinine_codelist))
    .where(clinical_events.numeric_value.is_not_null())
    .where(clinical_events.date.is_on_or_before(study_end_date))
)

scr_mgdl = creatinine_events.numeric_value / UMOL_TO_MGDL
age_years = (creatinine_events.date - patients.date_of_birth).years.as_float()

event_egfr_ckd_epi = ckd_epi_2021(scr_mgdl, patients.sex, age_years)
event_egfr_mdrd = mdrd_4var(scr_mgdl, patients.sex, age_years)
event_egfr_era = case(
    when(creatinine_events.date.is_before(FORMULA_CUTOFF_DATE)).then(event_egfr_mdrd),
    otherwise=event_egfr_ckd_epi,
)

has_ckd_epi, _, _ = apply_ckd_criteria(creatinine_events, event_egfr_ckd_epi)
has_mdrd, _, _ = apply_ckd_criteria(creatinine_events, event_egfr_mdrd)
has_era, index_creatinine_date, index_creatinine_umol = apply_ckd_criteria(
    creatinine_events, event_egfr_era
)

# --- Recorded eGFR strategy ---
egfr_events = (
    clinical_events.where(clinical_events.snomedct_code.is_in(egfr_codelist))
    .where(clinical_events.numeric_value.is_not_null())
    .where(clinical_events.date.is_on_or_before(study_end_date))
)

has_recorded_egfr, recorded_egfr_index_date, recorded_egfr_index_value = (
    apply_ckd_criteria(egfr_events, egfr_events.numeric_value)
)

# --- Population ---
is_registered = practice_registrations.for_patient_on(
    study_end_date
).exists_for_patient()

dataset.define_population(
    is_registered & (has_ckd_epi | has_mdrd | has_era | has_recorded_egfr)
)

# --- Output columns ---
dataset.sex = patients.sex
dataset.date_of_birth = patients.date_of_birth

# CKD flags — researcher picks the strategy appropriate for their study
dataset.ckd_by_ckd_epi_2021 = has_ckd_epi
dataset.ckd_by_mdrd = has_mdrd
dataset.ckd_by_era_appropriate = has_era
dataset.ckd_by_recorded_egfr = has_recorded_egfr

# Index measurement from the era-appropriate creatinine strategy
dataset.index_creatinine_date = index_creatinine_date
dataset.index_creatinine_umol = index_creatinine_umol

# Index measurement from the recorded-eGFR strategy (ml/min/1.73m²)
dataset.recorded_egfr_index_date = recorded_egfr_index_date
dataset.recorded_egfr_index_value = recorded_egfr_index_value
