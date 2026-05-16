"""
Shared CKD Stage 3+ phenotype components — reusable across dataset definitions.

Provides the CALIBER/PH42 codelists, constants, eGFR formula functions, and
the two-reading CKD criteria helper. Import from this module rather than
duplicating these definitions across studies.
"""

from ehrql import case, codelist_from_csv, days, when
from ehrql.tables.tpp import clinical_events, patients

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


def ckd_from_creatinine(study_end_date):
    """
    Run all three creatinine-derived CKD strategies up to study_end_date.

    Returns (has_ckd_epi, has_mdrd, has_era, index_creatinine_date, index_creatinine_umol)
    where index values are from the era-appropriate strategy.
    """
    creatinine_events = (
        clinical_events.where(clinical_events.snomedct_code.is_in(creatinine_codelist))
        .where(clinical_events.numeric_value.is_not_null())
        .where(clinical_events.date.is_on_or_before(study_end_date))
    )
    scr_mgdl = creatinine_events.numeric_value / UMOL_TO_MGDL
    age_years = (creatinine_events.date - patients.date_of_birth).years.as_float()

    egfr_ckd_epi = ckd_epi_2021(scr_mgdl, patients.sex, age_years)
    egfr_mdrd = mdrd_4var(scr_mgdl, patients.sex, age_years)
    egfr_era = case(
        when(creatinine_events.date.is_before(FORMULA_CUTOFF_DATE)).then(egfr_mdrd),
        otherwise=egfr_ckd_epi,
    )

    has_ckd_epi, _, _ = apply_ckd_criteria(creatinine_events, egfr_ckd_epi)
    has_mdrd, _, _ = apply_ckd_criteria(creatinine_events, egfr_mdrd)
    has_era, index_creatinine_date, index_creatinine_umol = apply_ckd_criteria(
        creatinine_events, egfr_era
    )

    return has_ckd_epi, has_mdrd, has_era, index_creatinine_date, index_creatinine_umol


def ckd_from_recorded_egfr(study_end_date):
    """
    Run the recorded-eGFR CKD strategy up to study_end_date.

    Returns (has_recorded_egfr, index_date, index_value) where index_value is in ml/min/1.73m².
    """
    egfr_events = (
        clinical_events.where(clinical_events.snomedct_code.is_in(egfr_codelist))
        .where(clinical_events.numeric_value.is_not_null())
        .where(clinical_events.date.is_on_or_before(study_end_date))
    )
    return apply_ckd_criteria(egfr_events, egfr_events.numeric_value)
