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

from ehrql import create_dataset
from ehrql.tables.tpp import patients, practice_registrations

from analysis.ckd_phenotype import ckd_from_creatinine, ckd_from_recorded_egfr

dataset = create_dataset()

study_end_date = "2020-03-31"

has_ckd_epi, has_mdrd, has_era, index_creatinine_date, index_creatinine_umol = (
    ckd_from_creatinine(study_end_date)
)
has_recorded_egfr, recorded_egfr_index_date, recorded_egfr_index_value = (
    ckd_from_recorded_egfr(study_end_date)
)

is_registered = practice_registrations.for_patient_on(study_end_date).exists_for_patient()

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
