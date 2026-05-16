"""CKD Stage 3+ (PH42) study. Algorithm details: analysis/ckd_phenotype.py"""

from ehrql import create_dataset
from ehrql.tables.tpp import patients, practice_registrations

from analysis.ckd_phenotype import ckd_stage3_plus

dataset = create_dataset()

study_end_date = "2020-03-31"

ckd = ckd_stage3_plus(study_end_date)
is_registered = practice_registrations.for_patient_on(study_end_date).exists_for_patient()

dataset.define_population(is_registered & ckd.any_strategy)

dataset.sex           = patients.sex
dataset.date_of_birth = patients.date_of_birth

dataset.ckd_by_ckd_epi_2021   = ckd.ckd_by_ckd_epi_2021
dataset.ckd_by_mdrd            = ckd.ckd_by_mdrd
dataset.ckd_by_era_appropriate = ckd.ckd_by_era_appropriate
dataset.ckd_by_recorded_egfr   = ckd.ckd_by_recorded_egfr

dataset.index_creatinine_date  = ckd.index_creatinine_date
dataset.index_creatinine_umol  = ckd.index_creatinine_umol

dataset.recorded_egfr_index_date  = ckd.recorded_egfr_index_date
dataset.recorded_egfr_index_value = ckd.recorded_egfr_index_value
