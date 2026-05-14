# CKD Stage 3+ Phenotype (PH42)

Implementation of the [Kuan et al. CALIBER CKD phenotype (PH42)](https://phenotypes.healthdatagateway.org/phenotypes/PH42/version/2933/detail/) in ehrQL for OpenSAFELY TPP backends

**🚨🚨🚨 !!!IMPORTANT: THIS IS AN EXAMPLE ONLY AND HAS NOT BEEN CLINICALLY VALIDATED!!!! DO NOT USE IN REAL RESEARCH!!!! 🚨🚨🚨 **

## Brief

Identifies patients with **Chronic Kidney Disease Stage 3 or above**. A patient is included if, under at least one of the four identification strategies below:

1. Their most recent qualifying eGFR ≤ 60 ml/min/1.73m² before the study end date (the "index measurement").
2. They also have at least one qualifying eGFR ≤ 60 recorded **strictly more than 90 days** before the index measurement — confirming persistent impairment rather than acute kidney injury.

The patient must also be registered at a practice on the study end date.

## Four identification strategies

| Output column | Source | Notes |
|---|---|---|
| `ckd_by_ckd_epi_2021` | Creatinine → 2021 CKD-EPI | Race-free, current recommended standard |
| `ckd_by_mdrd` | Creatinine → 4-var MDRD | Matches what clinicians saw pre-2014 |
| `ckd_by_era_appropriate` | Creatinine → MDRD before 2014, CKD-EPI from 2014 | Closest to original PH42 intent |
| `ckd_by_recorded_egfr` | eGFR recorded directly in EHR | Covers patients with no creatinine code; formula-agnostic |

All four flags are populated for every patient in the output. Researchers should choose the strategy appropriate for their study and filter accordingly.

### How eGFR is obtained under each strategy

**Creatinine-derived strategies** (`ckd_by_ckd_epi_2021`, `ckd_by_mdrd`, `ckd_by_era_appropriate`)

Raw serum creatinine values (in µmol/L, converted to mg/dL for the formulas) are taken from the UKRR creatinine codelist and passed through a formula to compute eGFR. This approach captures patients whose creatinine was measured but whose GP system did not separately record an eGFR code.

- **MDRD (4-variable, IDMS-standardised)** — the formula used by UK labs from ~2006 until NICE CG182 (2014) recommended switching. It tends to underestimate eGFR, particularly above 60 ml/min, so it classifies more patients as CKD Stage 3+ than CKD-EPI does. Use this strategy if your study concerns the pre-2014 era and you want to match the classification patients and clinicians were actually working with at the time.
- **2021 CKD-EPI (race-free)** — the current NICE-recommended formula, more accurate across the full GFR range and without a race coefficient. Use this for contemporary studies or where a single consistent formula across the full study period is preferable.
- **Era-appropriate** — applies MDRD to measurements before 2014-01-01 and CKD-EPI 2021 from that date onward, matching the formula switch recommended by NICE CG182. This is the closest approximation to the original PH42 algorithm and to what was actually reported in the EHR at the time of each measurement. Use this as the default when replicating the original phenotype.

Both formulas are sex-specific (male/female only) and use age at the date of each measurement, not current age. Patients with sex recorded as anything other than male or female return `NULL` for all creatinine-derived strategies and will only appear in the output if they qualify via `ckd_by_recorded_egfr`.

**Recorded eGFR strategy** (`ckd_by_recorded_egfr`)

eGFR values recorded directly in the EHR (using the UKRR eGFR codelist, which includes only codes confirmed to carry numeric values in TPP) are used without any formula. This reflects what the clinician actually saw: the GP system calculated and stored an eGFR at the time of the consultation, often alongside or instead of a separate creatinine code.

This strategy:
- Covers patients whose creatinine was not separately coded but whose eGFR was recorded (common in practices where the system reports eGFR as the primary result)
- Works for patients of any recorded sex, since no formula is applied
- Reflects the formula in use by the clinical system at the time, which may vary by practice and software version

The main limitation is that it cannot be applied consistently across time periods: the recorded eGFR value depends on whichever formula the clinical system used at that point, which researchers cannot audit or standardise retrospectively.

### Which strategy to choose

| Study scenario | Recommended strategy |
|---|---|
| Replicating the original PH42 phenotype | `ckd_by_era_appropriate` |
| Contemporary study, single consistent formula | `ckd_by_ckd_epi_2021` |
| Focus on pre-2014 period, matching clinical decisions of the time | `ckd_by_mdrd` |
| Maximising case ascertainment, including patients with no creatinine code | Union of all four flags |
| Sensitivity analysis | Compare all four; flag discordant cases |

## Generate a dataset

```
PYTHONHASHSEED=0 .venv/bin/ehrql generate-dataset analysis/dataset_definition.py \
  --dummy-tables dummy-tables/ \
  --output dataset.csv
```

## Run assurance tests

```
PYTHONHASHSEED=0 .venv/bin/ehrql assure analysis/test_dataset_definition.py
```

## Dummy data

The dummy tables represent a realistic UK GP-registered population of 5,000 patients
with approximately 6% CKD Stage 3+ prevalence (~300 patients), consistent with QOF
register rates. CKD patients are split across four profiles to exercise all identification
strategies:

| Profile | Count | Expected flags |
|---|---|---|
| `creatinine_and_egfr` | ~150 | All four True |
| `creatinine_only` | ~90 | `ckd_by_recorded_egfr` False |
| `recorded_egfr_only` | ~36 | Only `ckd_by_recorded_egfr` True |
| `borderline_scr` | ~36 | Some creatinine flags may differ (threshold boundary) |

Creatinine values are back-computed from a target eGFR using the inverse 2021 CKD-EPI
formula with ±8% biological noise, so values are physiologically plausible for each
patient's age and sex. Events are generated across the full history back to 2005.

### Regenerate dummy tables

```
uv run python scripts/generate_dummy_tables.py
```

## Codelists

| File | Source | Notes |
|---|---|---|
| `codelists/ukrr-creatinine-tests-level.csv` | [UKRR, 2022](https://www.opencodelists.org/codelist/ukrr/creatinine-tests-level/3ee9db71/) | 5 SNOMED codes confirmed to carry numeric values in TPP |
| `codelists/ukrr-egfr-tests-level.csv` | [UKRR, 2022](https://www.opencodelists.org/codelist/ukrr/egfr-tests-level/5293ceba/) | 6 SNOMED codes confirmed to carry numeric values in TPP |

**All codelists require clinical review before use in production.**
