# Youth labor-market transitions in Northeast Brazil (PNADc, 2012-2024)

This repository estimates quarter-to-quarter labor-market and school-to-work transitions of young people (aged 15 to 29) in three Northeastern Brazilian states: Pernambuco, Ceará and Bahia. It uses the longitudinal structure of Brazil's continuous household survey (PNADc, IBGE), in which the same dwellings are interviewed in five consecutive quarters.

The project started during my research assistantship at the Observatório do Trabalho, a labor-market research unit of the Government of Pernambuco. All data are public, so the results can be reproduced in full.

![Where unemployed youth are one quarter later, Pernambuco](figures/fig1_unemployed_PE.png)

## What the code does

1. **Data pipeline.** Reads 13 years of partitioned PNADc parquet files with `arrow`, pushing the state filter down to the file reader, and assembles about 3.8 million person-quarter records.
2. **Longitudinal weights.** Links individuals across consecutive quarters and builds two-step longitudinal weights following Monteiro (2019):
   1. an attrition adjustment within each PNADc estimation domain (`posest`);
   2. a calibration by domain, sex and age group to the full cross-sectional population.

   The 15-29 age restriction is applied only after the weights are built, so the calibration targets cover the whole population of each domain.
3. **Transition rates and figures.** Computes weighted transition rates (for example, where young people who were unemployed are one quarter later) and produces one figure per state.
4. **Built-in checks.** Reports match and attrition rates for every pair of quarters, verifies that the weighted pairs reproduce the calibration targets in every cell, summarizes the weight distribution, and reports unclassified observations instead of dropping them silently.

## Design choices

1. **Linking only within a calendar year.** The individual identifier (`id_rs`) is built separately for each annual download, so the same identifier in Q4 of one year and Q1 of the next does not refer to the same person. In the data, cross-year matches agree on sex in only about half of the cases, against more than 99.9% within the year. The code therefore uses only Q1-Q2, Q2-Q3 and Q3-Q4 transitions.
2. **Eligibility for re-interview.** The interview number (`V1016`) is not in the files, but the PNADc rotates whole primary sampling units (`UPA`) out of a panel after their fifth visit. A person counts as eligible at t+1 if their sampling unit is surveyed again at t+1. Eligible people who are not found are treated as attrition. A sampling unit with complete non-response at t+1 is treated as rotated out, which is rare.
3. **Unlinked records are kept.** About 12% of records have no longitudinal identifier. They stay in the cross-section used as calibration target and in the eligible population, rather than being discarded.

## Data

PNADc quarterly microdata (IBGE), downloaded one calendar year at a time with [`datazoom.social`](https://github.com/datazoompuc/datazoom.social):

```r
datazoom.social::load_pnadc(years = 2012, quarters = 1:4, panel = "advanced", raw_data = TRUE)
```

The code expects the files in the following layout, where `V1014` is the PNADc panel:

```
input/<year>/pnadc_panels/V1014=<panel>/part-0.parquet
```

The microdata are not included in this repository because of their size.


| File | Content |
|---|---|
| `longitudinal_pairs_PE_CE_BA.parquet` | person-level t to t+1 pairs with longitudinal weights |
| `diagnostics_pairs.csv` | eligible persons, matches and match rate for each pair of quarters |
| `diagnostics_weights.csv` | weight distribution |
| `transitions_unemployed.csv`, `transitions_students.csv` | weighted transition rates |
| `fig1_unemployed_<UF>.png`, `fig2_students_<UF>.png` | figures by state |

## Repository structure

```
.
├── youth_transitions_pnadc.R   # full pipeline: read, link, weight, check, plot
├── figures/                    # example output for Pernambuco
└── README.md
```

## Author

Julia Moura de Lima, Department of Economics, Universidade Federal de Pernambuco (UFPE).
