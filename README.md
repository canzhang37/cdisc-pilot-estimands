# CDISC pilot study: ADaM with admiral, estimands and reference-based imputation

A regulatory analysis of a clinical trial runs from raw SDTM domains to analysis datasets, from a statistical analysis plan to an estimand, and from a primary model to sensitivity analyses that ask how much the answer depends on what is assumed about patients who stop treatment. This project does that end to end in R on the public CDISC pilot study, a 24-week trial of two doses of xanomeline against placebo in mild to moderate Alzheimer's disease.

The ADSL, ADADAS (ADAS-Cog) and ADAE datasets are derived from SDTM with the admiral package and written as SAS transport (XPT) files. They are then checked by independent derivation against the pilot's own ADaM datasets, variable by variable, with every mismatch traced to its cause. The analysis follows a statistical analysis plan written before any result was computed: an ICH E9(R1) estimand with a hypothetical strategy for treatment discontinuation, an MMRM with unstructured covariance and Kenward-Roger degrees of freedom, a fixed-sequence test of the two doses, and sensitivity analyses with reference-based imputation (jump to reference, copy reference) from the rbmi package and with the last observation carried forward. Tables are produced with rtables.

The plan is [docs/sap.md](docs/sap.md) and the readout is [docs/readout.md](docs/readout.md). Neither dose separated from placebo: the high-dose difference in ADAS-Cog(11) change at Week 24 was -0.83 points (95% CI -3.21 to 1.55). Jump-to-reference imputation moved it 43% closer to zero and LOCF moved it away from zero, and no analysis excluded zero. All ADAS-Cog analysis records, their changes from baseline and all adverse-event flags matched the pilot ADaM; 18 of 19 ADSL variables agreed for at least 95% of patients.

```
$ Rscript scripts/06_verdicts.R
  id check                                                    value            met
1 P1 lowest ADSL agreement, tolerance rule (>= 95%)           88.6% (BMIBL)    FALSE
2 P1 lowest ADSL agreement, R rounding rule (>= 95%)          88.6% (BMIBL)    FALSE
3 P2 recomputed = collected total when 11 items present (all) 797 of 797       TRUE
4 P3 mmrm vs gls: max |d est| (< 1e-4), max |d se| (< 1e-3)   4.3e-06, 1.4e-05 TRUE
5 P4 |MAR - MMRM|, high dose (<= 0.3)                         0.02             TRUE
6 P5 J2R shrinkage toward 0, high dose (>= 20%)               43%              TRUE
7 P6 LOCF minus MMRM, high dose (< 0, more favorable)         -0.18            TRUE
```

![sensitivity analyses](figures/sensitivity_forest.png)

## Run

From the project folder, with R 4.4 and the packages admiral, safetyData, mmrm, emmeans, rbmi, rtables, haven, dplyr, tidyr, ggplot2 and testthat:

```
Rscript scripts/01_adam.R
Rscript scripts/02_qc.R
Rscript scripts/03_efficacy.R
Rscript scripts/04_tlf.R
Rscript scripts/05_figures.R
Rscript scripts/06_verdicts.R
Rscript -e "testthat::test_dir('tests')"
```

The whole pipeline runs in a few minutes on a laptop; the jackknife in `03_efficacy.R` refits the imputation model once per patient for each strategy. Results were produced with admiral 1.4.1, mmrm 0.3.17, rbmi 1.6.0 and rtables 0.6.16.

## Data

SDTM and ADaM data of the CDISC SDTM/ADaM pilot project (study CDISCPILOT01), as distributed in the R package safetyData 1.0.0 (MIT license). `adam/` holds RDS copies of the derived ADaM datasets; `01_adam.R` also writes them there as XPT files. The code is under the MIT license.

Can Zhang, zhan8639@umn.edu
