# Statistical analysis plan: ADAS-Cog(11) at Week 24 in the CDISC pilot study

This plan was written after reading the structure of the data (variable lists, test codes, visit schedule, analysis windows, population flags and the list of discontinuation reasons) and before any efficacy or safety result was computed. Changes made after results are seen go in the Amendments section at the end, with a reason, and are labeled in the readout.

## 1. Study and data

CDISCPILOT01 is the public CDISC SDTM/ADaM pilot study: a randomized, double-blind, placebo-controlled 24-week trial of xanomeline transdermal therapeutic system at a low dose (54 mg/day) and a high dose (81 mg/day) in mild to moderate Alzheimer's disease, with 254 randomized patients (placebo 86, low dose 84, high dose 84). The SDTM domains and the pilot's own ADaM datasets are taken from the R package safetyData 1.0.0 (MIT license). The analysis datasets in this project are derived again from SDTM; the pilot's ADaM datasets are used only as an independent reference for quality control.

## 2. Analysis datasets

Derived from SDTM with the admiral package and exported as SAS transport (XPT version 5) files.

ADSL, one record per patient:
- Treatment: TRT01P from ARM, TRT01A from ACTARM; TRTSDT and TRTEDT from the first and last EX records.
- Populations from SUPPDM: ITTFL (ITT), SAFFL (SAFETY), EFFFL (EFFICACY), COMP24FL (COMPLT24).
- Demographics from DM: AGE, AGEGR1 (<65, 65-80, >80), SEX, RACE, ETHNIC.
- Baseline characteristics: HEIGHTBL and WEIGHTBL from VS (last value on or before the first dose), BMIBL = weight / height^2, EDUCLVL from SC, MMSETOT (Mini-Mental State total at screening) from QS.
- Disposition from the DS disposition events: EOSSTT (completed or discontinued) and DCSREAS (reason).
- SITEGR1: sites pooled into one group when they enrolled fewer than 3 patients in any arm.

ADADAS, one record per patient, parameter and analysis visit:
- Parameters ACITM01 to ACITM14 and ACTOT (ADAS-Cog(11) total), AVAL from QSSTRESN.
- ADY from QSDTC and TRTSDT. Analysis visits by study-day windows: Baseline day <= 1, Week 8 days 2-84 (target 56), Week 16 days 85-140 (target 112), Week 24 day >= 141 (target 168).
- ANL01FL marks one record per patient, parameter and window: the one closest to the target day, the later one on a tie. Baseline is the last record on or before day 1.
- BASE and CHG = AVAL - BASE.
- A recomputed total, the sum of the 11 ADAS-Cog(11) items (ACITM01, 02, 04, 05, 06, 07, 08, 11, 12, 13, 14), is kept for checking against the collected total.

ADAE, one record per adverse event:
- ASTDT from AESTDTC with a missing day imputed as the first of the month and a missing month as January, but never before TRTSDT within the same month or year.
- TRTEMFL = Y when ASTDT is on or after TRTSDT.

## 3. Quality control by independent derivation

The derived datasets are compared with the pilot's own ADaM datasets:
- ADSL: TRT01P, TRT01A, TRTSDT, TRTEDT, AGE, AGEGR1, SEX, RACE, ITTFL, SAFFL, EFFFL, COMP24FL, reason for discontinuation, HEIGHTBL, WEIGHTBL, BMIBL, EDUCLVL, MMSETOT and SITEGR1; agreement per variable, with every mismatch explained.
- ADADAS: ACTOT values, BASE, CHG, analysis visit and ANL01FL against the observed (non-LOCF) records of the pilot ADQSADAS.
- ADAE: TRTEMFL.
- The recomputed ADAS-Cog(11) total against the collected total.
- MMRM software: the primary model fitted with `mmrm` (Satterthwaite degrees of freedom) and with `nlme::gls` (unstructured correlation, visit-specific variances, REML) must give the same treatment differences to 1e-4 and standard errors to 1e-3.

## 4. Estimands

Primary estimand (ICH E9(R1)):
- Population: patients with mild to moderate Alzheimer's disease, represented by the efficacy population (EFFFL = Y).
- Treatment condition: xanomeline high dose versus placebo for 24 weeks; low dose versus placebo as the second comparison.
- Variable: change from baseline in ADAS-Cog(11) total at Week 24 (higher scores are worse).
- Intercurrent event: discontinuation of study treatment for any reason. In this study, patients who stop treatment also leave the study, so no later measurements exist.
- Strategy: hypothetical, the effect if patients had stayed on their assigned treatment.
- Population-level summary: difference in mean change from baseline at Week 24, active minus placebo.

Supplementary estimand: the same, except that patients who discontinue an active arm are assumed to follow the placebo course after discontinuation (reference-based imputation), which is closer to a treatment-policy view of a treatment that stops.

## 5. Primary analysis

Mixed model for repeated measures (MMRM) on CHG at Weeks 8, 16 and 24 (ANL01FL = Y records), with fixed effects for treatment, visit, treatment by visit, BASE, BASE by visit and SITEGR1. Unstructured covariance, REML, Kenward-Roger degrees of freedom. If the unstructured model does not converge, the first of heterogeneous Toeplitz, heterogeneous AR(1) and compound symmetry that converges is used. Least-squares means at Week 24 by arm (proportional weights over site groups), differences from placebo with 95% confidence intervals and two-sided p values.

Multiplicity: fixed sequence. High dose versus placebo is tested at a two-sided 0.05 level; low dose versus placebo is declared significant only if the first test is.

## 6. Sensitivity analyses

With the `rbmi` package, conditional mean imputation with jackknife standard errors, the imputation model having the same fixed effects as the primary MMRM, and an ANCOVA at Week 24 on BASE and SITEGR1 as the analysis model:
- S1. Missing at random (MAR).
- S2. Jump to reference (J2R, reference = placebo) for visits after discontinuation in the active arms; intermittent missing visits stay MAR.
- S3. Copy reference (CR), otherwise as S2.
- S4. Tipping point: under MAR, a shift delta is added to every imputed value after discontinuation in the active arm, from 0 to 10 points in steps of 0.5; the result is the smallest delta at which the 95% confidence interval includes 0. Done only for comparisons significant in the primary analysis.
- S5. The pilot's original approach: ANCOVA at Week 24 on BASE and SITEGR1 with the last observation carried forward (LOCF).

## 7. Tables and figures

- T1 demographics and baseline characteristics (ITT).
- T2 disposition (ITT).
- T3 primary MMRM results (efficacy population).
- T4 sensitivity analyses.
- T5 treatment-emergent adverse events: an overview and events by system organ class and preferred term for preferred terms in at least 5% of any arm (safety population).
- F1 least-squares mean change from baseline over time.
- F2 tipping point.

## 8. Predictions

- P1. Each compared ADSL variable agrees for at least 95% of patients, and every mismatch has an explanation.
- P2. The recomputed ADAS-Cog(11) total equals the collected total whenever all 11 items are present.
- P3. `mmrm` and `nlme::gls` agree as stated in section 3.
- P4. MI under MAR (S1) is within 0.3 points of the MMRM estimate for high dose versus placebo.
- P5. J2R (S2) moves the high-dose effect at least 20% closer to zero than the MMRM estimate. This prediction is informed by the large number of discontinuations for adverse events seen in the ADSL structure.
- P6. LOCF (S5) gives a high-dose effect more favorable to the high dose than the MMRM estimate.

## Amendments

Amendment 1 (after the first QC run). Comparison rules in the QC: numeric ADSL variables agree when they differ by at most 0.05, because the pilot stores values rounded half up (as SAS does) while R rounds half to even; agreement under plain rounding is reported as well. CHG is compared after baseline only, because the pilot leaves CHG missing at baseline. TRTEMFL is compared as Y versus not Y, because the pilot uses Y/N and this project Y/blank.

Amendment 2 (after the first sensitivity run). The first run of the imputation model left out the treatment-by-visit interaction, which `rbmi` does not add by itself; the corrected model has the fixed effects this plan specifies. A check was added: without the site term, conditional mean imputation under MAR followed by a Week 24 ANCOVA must reproduce the MMRM estimate exactly.

Amendment 3 (after the primary analysis). No comparison was significant, so the tipping-point analysis (S4) was not run, as specified. Figure F2 shows a forest plot of all analyses instead.
