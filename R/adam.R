# ADSL, ADADAS, ADAE from SDTM with admiral
suppressPackageStartupMessages({
  library(admiral)
  library(dplyr)
})

ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
ITEMS11 <- sprintf("ACITM%02d", c(1, 2, 4, 5, 6, 7, 8, 11, 12, 13, 14))
WIN <- data.frame(AVISIT = c("Baseline", "Week 8", "Week 16", "Week 24"),
                  AVISITN = c(0, 8, 16, 24), AWTARGET = c(1, 56, 112, 168))

sdtm <- function(x) convert_blanks_to_na(getExportedValue("safetyData", paste0("sdtm_", x)))

## ADSL
derive_adsl <- function() {
  supp <- sdtm("suppdm") %>%
    select(USUBJID, QNAM, QVAL) %>%
    tidyr::pivot_wider(names_from = QNAM, values_from = QVAL)
  ex <- sdtm("ex") %>%
    derive_vars_dt(new_vars_prefix = "EXST", dtc = EXSTDTC) %>%
    derive_vars_dt(new_vars_prefix = "EXEN", dtc = EXENDTC)
  ds <- sdtm("ds") %>% filter(DSCAT == "DISPOSITION EVENT", DSDECOD != "SCREEN FAILURE")
  flag <- function(x) if_else(x %in% "Y", "Y", "N")

  adsl <- sdtm("dm") %>%
    filter(ARM %in% ARMS) %>%
    left_join(supp, by = "USUBJID") %>%
    mutate(TRT01P = ARM, TRT01A = ACTARM,
           TRT01PN = match(TRT01P, ARMS) - 1, TRT01AN = match(TRT01A, ARMS) - 1,
           ITTFL = flag(ITT), SAFFL = flag(SAFETY), EFFFL = flag(EFFICACY), COMP24FL = flag(COMPLT24),
           AGEGR1 = case_when(AGE < 65 ~ "<65", AGE <= 80 ~ "65-80", AGE > 80 ~ ">80")) %>%
    derive_vars_merged(dataset_add = ex, by_vars = exprs(STUDYID, USUBJID), new_vars = exprs(TRTSDT = EXSTDT),
                       order = exprs(EXSTDT, EXSEQ), mode = "first", filter_add = !is.na(EXSTDT)) %>%
    derive_vars_merged(dataset_add = ex, by_vars = exprs(STUDYID, USUBJID), new_vars = exprs(TRTEDT = EXENDT),
                       order = exprs(EXENDT, EXSEQ), mode = "last", filter_add = !is.na(EXENDT)) %>%
    derive_vars_merged(dataset_add = ds, by_vars = exprs(STUDYID, USUBJID), new_vars = exprs(DSDECOD),
                       order = exprs(DSSEQ), mode = "last") %>%
    mutate(EOSSTT = if_else(DSDECOD == "COMPLETED", "COMPLETED", "DISCONTINUED"),
           DCSREAS = if_else(DSDECOD == "COMPLETED", NA_character_, DSDECOD)) %>%
    select(-DSDECOD, -any_of(c("ITT", "SAFETY", "EFFICACY", "COMPLT8", "COMPLT16", "COMPLT24")))

  # baseline height/weight: last before 1st dose
  vs <- sdtm("vs") %>%
    filter(VSTESTCD %in% c("HEIGHT", "WEIGHT"), !is.na(VSSTRESN)) %>%
    derive_vars_dt(new_vars_prefix = "VS", dtc = VSDTC) %>%
    inner_join(select(adsl, USUBJID, TRTSDT), by = "USUBJID") %>%
    filter(!is.na(TRTSDT), VSDT <= TRTSDT)
  for (tc in c("HEIGHT", "WEIGHT")) {
    adsl <- derive_vars_merged(adsl, dataset_add = filter(vs, VSTESTCD == tc), by_vars = exprs(STUDYID, USUBJID),
                               new_vars = setNames(exprs(VSSTRESN), paste0(tc, "BL")), order = exprs(VSDT, VSSEQ), mode = "last")
  }
  mmse <- sdtm("qs") %>%
    filter(QSCAT == "MINI-MENTAL STATE") %>%
    group_by(STUDYID, USUBJID) %>%
    summarise(MMSETOT = sum(QSSTRESN), .groups = "drop")
  edu <- sdtm("sc") %>% filter(SCTESTCD == "EDLEVEL") %>% select(STUDYID, USUBJID, EDUCLVL = SCSTRESN)

  adsl <- adsl %>%
    mutate(BMIBL = WEIGHTBL / (HEIGHTBL / 100)^2) %>%
    derive_vars_merged(dataset_add = edu, by_vars = exprs(STUDYID, USUBJID)) %>%
    derive_vars_merged(dataset_add = mmse, by_vars = exprs(STUDYID, USUBJID))

  # pool sites with < 3 patients in any arm
  small <- adsl %>%
    count(SITEID, TRT01P) %>%
    tidyr::complete(SITEID, TRT01P, fill = list(n = 0)) %>%
    group_by(SITEID) %>%
    summarise(pool = min(n) < 3)
  adsl %>%
    left_join(small, by = "SITEID") %>%
    mutate(SITEGR1 = if_else(pool, "900", as.character(SITEID))) %>%
    select(-pool)
}

## ADADAS
derive_adadas <- function(adsl) {
  sdtm("qs") %>%
    filter(QSCAT == "ALZHEIMER'S DISEASE ASSESSMENT SCALE") %>%
    derive_vars_merged(dataset_add = adsl, by_vars = exprs(STUDYID, USUBJID),
                       new_vars = exprs(TRTSDT, TRTP = TRT01P, TRTPN = TRT01PN, SITEGR1, ITTFL, EFFFL)) %>%
    filter(!is.na(TRTSDT)) %>%
    mutate(PARAMCD = QSTESTCD, PARAM = QSTEST, AVAL = QSSTRESN) %>%
    derive_vars_dt(new_vars_prefix = "A", dtc = QSDTC) %>%
    derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT)) %>%
    mutate(AVISIT = case_when(ADY <= 1 ~ "Baseline", ADY <= 84 ~ "Week 8", ADY <= 140 ~ "Week 16", TRUE ~ "Week 24")) %>%
    left_join(WIN, by = "AVISIT") %>%
    mutate(AWTDIFF = abs(ADY - AWTARGET)) %>%
    restrict_derivation(derivation = derive_var_extreme_flag,
                        args = params(by_vars = exprs(STUDYID, USUBJID, PARAMCD), order = exprs(ADT, QSSEQ),
                                      new_var = ABLFL, mode = "last"),
                        filter = AVISIT == "Baseline" & !is.na(AVAL)) %>%
    restrict_derivation(derivation = derive_var_extreme_flag,
                        args = params(by_vars = exprs(STUDYID, USUBJID, PARAMCD, AVISIT), order = exprs(AWTDIFF, desc(ADY), QSSEQ),
                                      new_var = ANL01FL, mode = "first"),
                        filter = AVISIT != "Baseline" & !is.na(AVAL)) %>%
    mutate(ANL01FL = if_else(ABLFL %in% "Y", "Y", ANL01FL)) %>%
    derive_var_base(by_vars = exprs(STUDYID, USUBJID, PARAMCD), source_var = AVAL, new_var = BASE) %>%
    derive_var_chg() %>%
    select(STUDYID, USUBJID, SITEGR1, TRTP, TRTPN, TRTSDT, ITTFL, EFFFL, PARAMCD, PARAM, AVAL, BASE, CHG,
           ADT, ADY, AVISIT, AVISITN, AWTARGET, AWTDIFF, ABLFL, ANL01FL, VISIT, VISITNUM, QSSEQ)
}

## ADAE
derive_adae <- function(adsl) {
  sdtm("ae") %>%
    derive_vars_merged(dataset_add = adsl, by_vars = exprs(STUDYID, USUBJID),
                       new_vars = exprs(TRTSDT, TRTEDT, TRTA = TRT01A, TRTAN = TRT01AN, SAFFL)) %>%
    derive_vars_dt(new_vars_prefix = "AST", dtc = AESTDTC, highest_imputation = "M",
                   date_imputation = "first", min_dates = exprs(TRTSDT)) %>%
    mutate(TRTEMFL = if_else(!is.na(ASTDT) & !is.na(TRTSDT) & ASTDT >= TRTSDT, "Y", NA_character_)) %>%
    select(STUDYID, USUBJID, TRTA, TRTAN, SAFFL, TRTSDT, TRTEDT, AESEQ, AETERM, AEDECOD, AEBODSYS, AESOC,
           AESEV, AESER, AEREL, AESTDTC, ASTDT, ASTDTF, TRTEMFL)
}
