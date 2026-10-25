# double programming vs the pilot ADaM
suppressPackageStartupMessages(library(dplyr))
source("R/adam.R")
adsl <- readRDS("adam/adsl.rds")
adadas <- readRDS("adam/adadas.rds")
adae <- readRDS("adam/adae.rds")
ref_sl <- safetyData::adam_adsl
ref_as <- safetyData::adam_adqsadas
ref_ae <- safetyData::adam_adae

## ADSL; numbers within 0.05 (SAS rounding)
reason <- c("Adverse Event" = "ADVERSE EVENT", "Death" = "DEATH", "Lack of Efficacy" = "LACK OF EFFICACY",
            "Lost to Follow-up" = "LOST TO FOLLOW-UP", "Physician Decision" = "PHYSICIAN DECISION",
            "Protocol Violation" = "PROTOCOL VIOLATION", "Sponsor Decision" = "STUDY TERMINATED BY SPONSOR",
            "Withdrew Consent" = "WITHDRAWAL BY SUBJECT", "I/E Not Met" = "PROTOCOL VIOLATION")
ref <- ref_sl %>% mutate(DCSREAS = unname(reason[DCREASCD]), SITEGR1 = as.character(SITEGR1))
vars <- c("TRT01P", "TRT01A", "TRTSDT", "TRTEDT", "AGE", "AGEGR1", "SEX", "RACE", "ITTFL", "SAFFL", "EFFFL",
          "COMP24FL", "DCSREAS", "HEIGHTBL", "WEIGHTBL", "BMIBL", "EDUCLVL", "MMSETOT", "SITEGR1")
j <- inner_join(adsl, ref, by = "USUBJID", suffix = c("", ".ref"))
same <- function(a, b, rule = "tol") {
  if (is.numeric(a) && is.numeric(b)) {
    eq <- if (rule == "tol") abs(a - b) <= 0.05 + 1e-9 else round(a, 1) == round(b, 1)
    return(ifelse(is.na(a) & is.na(b), TRUE, eq))
  }
  a <- as.character(a); b <- as.character(b)
  ifelse(is.na(a) & is.na(b), TRUE, a == b)
}
qc_sl <- bind_rows(lapply(vars, function(v) {
  s <- same(j[[v]], j[[paste0(v, ".ref")]])
  r <- same(j[[v]], j[[paste0(v, ".ref")]], rule = "round")
  s[is.na(s)] <- FALSE; r[is.na(r)] <- FALSE
  data.frame(variable = v, n = nrow(j), agree = sum(s), pct = 100 * mean(s), pct_r_round = 100 * mean(r))
}))
mism <- bind_rows(lapply(vars, function(v) {
  s <- same(j[[v]], j[[paste0(v, ".ref")]]); s[is.na(s)] <- FALSE
  if (all(s)) return(NULL)
  data.frame(variable = v, USUBJID = j$USUBJID[!s], derived = as.character(j[[v]][!s]),
             reference = as.character(j[[paste0(v, ".ref")]][!s]))
}))
dir.create("results", showWarnings = FALSE)
write.csv(qc_sl, "results/qc_adsl.csv", row.names = FALSE)
write.csv(mism, "results/qc_adsl_mismatches.csv", row.names = FALSE)
cat("ADSL: patients matched", nrow(j), "of", nrow(adsl), "\n"); print(qc_sl)

## ADAS-Cog(11) total: items vs collected
tot <- adadas %>%
  filter(PARAMCD %in% c(ITEMS11, "ACTOT")) %>%
  group_by(USUBJID, ADT) %>%
  summarise(n_items = sum(PARAMCD %in% ITEMS11 & !is.na(AVAL)), recomputed = sum(AVAL[PARAMCD %in% ITEMS11]),
            collected = AVAL[PARAMCD == "ACTOT"][1], .groups = "drop")
full <- filter(tot, n_items == 11, !is.na(collected))
part <- filter(tot, n_items < 11, !is.na(collected))
qc_tot <- data.frame(check = c("all 11 items present", "fewer than 11 items"),
                     records = c(nrow(full), nrow(part)),
                     equal = c(sum(full$recomputed == full$collected), sum(part$recomputed == part$collected, na.rm = TRUE)))
write.csv(qc_tot, "results/qc_adas_total.csv", row.names = FALSE)
cat("\nADAS-Cog(11) total\n"); print(qc_tot)

## ADADAS ACTOT vs observed pilot records
mine <- adadas %>% filter(PARAMCD == "ACTOT", ANL01FL == "Y") %>% select(USUBJID, AVISIT, AVAL, BASE, CHG, ADY)
theirs <- ref_as %>% filter(PARAMCD == "ACTOT", DTYPE == "" | is.na(DTYPE), ANL01FL == "Y") %>%
  select(USUBJID, AVISIT, AVAL, BASE, CHG, ADY)
m <- full_join(mine, theirs, by = c("USUBJID", "AVISIT"), suffix = c("", ".ref"))
qc_as <- data.frame(check = c("analysis records, derived", "analysis records, pilot", "in both",
                              "AVAL equal", "BASE equal", "CHG equal", "same study day"),
                    n = c(nrow(mine), nrow(theirs), sum(!is.na(m$AVAL) & !is.na(m$AVAL.ref)),
                          sum(abs(m$AVAL - m$AVAL.ref) < 1e-6, na.rm = TRUE), sum(abs(m$BASE - m$BASE.ref) < 1e-6, na.rm = TRUE),
                          sum(abs(m$CHG - m$CHG.ref) < 1e-6, na.rm = TRUE), sum(m$ADY == m$ADY.ref, na.rm = TRUE)))
qc_as$n[qc_as$check == "CHG equal"] <- paste0(qc_as$n[qc_as$check == "CHG equal"], " (pilot CHG missing at baseline: ",
                                              sum(m$AVISIT == "Baseline"), " records)")
write.csv(qc_as, "results/qc_adadas.csv", row.names = FALSE)
write.csv(filter(m, is.na(AVAL) | is.na(AVAL.ref) | abs(AVAL - AVAL.ref) >= 1e-6 | ADY != ADY.ref), "results/qc_adadas_mismatches.csv", row.names = FALSE)
cat("\nADADAS ACTOT\n"); print(qc_as)

## ADAE TRTEMFL (pilot Y/N, here Y/blank)
ka <- adae %>% select(USUBJID, AESEQ, TRTEMFL)
kr <- ref_ae %>% select(USUBJID, AESEQ, TRTEMFL.ref = TRTEMFL) %>% mutate(TRTEMFL.ref = na_if(TRTEMFL.ref, ""))
ae <- inner_join(ka, kr, by = c("USUBJID", "AESEQ"))
qc_ae <- data.frame(records_derived = nrow(ka), records_pilot = nrow(kr), matched = nrow(ae),
                    TRTEMFL_equal = sum((ae$TRTEMFL %in% "Y") == (ae$TRTEMFL.ref %in% "Y")))
write.csv(qc_ae, "results/qc_adae.csv", row.names = FALSE)
cat("\nADAE\n"); print(qc_ae)
