# rtables layouts for T1-T5
suppressPackageStartupMessages({
  library(rtables)
  library(dplyr)
})

num_rows <- function(x, ...) {
  x <- x[!is.na(x)]
  in_rows("n" = rcell(length(x), format = "xx"),
          "Mean (SD)" = rcell(c(mean(x), sd(x)), format = "xx.x (xx.x)"),
          "Median" = rcell(median(x), format = "xx.x"),
          "Min, Max" = rcell(range(x), format = "xx.x, xx.x"))
}

cat_rows <- function(x, .N_col, ...) {
  tab <- table(x, useNA = "no")
  in_rows(.list = lapply(as.list(tab), function(n) rcell(c(n, n / .N_col), format = "xx (xx.x%)")), .names = names(tab))
}

## T1 demographics, ITT
t_dm <- function(adsl) {
  d <- adsl %>%
    filter(ITTFL == "Y") %>%
    mutate(TRT01P = factor(TRT01P, ARMS), AGEGR1 = factor(AGEGR1, c("<65", "65-80", ">80")),
           SEX = factor(SEX, c("F", "M")), RACE = factor(RACE))
  basic_table(title = "Table 1. Demographics and baseline characteristics (ITT population)") %>%
    split_cols_by("TRT01P") %>%
    add_overall_col("Total") %>%
    add_colcounts() %>%
    analyze("AGE", num_rows, var_labels = "Age (years)", show_labels = "visible") %>%
    analyze("AGEGR1", cat_rows, var_labels = "Age group", show_labels = "visible") %>%
    analyze("SEX", cat_rows, var_labels = "Sex", show_labels = "visible") %>%
    analyze("RACE", cat_rows, var_labels = "Race", show_labels = "visible") %>%
    analyze("BMIBL", num_rows, var_labels = "Baseline BMI (kg/m2)", show_labels = "visible") %>%
    analyze("EDUCLVL", num_rows, var_labels = "Years of education", show_labels = "visible") %>%
    analyze("MMSETOT", num_rows, var_labels = "MMSE total at screening", show_labels = "visible") %>%
    build_table(d)
}

## T2 disposition, ITT
t_ds <- function(adsl) {
  d <- adsl %>%
    filter(ITTFL == "Y") %>%
    mutate(TRT01P = factor(TRT01P, ARMS), EOSSTT = factor(EOSSTT), DCSREAS = factor(tolower(DCSREAS)))
  basic_table(title = "Table 2. Disposition (ITT population)") %>%
    split_cols_by("TRT01P") %>%
    add_overall_col("Total") %>%
    add_colcounts() %>%
    analyze("EOSSTT", cat_rows, var_labels = "End of study status", show_labels = "visible") %>%
    analyze("DCSREAS", cat_rows, var_labels = "Reason for discontinuation", show_labels = "visible") %>%
    build_table(d)
}

## model results as a manual rtable
fmt_ci <- function(e, l, u) sprintf("%.2f (%.2f, %.2f)", e, l, u)

t_primary <- function(lsm, prim) {
  w <- lsm[lsm$AVISIT == "Week 24", ]
  p <- prim[prim$AVISIT == "Week 24", ]
  rows <- lapply(seq_len(nrow(w)), function(i) {
    arm <- as.character(w$TRTP[i])
    k <- which(sub(" - Placebo", "", p$contrast) == arm)
    rrow(arm, sprintf("%.2f (%.2f)", w$emmean[i], w$SE[i]),
         if (length(k)) fmt_ci(p$estimate[k], p$lower.CL[k], p$upper.CL[k]) else "",
         if (length(k)) sprintf("%.3f", p$p.value[k]) else "")
  })
  tb <- do.call(rtable, c(list(header = c("LS mean change (SE)", "Difference vs placebo (95% CI)", "p value")), rows))
  main_title(tb) <- "Table 3. ADAS-Cog(11) change from baseline at Week 24, MMRM (efficacy population)"
  main_footer(tb) <- paste("Unstructured covariance, REML, Kenward-Roger df; treatment, visit, treatment x visit,",
                           "baseline, baseline x visit and pooled site. Higher scores are worse.")
  tb
}

t_sens <- function(sens) {
  lab <- c("MMRM (primary)" = "MMRM, hypothetical (primary)", MAR = "Imputation, MAR",
           JR = "Imputation, jump to reference", CR = "Imputation, copy reference", "LOCF ANCOVA" = "LOCF ANCOVA (pilot method)")
  rows <- lapply(names(lab), function(a) {
    s <- sens[sens$analysis == a, ]
    cell <- function(arm) {
      r <- s[s$contrast == arm, ]
      sprintf("%s, p = %.3f", fmt_ci(r$est, r$lci, r$uci), r$pval)
    }
    rrow(lab[[a]], cell(ARMS[2]), cell(ARMS[3]))
  })
  tb <- do.call(rtable, c(list(header = c("Low dose vs placebo", "High dose vs placebo")), rows))
  main_title(tb) <- "Table 4. Sensitivity analyses, difference in ADAS-Cog(11) change at Week 24 (95% CI)"
  main_footer(tb) <- "Imputation: conditional mean imputation with jackknife standard errors (rbmi), Week 24 ANCOVA."
  tb
}

## T5 TEAEs by SOC and PT, PTs in >= 5% of any arm
t_ae <- function(adsl, adae) {
  saf <- adsl %>% filter(SAFFL == "Y") %>% mutate(TRTA = factor(TRT01A, ARMS))
  te <- adae %>% filter(SAFFL == "Y", TRTEMFL == "Y") %>% mutate(TRTA = factor(TRTA, ARMS))
  n_arm <- table(saf$TRTA)
  keep <- te %>% distinct(TRTA, AEDECOD, USUBJID) %>% count(TRTA, AEDECOD) %>%
    mutate(pct = n / as.numeric(n_arm[as.character(TRTA)])) %>% filter(pct >= 0.05) %>% pull(AEDECOD) %>% unique()
  te5 <- te %>% filter(AEDECOD %in% keep) %>% mutate(AEBODSYS = factor(AEBODSYS), AEDECOD = factor(AEDECOD))
  subj_pct <- function(df, .N_col, ...) {
    lv <- levels(df$AEDECOD)
    in_rows(.list = lapply(lv, function(l) {
      n <- length(unique(df$USUBJID[df$AEDECOD == l]))
      rcell(c(n, n / .N_col), format = "xx (xx.x%)")
    }), .names = lv)
  }
  any_te <- function(df, .N_col, ...) {
    n <- length(unique(df$USUBJID))
    in_rows("Patients with any TEAE" = rcell(c(n, n / .N_col), format = "xx (xx.x%)"))
  }
  over <- basic_table() %>%
    split_cols_by("TRTA") %>%
    add_colcounts() %>%
    analyze("USUBJID", any_te) %>%
    build_table(te, alt_counts_df = saf)
  socs <- basic_table() %>%
    split_cols_by("TRTA") %>%
    add_colcounts() %>%
    split_rows_by("AEBODSYS", split_fun = trim_levels_in_group("AEDECOD")) %>%
    analyze("AEDECOD", subj_pct) %>%
    build_table(te5, alt_counts_df = saf)
  tb <- rbind(over, socs)
  main_title(tb) <- "Table 5. Treatment-emergent adverse events in >= 5% of patients in any arm (safety population)"
  main_footer(tb) <- "Patients counted once per preferred term; columns by actual treatment."
  tb
}
