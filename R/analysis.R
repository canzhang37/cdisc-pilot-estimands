# MMRM, gls check, rbmi sensitivity, LOCF ANCOVA
suppressPackageStartupMessages({
  library(dplyr)
  library(mmrm)
  library(emmeans)
})

VIS <- c("Week 8", "Week 16", "Week 24")
COVS <- c("us", "toeph", "ar1h", "cs")

## efficacy population, patient x visit
eff_data <- function(adsl, adadas, arms) {
  pts <- adsl %>% filter(EFFFL == "Y") %>% select(USUBJID, TRTP = TRT01P, SITEGR1, EOSSTT, DCSREAS)
  obs <- adadas %>%
    filter(PARAMCD == "ACTOT", ANL01FL == "Y", AVISIT %in% VIS) %>%
    select(USUBJID, AVISIT, AVAL, BASE, CHG)
  base <- adadas %>% filter(PARAMCD == "ACTOT", ABLFL == "Y") %>% distinct(USUBJID, BASE = AVAL)
  tidyr::expand_grid(pts, AVISIT = VIS) %>%
    left_join(base, by = "USUBJID") %>%
    left_join(select(obs, -BASE), by = c("USUBJID", "AVISIT")) %>%
    mutate(TRTP = factor(TRTP, arms), AVISIT = factor(AVISIT, VIS), SITEGR1 = factor(SITEGR1),
           USUBJID = factor(USUBJID))
}

## primary MMRM, covariance fallback order
fit_mmrm <- function(d, method = "Kenward-Roger") {
  for (cv in COVS) {
    f <- as.formula(paste0("CHG ~ TRTP * AVISIT + BASE * AVISIT + SITEGR1 + ", cv, "(AVISIT | USUBJID)"))
    fit <- tryCatch(mmrm(f, data = d, method = method), error = function(e) NULL)
    if (!is.null(fit) && isTRUE(component(fit, "convergence") == 0)) return(list(fit = fit, cov = cv))
  }
  stop("no covariance structure converged")
}

lsm <- function(fit, ...) emmeans(fit, ~ TRTP | AVISIT, weights = "proportional", ...)

diffs <- function(em) {
  ct <- contrast(em, "trt.vs.ctrl", ref = 1)
  # fixed sequence, so no Dunnett adjustment
  s <- as.data.frame(summary(ct, infer = TRUE, adjust = "none"))
  s[, c("contrast", "AVISIT", "estimate", "SE", "df", "lower.CL", "upper.CL", "p.value")]
}

fit_gls <- function(d) {
  dd <- filter(d, !is.na(CHG)) %>% mutate(vi = as.integer(AVISIT))
  fit <- nlme::gls(CHG ~ TRTP * AVISIT + BASE * AVISIT + SITEGR1, data = dd, method = "REML",
                   correlation = nlme::corSymm(form = ~ vi | USUBJID), weights = nlme::varIdent(form = ~ 1 | AVISIT))
  list(fit = fit, data = dd)
}

## ICE: first visit after the last observed
ice_data <- function(d, strategy) {
  d %>%
    group_by(USUBJID) %>%
    summarise(last = suppressWarnings(max(as.integer(AVISIT)[!is.na(CHG)])), .groups = "drop") %>%
    filter(is.finite(last), last < length(VIS)) %>%
    transmute(USUBJID, AVISIT = factor(VIS[last + 1], VIS), strategy = strategy)
}

## Week 24 ANCOVA per imputed set
ancova24 <- function(data, ...) {
  w <- filter(data, AVISIT == "Week 24")
  fit <- lm(CHG ~ TRTP + BASE + SITEGR1, data = w)
  em <- emmeans(fit, ~ TRTP, weights = "proportional")
  ct <- as.data.frame(summary(contrast(em, "trt.vs.ctrl", ref = 1)))
  out <- lapply(seq_len(nrow(ct)), function(i) list(est = ct$estimate[i], se = ct$SE[i], df = ct$df[i]))
  setNames(out, paste0("trt_", seq_len(nrow(ct)) + 1))
}

## LOCF at Week 24, ANCOVA
locf_ancova <- function(d) {
  w <- d %>%
    filter(!is.na(CHG)) %>%
    group_by(USUBJID) %>%
    slice_max(as.integer(AVISIT), n = 1) %>%
    ungroup()
  fit <- lm(CHG ~ TRTP + BASE + SITEGR1, data = w)
  em <- emmeans(fit, ~ TRTP, weights = "proportional")
  list(n = nrow(w), em = as.data.frame(em),
       diff = as.data.frame(summary(contrast(em, "trt.vs.ctrl", ref = 1), infer = TRUE, adjust = "none")))
}
