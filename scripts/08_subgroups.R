# subgroup MMRM: subgroup x trt x visit
source("R/adam.R")
source("R/analysis.R")
adsl <- readRDS("adam/adsl.rds")
adadas <- readRDS("adam/adadas.rds")
d0 <- eff_data(adsl, adadas, ARMS) %>%
  left_join(adsl %>% select(USUBJID, AGEGR1, SEX, MMSETOT), by = "USUBJID") %>%
  mutate(age_group = factor(AGEGR1, c("<65", "65-80", ">80")),
         sex = factor(SEX, c("F", "M"), c("Female", "Male")),
         severity = factor(if_else(MMSETOT <= 19, "Moderate (MMSE 10-19)", "Mild (MMSE 20-24)"),
                           c("Moderate (MMSE 10-19)", "Mild (MMSE 20-24)")))

fit_sub <- function(d) {
  for (cv in COVS) {
    f <- as.formula(paste0("CHG ~ TRTP * AVISIT * S + BASE * AVISIT + SITEGR1 + ", cv, "(AVISIT | USUBJID)"))
    fit <- tryCatch(mmrm(f, data = d, method = "Kenward-Roger"), error = function(e) NULL)
    if (!is.null(fit) && isTRUE(component(fit, "convergence") == 0)) return(list(fit = fit, cov = cv))
  }
  stop("no covariance structure converged")
}

res <- list(); tests <- list()
for (v in c("age_group", "sex", "severity")) {
  d <- mutate(d0, S = .data[[v]])
  m <- fit_sub(d)
  em <- emmeans(m$fit, ~ TRTP | S, at = list(AVISIT = "Week 24"), weights = "proportional")
  dd <- contrast(em, "trt.vs.ctrl", ref = 1)
  n <- d %>% filter(AVISIT == "Week 24") %>% count(S, TRTP) %>%
    tidyr::pivot_wider(names_from = TRTP, values_from = n)
  s <- as.data.frame(summary(dd, infer = TRUE, adjust = "none")) %>% left_join(n, by = "S")
  res[[v]] <- data.frame(subgroup = v, level = s$S, contrast = s$contrast, n_placebo = s$Placebo,
                         n_active = ifelse(grepl("High", s$contrast), s$`Xanomeline High Dose`, s$`Xanomeline Low Dose`),
                         est = s$estimate, se = s$SE, df = s$df, lower = s$lower.CL, upper = s$upper.CL,
                         covariance = m$cov)
  ii <- contrast(dd, "trt.vs.ctrl", by = "contrast", simple = "S")
  jt <- as.data.frame(test(ii, joint = TRUE, by = "contrast"))
  tests[[v]] <- data.frame(subgroup = v, contrast = jt$contrast, df1 = jt$df1, df2 = jt$df2, F = jt$F.ratio,
                           p_interaction = jt$p.value)
}
sub <- bind_rows(res); itx <- bind_rows(tests)
write.csv(sub, "results/subgroups.csv", row.names = FALSE)
write.csv(itx, "results/subgroup_interactions.csv", row.names = FALSE)
print(sub, digits = 3); print(itx, digits = 3)
