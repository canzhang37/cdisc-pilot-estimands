# MMRM, software check, sensitivity
source("R/adam.R")
source("R/analysis.R")
suppressPackageStartupMessages(library(rbmi))
adsl <- readRDS("adam/adsl.rds")
adadas <- readRDS("adam/adadas.rds")
d <- eff_data(adsl, adadas, ARMS)
out <- function(x, f) write.csv(x, file.path("results", f), row.names = FALSE)

## observed data and dropout
obs <- d %>%
  group_by(TRTP, AVISIT) %>%
  summarise(n = n(), observed = sum(!is.na(CHG)), mean_chg = mean(CHG, na.rm = TRUE), sd_chg = sd(CHG, na.rm = TRUE), .groups = "drop")
out(obs, "observed_by_visit.csv")
ice <- ice_data(d, "MAR")
drop <- d %>% distinct(USUBJID, TRTP) %>% mutate(dropout = USUBJID %in% ice$USUBJID) %>%
  count(TRTP, dropout) %>% tidyr::pivot_wider(names_from = dropout, values_from = n, names_prefix = "dropout_")
out(drop, "dropout_by_arm.csv")
print(obs, n = 20); print(drop)

## primary
pm <- fit_mmrm(d, "Kenward-Roger")
em <- lsm(pm$fit)
lsmeans <- as.data.frame(summary(em, infer = TRUE))
prim <- diffs(em)
cat("\ncovariance:", pm$cov, "\n"); print(prim)
out(cbind(covariance = pm$cov, lsmeans), "lsmeans_by_visit.csv")
out(cbind(covariance = pm$cov, prim), "primary_mmrm.csv")

## software check: mmrm (Satterthwaite) vs gls
ms <- fit_mmrm(d, "Satterthwaite")
g <- fit_gls(d)
a <- diffs(lsm(ms$fit))
b <- diffs(emmeans(g$fit, ~ TRTP | AVISIT, weights = "proportional", mode = "df.error", data = g$data))
chk <- data.frame(contrast = a$contrast, AVISIT = a$AVISIT, est_mmrm = a$estimate, est_gls = b$estimate,
                  se_mmrm = a$SE, se_gls = b$SE)
chk$d_est <- abs(chk$est_mmrm - chk$est_gls); chk$d_se <- abs(chk$se_mmrm - chk$se_gls)
out(chk, "check_mmrm_gls.csv")
cat("\nmmrm vs gls: max |d est|", signif(max(chk$d_est), 3), " max |d se|", signif(max(chk$d_se), 3), "\n")

## rbmi condmean + jackknife, per strategy
vars <- set_vars(subjid = "USUBJID", visit = "AVISIT", group = "TRTP", outcome = "CHG",
                 covariates = c("TRTP*AVISIT", "BASE*AVISIT", "SITEGR1"), strategy = "strategy")
rd <- as.data.frame(select(d, USUBJID, TRTP, AVISIT, CHG, BASE, SITEGR1))
refs <- setNames(rep("Placebo", 3), ARMS)
imps <- lapply(setNames(nm = c("MAR", "JR", "CR")), function(s) {
  dr <- draws(data = rd, data_ice = as.data.frame(mutate(ice, strategy = s)), vars = vars,
              method = method_condmean(type = "jackknife"), quiet = TRUE)
  impute(dr, references = refs)
})
run <- function(imp, delta_arm = NULL, delta = 0) {
  dl <- NULL
  if (!is.null(delta_arm)) {
    dl <- delta_template(imp)
    dl$delta <- ifelse(dl$TRTP == delta_arm & dl$is_post_ice & dl$is_missing, delta, 0)
  }
  po <- as.data.frame(pool(analyse(imp, ancova24, delta = dl)))
  po$contrast <- c(trt_2 = ARMS[2], trt_3 = ARMS[3])[po$parameter]
  po
}
sens <- bind_rows(lapply(names(imps), function(s) cbind(analysis = s, run(imps[[s]]))))

## check: no site, condmean MAR = MMRM
v0 <- set_vars(subjid = "USUBJID", visit = "AVISIT", group = "TRTP", outcome = "CHG",
               covariates = c("TRTP*AVISIT", "BASE*AVISIT"), strategy = "strategy")
i0 <- impute(draws(rd, as.data.frame(ice), v0, method_condmean(type = "jackknife"), quiet = TRUE), references = refs)
a0 <- function(data, ...) {
  f <- lm(CHG ~ TRTP + BASE, data = filter(data, AVISIT == "Week 24"))
  cf <- summary(f)$coefficients
  list(trt_2 = list(est = cf[2, 1], se = cf[2, 2], df = f$df.residual), trt_3 = list(est = cf[3, 1], se = cf[3, 2], df = f$df.residual))
}
p0 <- as.data.frame(pool(analyse(i0, a0)))
f0 <- mmrm(CHG ~ TRTP * AVISIT + BASE * AVISIT + us(AVISIT | USUBJID), data = d, method = "Kenward-Roger")
e0 <- diffs(emmeans(f0, ~ TRTP | AVISIT)); e0 <- e0[e0$AVISIT == "Week 24", ]
cm <- data.frame(contrast = sub(" - Placebo", "", e0$contrast), mmrm_est = e0$estimate, rbmi_est = p0$est,
                 mmrm_se = e0$SE, rbmi_se = p0$se)
out(cm, "check_condmean_mmrm.csv")
cat("
condmean vs MMRM (no site): max |d est|", signif(max(abs(cm$mmrm_est - cm$rbmi_est)), 3), "
")

## LOCF ANCOVA, the pilot's approach
lc <- locf_ancova(d)
locf <- data.frame(analysis = "LOCF ANCOVA", parameter = NA, est = lc$diff$estimate, se = lc$diff$SE,
                   lci = lc$diff$lower.CL, uci = lc$diff$upper.CL, pval = lc$diff$p.value,
                   contrast = sub(" - Placebo", "", lc$diff$contrast))
w24 <- prim[prim$AVISIT == "Week 24", ]
mm <- data.frame(analysis = "MMRM (primary)", parameter = NA, est = w24$estimate, se = w24$SE, lci = w24$lower.CL,
                 uci = w24$upper.CL, pval = w24$p.value, contrast = sub(" - Placebo", "", w24$contrast))
sens_all <- bind_rows(mm, sens, locf)
out(sens_all, "sensitivity.csv")
print(sens_all)

## tipping point, significant comparisons only
tp <- NULL
for (arm in ARMS[2:3]) {
  if (w24$p.value[sub(" - Placebo", "", w24$contrast) == arm] >= 0.05) next
  for (dv in seq(0, 10, by = 0.5)) {
    r <- run(imps$MAR, delta_arm = arm, delta = dv)
    r <- r[r$contrast == arm, ]
    tp <- rbind(tp, data.frame(contrast = arm, delta = dv, est = r$est, lci = r$lci, uci = r$uci, pval = r$pval))
  }
}
if (!is.null(tp)) out(tp, "tipping_point.csv")
print(tp)
