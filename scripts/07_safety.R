# safety signals, app-site onset
source("R/adam.R")
source("R/safety.R")
suppressPackageStartupMessages({
  library(cmdstanr)
  library(survival)
})
adsl <- readRDS("adam/adsl.rds")
adae <- readRDS("adam/adae.rds")
te <- filter(adae, SAFFL == "Y", TRTEMFL == "Y")
out <- function(x, f) write.csv(x, file.path("results", f), row.names = FALSE)
bb <- cmdstan_model("stan/berry_berry.stan", quiet = TRUE)

## V3: planted signals, trial-shaped data
set.seed(2027)
J <- 63; n <- 80
sim <- data.frame(AEBODSYS = rep(c("S1", "S2", "S3"), each = 21), AEDECOD = paste0("PT", 1:J),
                  planted = seq_len(J) %in% c(1, 22, 43))
pc <- plogis(rnorm(J, qlogis(0.05), 0.8))
pc[sim$planted] <- 0.10
pt <- plogis(qlogis(pc) + log(4) * sim$planted)
sim <- mutate(sim, xc = rbinom(J, n, pc), xt = rbinom(J, n, pt), nc = n, nt = n)
v3fit <- fit_bb(bb, sim)
v3 <- bb_summary(v3fit)
sim$p_raised <- v3$p_raised
out(sim, "v3_simulation.csv")
cat(sprintf("V3: planted flagged %d of 3, null flagged %d of 60, divergences %d, max R-hat %.3f\n",
            sum(sim$p_raised[sim$planted] >= 0.8), sum(sim$p_raised[!sim$planted] >= 0.8), v3$divergences, v3$max_rhat))

## A1-A3 for each dose vs placebo
diag <- list()
scr <- lapply(ARMS[3:2], function(a) {
  s <- screen(pt_counts(te, adsl, a))
  s$flag_dfdr <- double_fdr(s$p_fisher, s$AEBODSYS)
  fit <- fit_bb(bb, s)
  b <- bb_summary(fit)
  diag[[a]] <<- data.frame(comparison = a, pts = nrow(s), divergences = b$divergences, max_rhat = b$max_rhat)
  s %>% mutate(p_raised = b$p_raised, or_median = b$or_med, or_lower = b$or_lo, or_upper = b$or_hi,
               flag_bayes = p_raised >= 0.8)
}) %>% bind_rows()
out(scr, "safety_screen.csv")
out(bind_rows(diag), "stan_diagnostics.csv")
print(bind_rows(diag))
flagged <- scr %>% filter(flag_crude | flag_dfdr | flag_bayes) %>%
  arrange(comparison, p_fisher) %>%
  select(comparison, AEDECOD, xt, nt, xc, nc, rd, p_fisher, flag_crude, flag_dfdr, p_raised, flag_bayes)
print(as.data.frame(flagged), digits = 3)
cat("\nflags per method:\n")
print(scr %>% group_by(comparison) %>% summarise(pts = n(), crude = sum(flag_crude), double_fdr = sum(flag_dfdr), bayes = sum(flag_bayes)))

## A4: time to the first application-site reaction
od <- onset_data(adsl, adae)
km <- survfit(Surv(time, event) ~ TRTA, data = od)
s <- summary(km, times = c(14, 28, 56, 168), extend = TRUE)
cum <- data.frame(arm = sub("TRTA=", "", s$strata), day = s$time, cum_inc = 1 - s$surv, at_risk = s$n.risk)
cox <- coxph(Surv(time, event) ~ TRTA, data = od)
hr <- data.frame(arm = sub("TRTA", "", names(coef(cox))), hr = exp(coef(cox)), lower = exp(confint(cox))[, 1],
                 upper = exp(confint(cox))[, 2])
med <- od %>% filter(event) %>% group_by(TRTA) %>%
  summarise(patients = n(), median_day = median(time), q1 = quantile(time, 0.25), q3 = quantile(time, 0.75), .groups = "drop")
out(cum, "onset_cuminc.csv"); out(hr, "onset_hr.csv"); out(med, "onset_median.csv")
print(cum); print(hr); print(med)

dc <- od %>%
  group_by(TRTA) %>%
  summarise(n = n(), disc_ae = sum(DCSREAS %in% "ADVERSE EVENT"),
            disc_ae_after_app_site = sum(DCSREAS %in% "ADVERSE EVENT" & event & ASTDT <= TRTEDT), .groups = "drop")
out(dc, "ae_discontinuation.csv")
print(dc)
