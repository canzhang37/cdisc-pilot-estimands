# AE signal screen, double FDR, Berry-Berry, onset
suppressPackageStartupMessages(library(dplyr))

wilson <- function(x, n, z = qnorm(0.975)) {
  p <- x / n
  mid <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(mid - half, mid + half)
}

## Newcombe (1998) hybrid score CI, p1 - p2
newcombe <- function(x1, n1, x2, n2) {
  p1 <- x1 / n1; p2 <- x2 / n2
  a <- wilson(x1, n1); b <- wilson(x2, n2)
  d <- p1 - p2
  c(est = d, lower = d - sqrt((p1 - a[1])^2 + (b[2] - p2)^2), upper = d + sqrt((a[2] - p1)^2 + (p2 - b[1])^2))
}

## patients with >= 1 TEAE per PT, arm vs placebo
pt_counts <- function(te, adsl, active, ref = "Placebo") {
  N <- adsl %>% filter(SAFFL == "Y") %>% count(TRT01A)
  n_of <- function(a) N$n[N$TRT01A == a]
  te %>%
    filter(TRTA %in% c(active, ref)) %>%
    distinct(USUBJID, TRTA, AEBODSYS, AEDECOD) %>%
    count(AEBODSYS, AEDECOD, TRTA) %>%
    tidyr::pivot_wider(names_from = TRTA, values_from = n, values_fill = 0) %>%
    rename(xt = all_of(active), xc = all_of(ref)) %>%
    mutate(comparison = active, nt = n_of(active), nc = n_of(ref))
}

screen <- function(cnt) {
  r <- t(mapply(newcombe, cnt$xt, cnt$nt, cnt$xc, cnt$nc))
  cnt %>%
    mutate(rd = r[, "est"], rd_lower = r[, "lower"], rd_upper = r[, "upper"],
           p_fisher = mapply(function(a, b, c, d) fisher.test(matrix(c(a, b - a, c, d - c), 2))$p.value, xt, nt, xc, nc),
           flag_crude = p_fisher < 0.05)
}

## Mehrotra & Heyse (2004): BH on SOCs, then PTs
double_fdr <- function(p, soc, a1 = 0.1, a2 = 0.1) {
  soc_p <- tapply(p, soc, min)
  keep <- names(soc_p)[p.adjust(soc_p, "BH") <= a1]
  sel <- rep(FALSE, length(p))
  for (s in keep) {
    i <- which(soc == s)
    sel[i] <- p.adjust(p[i], "BH") <= a2
  }
  sel
}

stan_input <- function(cnt) {
  soc <- as.integer(factor(cnt$AEBODSYS))
  list(J = nrow(cnt), B = max(soc), soc = soc, xc = cnt$xc, xt = cnt$xt, Nc = cnt$nc[1], Nt = cnt$nt[1])
}

fit_bb <- function(model, cnt, seed = 2027, adapt_delta = 0.95) {
  model$sample(data = stan_input(cnt), seed = seed, chains = 4, parallel_chains = 4, iter_warmup = 2000,
               iter_sampling = 2000, adapt_delta = adapt_delta, refresh = 0, show_messages = FALSE)
}

bb_summary <- function(fit) {
  pp <- fit$summary("p_raised", "mean")$mean
  or <- fit$summary("effect", ~ quantile(exp(.x), c(0.5, 0.025, 0.975)))
  d <- fit$diagnostic_summary(quiet = TRUE)
  rh <- fit$summary(c("mu_gamma0", "tau_gamma", "sigma_gamma", "mu_theta0", "tau_theta", "sigma_theta", "mu_pi",
                      "sigma_pi", "gamma", "theta", "lpi"), "rhat")$rhat
  list(p_raised = pp, or_med = or[[2]], or_lo = or[[3]], or_hi = or[[4]],
       divergences = sum(d$num_divergent), max_rhat = max(rh, na.rm = TRUE))
}

## first app-site TEAE, censored at last dose
onset_data <- function(adsl, adae) {
  first <- adae %>%
    filter(TRTEMFL == "Y", grepl("^APPLICATION SITE", AEDECOD)) %>%
    group_by(USUBJID) %>%
    summarise(ASTDT = min(ASTDT), .groups = "drop")
  adsl %>%
    filter(SAFFL == "Y") %>%
    select(USUBJID, TRTA = TRT01A, TRTSDT, TRTEDT, DCSREAS) %>%
    left_join(first, by = "USUBJID") %>%
    mutate(event = !is.na(ASTDT),
           time = if_else(event, as.numeric(ASTDT - TRTSDT) + 1, as.numeric(TRTEDT - TRTSDT) + 1),
           TRTA = factor(TRTA, ARMS))
}
