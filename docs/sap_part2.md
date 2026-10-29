# Plan, part 2: safety signal detection and subgroup analyses

This plan extends [sap.md](sap.md). It was written after the part 1 results, including table T5 of treatment-emergent adverse events (TEAEs) by preferred term, had been seen, and after reading the structure of ADAE and ADSL (230 preferred terms in 23 system organ classes; baseline MMSE median 19.5 in the efficacy population). No signal-detection or subgroup result had been computed. Predictions that lean on T5 say so. Changes made after results are seen go in the Amendments section and are labeled in the readout.

## 1. Why

A trial with 230 adverse event terms will show some differences between arms by chance, and a trial with a null primary result invites a search for the patients who did benefit. Both questions need the same discipline: a method that accounts for how many comparisons were made, decided before looking.

## 2. Safety signal detection

Safety population (SAFFL = Y) by actual treatment, TEAEs only (TRTEMFL = Y). For each preferred term (PT), the number of patients with at least one TEAE of that term. Comparisons: high dose versus placebo (main) and low dose versus placebo. A PT enters a comparison when at least one patient in either of its two arms had it.

- A1. Crude screen. Risk difference with the Newcombe hybrid score 95% interval (Newcombe 1998, method 10) and Fisher's exact two-sided p value. Flag: p < 0.05.
- A2. Double false discovery rate (Mehrotra and Heyse 2004). The p value of a system organ class (SOC) is the smallest Fisher p value among its PTs. Benjamini-Hochberg at 0.1 across SOCs; then, within each selected SOC, Benjamini-Hochberg at 0.1 across its PTs. Flag: PT selected in both steps.
- A3. Bayesian hierarchical mixture model (after Berry and Berry 2004), fitted in Stan. For PT j in SOC b, placebo and active counts are binomial with logit rates gamma_bj and gamma_bj + theta_bj. gamma_bj ~ N(mu_gamma_b, sigma_gamma), mu_gamma_b ~ N(mu_gamma0, tau_gamma). The log odds ratio theta_bj is exactly 0 with probability pi_b and N(mu_theta_b, sigma_theta) otherwise, with mu_theta_b ~ N(mu_theta0, tau_theta) and logit(pi_b) ~ N(mu_pi, sigma_pi). Priors: mu_gamma0 ~ N(-3, 2), mu_theta0 ~ N(0, 1), mu_pi ~ N(1.4, 1) (prior probability of no effect about 0.8), half-normal(0, 1) for every scale except tau_theta ~ half-normal(0, 0.5). The point mass is handled by summing over the indicator in the likelihood; the posterior probability that a PT has a raised risk, P(theta_bj > 0 and not in the point mass | data), is averaged over draws from the indicator's conditional probability. Four chains, 2,000 warmup and 2,000 sampling iterations each, adapt_delta 0.95. Flag: posterior probability at least 0.8.
- A4. Leading signal over time. Application-site reactions are grouped as all PTs whose name starts with "APPLICATION SITE". Time to the first such TEAE from the first dose (ASTDT - TRTSDT + 1); patients without one are censored at their last dose (TRTEDT - TRTSDT + 1). Kaplan-Meier by arm, cumulative incidence at Weeks 2, 4, 8 and 24, median onset among patients with an event, and a Cox hazard ratio for each dose versus placebo. Discontinuation for adverse events by arm (DCSREAS), and how many of those patients had an application-site reaction first.

## 3. Subgroup analyses of the primary endpoint

Efficacy population, change from baseline in ADAS-Cog(11), the primary MMRM of sap.md with a subgroup factor S added as S x visit and S x treatment x visit (same covariance fallback, Kenward-Roger). Subgroups:
- age group (AGEGR1: <65, 65-80, >80);
- sex;
- baseline disease severity by screening MMSE: moderate (10-19) versus mild (20-24), which is also the median split.

For each level, the Week 24 difference from placebo for each dose with its 95% interval; for each subgroup factor and dose, the joint test that the Week 24 differences are equal across levels (2 degrees of freedom for age group, 1 otherwise). Six interaction tests in all, no multiplicity adjustment, read as exploratory (EMA guideline on subgroup analyses in confirmatory trials, 2019).

## 4. Checks of the code

- V1. The Newcombe interval reproduces the worked example of Newcombe (1998): 56/70 versus 48/80 gives 0.0524 to 0.3339.
- V2. With every PT placed in one SOC, the double FDR selects the same PTs as Benjamini-Hochberg at 0.1.
- V3. On simulated data shaped like this trial (60 null PTs and 3 PTs with odds ratio 4 at a 10% placebo rate, 3 SOCs, 80 patients per arm), the Bayesian model gives posterior probability at least 0.8 to at least 2 of the 3 planted signals and to at most 3 of the 60 null PTs.

## 5. Predictions

- Q1. For high dose versus placebo, application site pruritus and application site erythema are flagged by A1, A2 and A3. Informed by T5.
- Q2. For high dose versus placebo, A2 flags fewer PTs than A1.
- Q3. For high dose versus placebo, A3 flags fewer PTs than A1.
- Q4. Among high-dose patients with an application-site reaction, the median time to the first one is at most 28 days.
- Q5. None of the six subgroup interaction tests has p < 0.01, and at most one has p < 0.05.
- Q6. The Stan fit of A3 has no divergent transitions and every R-hat below 1.01. If not, it is rerun once with adapt_delta 0.99, and the rerun is reported as an amendment.

## 6. Tables and figures

- T6 safety signals: PTs flagged by any method in either comparison, with counts, risk difference, Fisher p and each method's flag.
- T7 subgroup results at Week 24.
- F3 risk differences with intervals for the flagged PTs.
- F4 Kaplan-Meier curves of time to the first application-site reaction.
- F5 subgroup forest plot.

## Amendments

None yet.
