# Q1-Q6 and V3 as pre-registered (V1, V2 in tests)
rd <- function(f) read.csv(file.path("results", f))
scr <- rd("safety_screen.csv"); dg <- rd("stan_diagnostics.csv"); med <- rd("onset_median.csv")
itx <- rd("subgroup_interactions.csv"); v3 <- rd("v3_simulation.csv")
hi <- subset(scr, comparison == "Xanomeline High Dose")
all3 <- function(pt) with(hi[hi$AEDECOD == pt, ], flag_crude & flag_dfdr & flag_bayes)
m_hi <- med$median_day[med$TRTA == "Xanomeline High Dose"]
v <- data.frame(
  id = c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "V3"),
  check = c("application site pruritus and erythema flagged by A1, A2, A3 (high dose)",
            "double FDR flags fewer PTs than the crude screen (high dose)",
            "Bayesian model flags fewer PTs than the crude screen (high dose)",
            "median days to first application-site reaction, high dose (<= 28)",
            "subgroup interactions: none p < 0.01, at most one p < 0.05",
            "Stan: no divergences, max R-hat < 1.01 (both comparisons)",
            "simulation: >= 2 of 3 planted and <= 3 of 60 null PTs flagged"),
  value = c(sprintf("pruritus %s, erythema %s", all3("APPLICATION SITE PRURITUS"), all3("APPLICATION SITE ERYTHEMA")),
            sprintf("%d vs %d", sum(hi$flag_dfdr), sum(hi$flag_crude)),
            sprintf("%d vs %d", sum(hi$flag_bayes), sum(hi$flag_crude)),
            sprintf("%.1f", m_hi),
            sprintf("min p %.2f, %d below 0.05", min(itx$p_interaction), sum(itx$p_interaction < 0.05)),
            sprintf("%d divergences, max R-hat %.4f", sum(dg$divergences), max(dg$max_rhat)),
            sprintf("%d of 3, %d of 60", sum(v3$p_raised[v3$planted] >= 0.8), sum(v3$p_raised[!v3$planted] >= 0.8))),
  met = c(all3("APPLICATION SITE PRURITUS") & all3("APPLICATION SITE ERYTHEMA"),
          sum(hi$flag_dfdr) < sum(hi$flag_crude), sum(hi$flag_bayes) < sum(hi$flag_crude), m_hi <= 28,
          min(itx$p_interaction) >= 0.01 & sum(itx$p_interaction < 0.05) <= 1,
          sum(dg$divergences) == 0 & max(dg$max_rhat) < 1.01,
          sum(v3$p_raised[v3$planted] >= 0.8) >= 2 & sum(v3$p_raised[!v3$planted] >= 0.8) <= 3)
)
write.csv(v, "results/predictions_part2.csv", row.names = FALSE)
print(v, right = FALSE)
