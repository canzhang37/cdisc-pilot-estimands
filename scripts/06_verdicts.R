# P1-P6 as pre-registered
rd <- function(f) read.csv(file.path("results", f))
sl <- rd("qc_adsl.csv"); tot <- rd("qc_adas_total.csv"); chk <- rd("check_mmrm_gls.csv"); s <- rd("sensitivity.csv")
hi <- "Xanomeline High Dose"
g <- function(a) s$est[s$analysis == a & s$contrast == hi]
v <- data.frame(
  id = c("P1", "P1", "P2", "P3", "P4", "P5", "P6"),
  check = c("lowest ADSL agreement, tolerance rule (>= 95%)", "lowest ADSL agreement, R rounding rule (>= 95%)",
            "recomputed = collected total when 11 items present (all)", "mmrm vs gls: max |d est| (< 1e-4), max |d se| (< 1e-3)",
            "|MAR - MMRM|, high dose (<= 0.3)", "J2R shrinkage toward 0, high dose (>= 20%)",
            "LOCF minus MMRM, high dose (< 0, more favorable)"),
  value = c(sprintf("%.1f%% (%s)", min(sl$pct), sl$variable[which.min(sl$pct)]),
            sprintf("%.1f%% (%s)", min(sl$pct_r_round), sl$variable[which.min(sl$pct_r_round)]),
            sprintf("%d of %d", tot$equal[1], tot$records[1]),
            sprintf("%.1e, %.1e", max(chk$d_est), max(chk$d_se)),
            sprintf("%.2f", abs(g("MAR") - g("MMRM (primary)"))),
            sprintf("%.0f%%", 100 * (1 - g("JR") / g("MMRM (primary)"))),
            sprintf("%.2f", g("LOCF ANCOVA") - g("MMRM (primary)"))),
  met = c(min(sl$pct) >= 95, min(sl$pct_r_round) >= 95, tot$equal[1] == tot$records[1],
          max(chk$d_est) < 1e-4 & max(chk$d_se) < 1e-3, abs(g("MAR") - g("MMRM (primary)")) <= 0.3,
          1 - g("JR") / g("MMRM (primary)") >= 0.20, g("LOCF ANCOVA") - g("MMRM (primary)") < 0)
)
write.csv(v, "results/predictions.csv", row.names = FALSE)
print(v, right = FALSE)
