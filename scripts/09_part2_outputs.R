# tables T6-T7, figures F3-F5
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(rtables)
  library(survival)
})
source("R/adam.R")
source("R/safety.R")
COL <- setNames(c("#2a78d6", "#eb6834", "#1baf7a"), ARMS)
th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#e6e6e3", linewidth = 0.3),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA),
        legend.position = "top", legend.title = element_blank(), axis.text = element_text(colour = "#5f5e5a"),
        strip.text = element_text(face = "bold", hjust = 0))
short <- function(a) sub("Xanomeline ", "", a)
pct <- function(x, n) sprintf("%d/%d (%.1f%%)", x, n, 100 * x / n)

## T6
scr <- read.csv("results/safety_screen.csv")
fl <- scr %>%
  filter(flag_crude | flag_dfdr | flag_bayes) %>%
  mutate(comparison = factor(comparison, ARMS[3:2])) %>%
  arrange(comparison, p_fisher)
t6 <- data.frame(
  Comparison = paste(short(as.character(fl$comparison)), "vs placebo"),
  `Preferred term` = tools::toTitleCase(tolower(fl$AEDECOD)),
  Active = pct(fl$xt, fl$nt), Placebo = pct(fl$xc, fl$nc),
  `Risk difference (95% CI)` = sprintf("%.1f (%.1f, %.1f)", 100 * fl$rd, 100 * fl$rd_lower, 100 * fl$rd_upper),
  `Fisher p` = sprintf("%.4f", fl$p_fisher),
  `Double FDR` = ifelse(fl$flag_dfdr, "flag", ""),
  `Bayes P(raised)` = sprintf("%.2f%s", fl$p_raised, ifelse(fl$flag_bayes, " flag", "")),
  check.names = FALSE)
cat("Table 6. Preferred terms flagged by any method (safety population, TEAEs)\n",
    "Crude screen flags Fisher p < 0.05; double FDR at 0.1/0.1 by system organ class; Bayesian hierarchical\n",
    "mixture model flags posterior probability of a raised risk >= 0.8. Risk differences in percentage points.\n\n",
    file = "tlf/t6_signals.txt")
cat(export_as_txt(df_to_tt(t6), paginate = FALSE), file = "tlf/t6_signals.txt", append = TRUE)

## T7
sub <- read.csv("results/subgroups.csv")
itx <- read.csv("results/subgroup_interactions.csv")
lab <- c(age_group = "Age group", sex = "Sex", severity = "Baseline severity")
t7 <- sub %>%
  left_join(itx %>% select(subgroup, contrast, p_interaction), by = c("subgroup", "contrast")) %>%
  transmute(Subgroup = lab[subgroup], Level = level, Comparison = paste(short(sub(" - Placebo", "", contrast)), "vs placebo"),
            `n (active / placebo)` = paste(n_active, n_placebo, sep = " / "),
            `Difference (95% CI)` = sprintf("%.2f (%.2f, %.2f)", est, lower, upper),
            `Interaction p` = sprintf("%.3f", p_interaction))
cat("Table 7. ADAS-Cog(11) change from baseline at Week 24 by subgroup (efficacy population)\n",
    "MMRM with subgroup x treatment x visit; unadjusted 95% intervals; interaction p is the joint test that the\n",
    "Week 24 difference is the same in every level of the subgroup. Exploratory.\n\n", file = "tlf/t7_subgroups.txt")
cat(export_as_txt(df_to_tt(t7), paginate = FALSE), file = "tlf/t7_subgroups.txt", append = TRUE)

## F3
f3 <- fl %>%
  mutate(n_methods = flag_crude + flag_dfdr + flag_bayes,
         methods = factor(c("1 method", "2 methods", "3 methods")[n_methods], c("1 method", "2 methods", "3 methods")),
         term = tools::toTitleCase(tolower(AEDECOD)),
         panel = paste(short(as.character(comparison)), "vs placebo"))
f3$term <- factor(f3$term, unique(f3$term[order(f3$rd)]))
p3 <- ggplot(f3, aes(100 * rd, term, colour = methods)) +
  geom_vline(xintercept = 0, colour = "#b4b2a9", linewidth = 0.4) +
  geom_errorbar(aes(xmin = 100 * rd_lower, xmax = 100 * rd_upper), width = 0, orientation = "y", linewidth = 0.6) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c("1 method" = "#b4b2a9", "2 methods" = "#eb9a34", "3 methods" = "#c4351e"), drop = FALSE) +
  facet_wrap(~panel) +
  labs(x = "Risk difference, percentage points (Newcombe 95% CI)", y = NULL,
       caption = "Terms flagged by at least one of: Fisher p < 0.05, double FDR, Bayesian P(raised) >= 0.8") + th
ggsave("figures/safety_signals.png", p3, width = 9, height = 5.2, dpi = 150)

## F4
od <- onset_data(readRDS("adam/adsl.rds"), readRDS("adam/adae.rds"))
km <- survfit(Surv(time, event) ~ TRTA, data = od)
kd <- data.frame(time = km$time, cum = 1 - km$surv, arm = rep(sub("TRTA=", "", names(km$strata)), km$strata))
kd <- bind_rows(data.frame(time = 0, cum = 0, arm = ARMS), kd) %>% mutate(arm = factor(arm, ARMS))
p4 <- ggplot(kd, aes(time / 7, 100 * cum, colour = arm)) +
  geom_step(linewidth = 0.9) +
  scale_colour_manual(values = COL) +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24, 28)) +
  labs(x = "Weeks since first dose", y = "Patients with an application-site reaction (%)",
       caption = "Kaplan-Meier, safety population; censored at last dose") +
  th + theme(panel.grid.major.y = element_line(colour = "#e6e6e3", linewidth = 0.3))
ggsave("figures/onset_application_site.png", p4, width = 7, height = 4.3, dpi = 150)

## F5
prim <- read.csv("results/primary_mmrm.csv") %>% filter(AVISIT == "Week 24")
f5 <- sub %>%
  left_join(itx %>% select(subgroup, contrast, p_interaction), by = c("subgroup", "contrast")) %>%
  mutate(panel = paste(short(sub(" - Placebo", "", contrast)), "vs placebo"),
         row = paste0(lab[subgroup], ": ", level),
         row = factor(row, rev(unique(row))))
ov <- prim %>% transmute(panel = paste(short(sub(" - Placebo", "", contrast)), "vs placebo"), overall = estimate)
p5 <- ggplot(f5, aes(est, row)) +
  geom_vline(xintercept = 0, colour = "#b4b2a9", linewidth = 0.4) +
  geom_vline(data = ov, aes(xintercept = overall), colour = "#2a78d6", linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, orientation = "y", linewidth = 0.6, colour = "#2a78d6") +
  geom_point(size = 2.4, colour = "#2a78d6") +
  facet_wrap(~panel) +
  labs(x = "Difference in ADAS-Cog(11) change at Week 24 (95% CI); below 0 favors xanomeline", y = NULL,
       caption = "Dashed line: overall MMRM estimate. All six interaction p values >= 0.22.") + th
ggsave("figures/subgroup_forest.png", p5, width = 9, height = 4.4, dpi = 150)
cat("T6, T7, F3-F5 written\n")
