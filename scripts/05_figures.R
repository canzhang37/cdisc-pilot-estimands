# F1 LS means over time, F2 forest of all analyses
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("R/adam.R")
COL <- setNames(c("#2a78d6", "#eb6834", "#1baf7a"), ARMS)
th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e6e3", linewidth = 0.3),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA),
        legend.position = "top", legend.title = element_blank(), axis.text = element_text(colour = "#5f5e5a"))

## F1
lsm <- read.csv("results/lsmeans_by_visit.csv") %>%
  mutate(TRTP = factor(TRTP, ARMS), wk = c("Week 8" = 8, "Week 16" = 16, "Week 24" = 24)[AVISIT])
base <- data.frame(TRTP = factor(ARMS, ARMS), wk = 0, emmean = 0, lower.CL = 0, upper.CL = 0)
dd <- bind_rows(base, lsm) %>% mutate(x = wk + (as.integer(TRTP) - 2) * 0.5)
p1 <- ggplot(dd, aes(x, emmean, colour = TRTP)) +
  geom_hline(yintercept = 0, colour = "#b4b2a9", linewidth = 0.4) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0, linewidth = 0.6) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
  scale_colour_manual(values = COL) +
  scale_x_continuous(breaks = c(0, 8, 16, 24), labels = c("Baseline", "Week 8", "Week 16", "Week 24")) +
  labs(x = NULL, y = "LS mean change in ADAS-Cog(11), 95% CI",
       caption = "MMRM, efficacy population; higher is worse") + th
ggsave("figures/lsmeans_over_time.png", p1, width = 7, height = 4.3, dpi = 150)

## F2
sens <- read.csv("results/sensitivity.csv")
lab <- c("MMRM (primary)" = "MMRM (primary)", MAR = "Imputation, MAR", JR = "Imputation, J2R",
         CR = "Imputation, CR", "LOCF ANCOVA" = "LOCF ANCOVA")
sens <- sens %>% mutate(analysis = factor(lab[analysis], rev(lab)),
                        contrast = factor(paste(contrast, "vs placebo"), paste(ARMS[2:3], "vs placebo")))
p2 <- ggplot(sens, aes(est, analysis)) +
  geom_vline(xintercept = 0, colour = "#b4b2a9", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0, linewidth = 0.7, colour = "#2a78d6") +
  geom_point(size = 2.6, colour = "#2a78d6") +
  facet_wrap(~contrast) +
  labs(x = "Difference in ADAS-Cog(11) change at Week 24 (95% CI); below 0 favors xanomeline", y = NULL) +
  th + theme(strip.text = element_text(face = "bold", hjust = 0))
ggsave("figures/sensitivity_forest.png", p2, width = 8, height = 3.6, dpi = 150)
cat("figures written\n")
