# tables T1-T5 as text
source("R/adam.R")
source("R/tlf.R")
adsl <- readRDS("adam/adsl.rds")
adae <- readRDS("adam/adae.rds")
lsm <- read.csv("results/lsmeans_by_visit.csv")
prim <- read.csv("results/primary_mmrm.csv")
sens <- read.csv("results/sensitivity.csv")

dir.create("tlf", showWarnings = FALSE)
tabs <- list(t1_demographics = t_dm(adsl), t2_disposition = t_ds(adsl), t3_primary = t_primary(lsm, prim),
             t4_sensitivity = t_sens(sens), t5_teae = t_ae(adsl, adae))
for (nm in names(tabs)) {
  formatters::export_as_txt(tabs[[nm]], file = file.path("tlf", paste0(nm, ".txt")), paginate = FALSE)
  cat("\n"); print(tabs[[nm]])
}
