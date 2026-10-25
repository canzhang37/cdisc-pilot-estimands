# derive ADaM, write XPT + RDS
source("R/adam.R")

adsl <- derive_adsl()
adadas <- derive_adadas(adsl)
adae <- derive_adae(adsl)

dir.create("adam", showWarnings = FALSE)
for (nm in c("adsl", "adadas", "adae")) {
  d <- get(nm)
  saveRDS(d, file.path("adam", paste0(nm, ".rds")))
  haven::write_xpt(d, file.path("adam", paste0(nm, ".xpt")), version = 5, name = toupper(nm))
}
cat("ADSL", nrow(adsl), "| ADADAS", nrow(adadas), "| ADAE", nrow(adae), "\n")
print(table(adsl$TRT01P, adsl$EFFFL))
