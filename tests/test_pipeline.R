library(testthat)
suppressPackageStartupMessages(library(dplyr))
source(file.path("..", "R", "adam.R"))
source(file.path("..", "R", "analysis.R"))

test_that("ADSL has one row per randomized patient and the pilot's populations", {
  adsl <- derive_adsl()
  expect_equal(nrow(adsl), 254)
  expect_false(anyDuplicated(adsl$USUBJID) > 0)
  expect_equal(sum(adsl$EFFFL == "Y"), 234)
  expect_true(all(adsl$AGEGR1[adsl$AGE < 65] == "<65"))
  expect_true(all(adsl$AGEGR1[adsl$AGE == 80] == "65-80"))
})

test_that("ADADAS keeps one analysis record per patient, parameter and visit", {
  adsl <- derive_adsl()
  a <- derive_adadas(adsl) %>% filter(ANL01FL == "Y")
  expect_equal(nrow(distinct(a, USUBJID, PARAMCD, AVISIT)), nrow(a))
  expect_true(all(a$CHG == a$AVAL - a$BASE, na.rm = TRUE))
  expect_true(all(a$AVISIT[a$ADY >= 141] == "Week 24"))
})

test_that("ICE is the first visit after the last observed one, intermittent gaps are not ICEs", {
  d <- data.frame(USUBJID = rep(c("a", "b", "c"), each = 3), AVISIT = factor(rep(VIS, 3), VIS),
                  CHG = c(1, NA, NA, 1, NA, 2, 1, 2, 3))
  ice <- ice_data(d, "JR")
  expect_equal(as.character(ice$USUBJID), "a")
  expect_equal(as.character(ice$AVISIT), "Week 16")
})

test_that("Week 24 ANCOVA returns one effect per active arm", {
  set.seed(1)
  d <- data.frame(AVISIT = "Week 24", TRTP = factor(rep(ARMS, each = 20), ARMS), BASE = rnorm(60, 20, 5),
                  SITEGR1 = factor(rep(c("1", "2"), 30)))
  d$CHG <- 0.1 * d$BASE + rnorm(60)
  r <- ancova24(d)
  expect_named(r, c("trt_2", "trt_3"))
  expect_true(all(sapply(r, function(x) is.finite(x$est) && x$se > 0)))
})
