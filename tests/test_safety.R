library(testthat)
ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
source(file.path("..", "R", "safety.R"))

test_that("V1: Newcombe hybrid score interval matches Newcombe (1998), 56/70 vs 48/80", {
  r <- newcombe(56, 70, 48, 80)
  expect_equal(unname(r["est"]), 0.2)
  expect_equal(unname(round(r[c("lower", "upper")], 4)), c(0.0524, 0.3339))
})

test_that("Newcombe interval stays inside [-1, 1] at the boundaries", {
  r <- newcombe(0, 72, 0, 86)
  expect_equal(unname(r["est"]), 0)
  expect_true(r["lower"] < 0 && r["upper"] > 0)
  expect_true(all(abs(newcombe(72, 72, 0, 86)) <= 1))
})

test_that("V2: double FDR with one SOC equals Benjamini-Hochberg at 0.1", {
  set.seed(3)
  p <- c(runif(40), 1e-4, 0.002, 0.03)
  expect_equal(double_fdr(p, rep("one", length(p))), p.adjust(p, "BH") <= 0.1)
})

test_that("double FDR drops every PT of a SOC that fails the first step", {
  p <- c(0.001, 0.5, 0.6, 0.9)
  soc <- c("a", "a", "b", "b")
  expect_equal(double_fdr(p, soc), c(TRUE, FALSE, FALSE, FALSE))
})
