# Calibration regression gate. Seeded and short; each check is within three
# Monte Carlo standard errors of the nominal rate.

within_3se <- function(rate, target, reps) {
  abs(rate - target) <= 3 * sqrt(target * (1 - target) / reps)
}

test_that("Wilson covers 0.95 at n = 50, p = 0.90", {
  reps <- 500L
  set.seed(101)
  k <- rbinom(reps, 50, 0.90)
  w <- structeval:::.wilson_ci(k, 50, 0.95)
  cover <- mean(w[, "lower"] <= 0.90 & 0.90 <= w[, "upper"])
  expect_true(within_3se(cover, 0.95, reps))
})

list_f1_coverage <- function(method, reps, seed) {
  set.seed(seed)
  hits <- vapply(seq_len(reps), function(i) {
    counts <- simulate_counts(100, 0.90, m = 2, rho = 0.2, type = "list")
    m <- metrics(counts, metric = "f1", method = method, B = 1000L)
    m$lower <= 0.90 && 0.90 <= m$upper
  }, logical(1L))
  mean(hits)
}

test_that("routed list F1 is cluster_bca and covers 0.95; naive Wilson does not", {
  counts <- simulate_counts(100, 0.90, m = 2, rho = 0.2, type = "list",
    seed = 1)
  expect_equal(metrics(counts, metric = "f1", seed = 1)$method, "cluster_bca")
  reps <- 300L
  auto <- list_f1_coverage("auto", reps, seed = 202)
  expect_true(within_3se(auto, 0.95, reps))
  naive <- list_f1_coverage("naive_wilson", reps, seed = 202)
  expect_lt(naive, 0.92)
})

test_that("cluster_delta paired test holds its size at delta = 0", {
  reps <- 200L
  set.seed(303)
  rej <- vapply(seq_len(reps), function(i) {
    pair <- simulate_counts(150, 0.90, m = 2, rho = 0.2, type = "list",
      delta_logit = 0, kappa = 0.5)
    test_paired(pair$a, pair$b, "field", "f1")$reject
  }, logical(1L))
  expect_true(within_3se(mean(rej), 0.05, reps))
})
