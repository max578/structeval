test_that("min_certifiable(50, 0.90) is 50 because 49/50 does not certify", {
  w49 <- structeval:::.wilson_ci(49, 50, 0.95)
  w50 <- structeval:::.wilson_ci(50, 50, 0.95)
  expect_lt(w49[, "lower"], 0.90)
  expect_equal(round(unname(w49[, "lower"]), 3), 0.895)
  expect_gt(w50[, "lower"], 0.90)
  expect_identical(min_certifiable(50, 0.90), 50L)
  expect_identical(min_certifiable(10, 0.99), NA_integer_)
})

test_that("n_threshold gives 239 units and DEFF 1.6 at m = 4, rho = 0.2", {
  nt <- n_threshold(0.90, 0.95)
  expect_identical(nt$n_units, 239L)
  expect_identical(nt$n_notes, 239L)
  nt4 <- n_threshold(0.90, 0.95, m = 4, rho = 0.2, prevalence = 0.5)
  expect_equal(nt4$deff, 1.6)
  expect_identical(nt4$n_items, as.integer(ceiling(239 * 1.6)))
  expect_identical(nt4$n_present, as.integer(ceiling(ceiling(239 * 1.6) / 4)))
  expect_identical(nt4$n_notes, 2L * nt4$n_present)
  expect_error(n_threshold(0.95, 0.90), "greater than")
})

test_that("n_paired follows the Connor formula", {
  np <- n_paired(0.05, 0.10)
  z <- qnorm(0.975) * sqrt(0.10) + qnorm(0.80) * sqrt(0.10 - 0.05^2)
  expect_identical(np$n_present, as.integer(ceiling(z^2 / 0.05^2)))
  expect_error(n_paired(0.4, 0.1), "exceed")
})

test_that("power functions return rates with Monte Carlo standard errors", {
  p1 <- power_threshold(60, 0.80, 0.95, nsim = 50, seed = 1)
  expect_equal(p1$power_se, sqrt(p1$power * (1 - p1$power) / 50))
  p3 <- power_asymmetry(60, 0.9, m = 3, nsim = 30, n_perm = 200, seed = 1)
  expect_true(p3$reject >= 0 && p3$reject <= 1)
  p2 <- power_paired(60, delta_logit = 0, nsim = 30, seed = 1)
  expect_true(p2$reject >= 0 && p2$reject <= 1)
})
