# Exact power of the package's McNemar test (exact conditional binomial on
# the discordant notes, two-sided at 2 * alpha), enumerated with binom.test.
mcnemar_power_brute <- function(n, delta, psi, alpha) {
  p10 <- (psi + delta) / 2
  total <- 0
  for (d in seq_len(n)) {
    pv <- vapply(0:d, function(x) binom.test(x, d)$p.value, numeric(1L))
    total <- total + dbinom(d, n, psi) *
      sum(dbinom(0:d, d, p10 / psi)[pv < 2 * alpha])
  }
  total
}

test_that("n_paired delivers the target power of the exact McNemar test", {
  for (cfg in list(c(0.20, 0.30), c(0.15, 0.25), c(0.10, 0.20))) {
    np <- n_paired(cfg[1], cfg[2])
    z <- qnorm(0.975) * sqrt(cfg[2]) + qnorm(0.80) * sqrt(cfg[2] - cfg[1]^2)
    n_connor <- ceiling(z^2 / cfg[1]^2)
    expect_identical(np$n_connor, as.integer(n_connor))
    expect_gte(np$n_present, n_connor)
    exact <- mcnemar_power_brute(np$n_present, cfg[1], cfg[2], 0.025)
    expect_equal(np$exact_power, exact, tolerance = 1e-10)
    expect_gte(exact, 0.80)
    # The first n from Connor's n upward that reaches the target.
    if (np$n_present > n_connor) {
      below <- vapply(n_connor:(np$n_present - 1L), function(n) {
        mcnemar_power_brute(n, cfg[1], cfg[2], 0.025)
      }, numeric(1L))
      expect_true(all(below < 0.80))
    }
  }
  expect_identical(n_paired(0.20, 0.30, prevalence = 0.5)$n_notes,
    2L * n_paired(0.20, 0.30)$n_present)
})
