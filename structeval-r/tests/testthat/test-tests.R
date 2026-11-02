test_that("T1 Wilson decision matches the score-test p-value", {
  sc <- simulate_counts(150, 0.95, type = "proportion", seed = 2)
  t1 <- test_threshold(sc, "field", "recall", theta0 = 0.90)
  expect_equal(t1$method, "wilson")
  expect_equal(t1$reject, t1$p_value < 0.025)
  w <- structeval:::.wilson_ci(t1$n * t1$estimate, t1$n, 0.95)
  expect_equal(t1$lower, unname(w[, 1]))
})

test_that("T1 labels the p-value source for cluster_bca and applies Holm", {
  li <- simulate_counts(120, 0.92, m = 3, rho = 0.2, type = "list", seed = 3)
  li2 <- li
  li2$field <- "other"
  both <- rbind(li, li2)
  t1 <- test_threshold(both, c("field", "other"), "f1", theta0 = 0.85,
    seed = 1, adjust = "holm")
  expect_match(t1$method, "cluster_bca \\(p: cluster_delta\\)")
  expect_equal(t1$p_adjusted, p.adjust(t1$p_value, "holm"))
})

test_that("T2 mcnemar equals the exact binomial test on discordant notes", {
  pair <- simulate_counts(200, 0.85, type = "proportion", delta_logit = 0.6,
    kappa = 0.3, seed = 8)
  t2 <- test_paired(pair$a, pair$b, "field", "recall", method = "mcnemar")
  n10 <- sum(pair$a$cor == 1 & pair$b$cor == 0)
  n01 <- sum(pair$a$cor == 0 & pair$b$cor == 1)
  expect_equal(t2$p_value, binom.test(n10, n10 + n01)$p.value)
  tab <- table(factor(pair$a$cor, 0:1), factor(pair$b$cor, 0:1))
  expect_equal(c(n10, n01), c(tab["1", "0"], tab["0", "1"]))
  expect_error(
    test_paired(pair$a, pair$b, "field", "f1", method = "mcnemar"),
    "same note-level 0/1 denominator"
  )
})

test_that("T2 sign-flip statistic is the difference itself", {
  pair <- simulate_counts(150, 0.9, m = 2, rho = 0.2, type = "list",
    delta_logit = 0.8, kappa = 0.5, seed = 9)
  sf <- test_paired(pair$a, pair$b, "field", "f1", method = "signflip",
    seed = 1)
  cd <- test_paired(pair$a, pair$b, "field", "f1")
  expect_equal(sf$estimate, cd$estimate)
  expect_lt(sf$p_value, 0.05)
  expect_lt(cd$p_value, 0.05)
})

test_that("T3 sign-flip p-value is invariant to added substitutions", {
  li <- simulate_counts(120, 0.85, m = 4, rho = 0.2, spu_ratio = 0.6,
    type = "list", seed = 12)
  base <- test_asymmetry(li, "field", seed = 21)
  set.seed(5)
  hit <- sample(nrow(li), 40, replace = TRUE)
  k <- tabulate(hit, nrow(li))
  sub <- li
  sub$mis <- sub$mis + k
  sub$spu <- sub$spu + k
  moved <- test_asymmetry(sub, "field", seed = 21)
  expect_identical(moved$p_value, base$p_value)
  bin_base <- test_asymmetry(li, "field", method = "binomial")
  bin_sub <- test_asymmetry(sub, "field", method = "binomial")
  expect_gt(bin_sub$p_value, bin_base$p_value)
  expect_lt(abs(bin_sub$estimate - 0.5), abs(bin_base$estimate - 0.5))
})

test_that("T3 cluster_t equals t.test on the note differences", {
  li <- simulate_counts(80, 0.85, m = 3, rho = 0.2, spu_ratio = 0.7,
    type = "list", seed = 13)
  t3 <- test_asymmetry(li, "field", method = "cluster_t")
  expect_equal(t3$p_value, t.test(li$mis - li$spu)$p.value)
})
