test_that("F1 as a ratio equals 2PR/(P+R) on random counts", {
  set.seed(11)
  for (i in seq_len(200)) {
    cor <- rpois(1, 30) + 1
    inc <- rpois(1, 3)
    mis <- rpois(1, 4)
    spu <- rpois(1, 4)
    p <- cor / (cor + spu + inc)
    r <- cor / (cor + mis + inc)
    f1 <- 2 * cor / (2 * cor + mis + spu + 2 * inc)
    expect_equal(f1, 2 * p * r / (p + r), tolerance = 1e-12)
  }
})

test_that("the router sends note-level proportions to Wilson and ratios to BCa", {
  sc <- simulate_counts(80, 0.9, type = "proportion", seed = 1)
  m <- metrics(sc, metric = c("recall", "f1"), seed = 1)
  expect_equal(m$method, c("wilson", "cluster_bca"))
  li <- simulate_counts(80, 0.9, m = 3, rho = 0.2, type = "list", seed = 1)
  m <- metrics(li, metric = c("precision", "recall", "f1"), seed = 1)
  expect_true(all(m$method == "cluster_bca"))
})

test_that("an all-correct list field takes the boundary rule, never zero width", {
  counts <- data.frame(
    note = as.character(1:30), field = "drugs", type = "list",
    cor = rep(1:3, 10), inc = 0L, mis = 0L, spu = 0L, tn = 0L
  )
  for (meth in c("auto", "cluster_bca", "cluster_delta")) {
    m <- metrics(counts, metric = "f1", method = meth, seed = 1)
    expect_true(m$boundary)
    expect_equal(m$estimate, 1)
    expect_gt(m$upper - m$lower, 0)
    # Wilson on notes, the cap for every sample on these notes (0.1.1).
    w <- structeval:::.wilson_ci(30, 30, 0.95)
    expect_equal(m$lower, unname(w[, 1]))
    expect_equal(m$upper, 1)
  }
})

test_that("bootstrap intervals are reproducible and leave the RNG alone", {
  li <- simulate_counts(60, 0.85, m = 3, rho = 0.3, type = "list", seed = 5)
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  a <- metrics(li, metric = "f1", seed = 7)
  after <- runif(1)
  b <- metrics(li, metric = "f1", seed = 7)
  expect_identical(a, b)
  expect_identical(before, after)
})

test_that("cluster bootstrap meets the speed target at n = 400, B = 2000", {
  set.seed(3)
  a <- rpois(400, 3)
  b <- a + rpois(400, 1)
  elapsed <- system.time(
    structeval:::ratio_ci(a, b, "cluster_bca", B = 2000L)
  )[["elapsed"]]
  expect_lt(elapsed, 0.5)
})

test_that("Table 1 point estimates reproduce to the printed 0.1 percent", {
  counts <- table1_counts()
  m <- metrics(counts, method = "naive_wilson")
  key <- paste(m$field, m$metric)
  got <- round(100 * m$estimate[match(
    paste(table1_published$field, table1_published$metric), key
  )], 1)
  expect_false(anyNA(got))
  expect_equal(max(abs(got - table1_published$value)), 0)
})
