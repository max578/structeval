# Certification is monotone in errors at the boundary: a perfect record keeps
# Wilson on notes, and no sample on the same notes with fewer correct items
# gets a lower bound above it (symmetrically for the upper bound at 0).

list_counts <- function(cor, mis, inc = 0L * cor, spu = 0L * cor) {
  data.frame(
    note = as.character(seq_along(cor)), field = "drugs", type = "list",
    cor = cor, inc = inc, mis = mis, spu = spu, tn = 0L
  )
}

test_that("200/200 certifies whenever 199/200 does (50 notes x 4 items)", {
  perfect <- list_counts(rep(4L, 50), rep(0L, 50))
  one <- list_counts(c(3L, rep(4L, 49)), c(1L, rep(0L, 49)))
  for (meth in c("auto", "cluster_bca", "cluster_delta")) {
    tp <- test_threshold(perfect, "drugs", "recall", theta0 = 0.95,
      method = meth, seed = 1)
    to <- test_threshold(one, "drugs", "recall", theta0 = 0.95,
      method = meth, seed = 1)
    expect_gte(tp$lower, to$lower)
    expect_gte(as.numeric(tp$reject), as.numeric(to$reject))
    expect_lte(tp$p_value, to$p_value)
    expect_true(tp$boundary)
    expect_lt(tp$lower, 1)
  }
})

test_that("the boundary interval is never below any one-step neighbour", {
  set.seed(42)
  designs <- list(
    rep(1L, 20), rep(2L, 10), rep(4L, 50), rep(4L, 5),
    1L + rpois(30, 3), 1L + rpois(100, 1), c(1L, 8L), rep(3L, 200)
  )
  for (b in designs) {
    for (meth in c("cluster_bca", "cluster_delta")) {
      for (seed in 1:3) {
        top <- structeval:::ratio_ci(b, b, meth, B = 1000L, seed = seed)
        bottom <- structeval:::ratio_ci(0 * b, b, meth, B = 1000L,
          seed = seed)
        expect_true(top$boundary && bottom$boundary)
        expect_gt(top$upper - top$lower, 0)
        expect_gt(bottom$upper - bottom$lower, 0)
        # The boundary keeps Wilson on notes; the neighbours are capped.
        w1 <- structeval:::.wilson_ci(length(b), length(b), 0.95)
        w0 <- structeval:::.wilson_ci(0, length(b), 0.95)
        expect_equal(top$lower, unname(w1[, 1]))
        expect_equal(bottom$upper, unname(w0[, 2]))
        for (j in seq_along(b)) {
          a1 <- b
          a1[j] <- a1[j] - 1
          nb <- structeval:::ratio_ci(a1, b, meth, B = 1000L, seed = seed)
          expect_gte(top$lower, nb$lower)
          b1 <- b
          b1[j] <- b1[j] + 1
          nb <- structeval:::ratio_ci(b, b1, meth, B = 1000L, seed = seed)
          expect_gte(top$lower, nb$lower)
          a0 <- 0 * b
          a0[j] <- 1
          nb0 <- structeval:::ratio_ci(a0, b, meth, B = 1000L, seed = seed)
          expect_lte(bottom$upper, nb0$upper)
          nb0 <- structeval:::ratio_ci(a0, b1, meth, B = 1000L, seed = seed)
          expect_lte(bottom$upper, nb0$upper)
        }
      }
    }
  }
})

test_that("certification decisions are monotone over a threshold sweep", {
  set.seed(7)
  for (b in list(rep(4L, 50), 1L + rpois(40, 2), rep(2L, 25))) {
    for (meth in c("cluster_bca", "cluster_delta")) {
      for (theta0 in c(0.80, 0.90, 0.93, 0.95, 0.97, 0.98)) {
        tp <- structeval:::.t1_core(b, b, theta0, 0.025, meth, 1000L)
        for (j in seq_along(b)) {
          a1 <- b
          a1[j] <- a1[j] - 1
          to <- structeval:::.t1_core(a1, b, theta0, 0.025, meth, 1000L)
          if (isTRUE(to$reject) && meth == "cluster_delta") {
            expect_true(tp$reject)
          }
          if (meth == "cluster_delta") {
            expect_lte(tp$p_value, to$p_value)
          }
        }
        # The delta p-value agrees with the decision at the boundary.
        if (meth == "cluster_delta") {
          expect_identical(tp$reject, tp$p_value < 0.025)
        }
      }
    }
  }
})

test_that("one missed, spurious or incorrect item never certifies more", {
  set.seed(3)
  cor <- 1L + rpois(40, 2)
  zero <- 0L * cor
  perfect <- list_counts(cor, zero)
  for (j in c(1L, which.max(cor), which.min(cor))) {
    hit <- replace(zero, j, 1L)
    worse <- list(
      missed = list_counts(cor - hit, hit),
      spurious = list_counts(cor, zero, spu = hit),
      incorrect = list_counts(cor - hit, zero, inc = hit)
    )
    for (metric in c("precision", "recall", "f1")) {
      for (meth in c("cluster_bca", "cluster_delta")) {
        tp <- test_threshold(perfect, "drugs", metric, theta0 = 0.9,
          method = meth, seed = 2)
        for (w in worse) {
          tw <- test_threshold(w, "drugs", metric, theta0 = 0.9,
            method = meth, seed = 2)
          # A missed item is no error for precision: still a perfect record.
          if (tw$estimate == 1) next
          expect_gte(tp$lower, tw$lower)
          expect_gte(as.numeric(tp$reject), as.numeric(tw$reject))
        }
      }
    }
  }
})

test_that("no sample with errors gets a lower bound above the perfect record", {
  set.seed(11)
  for (b in list(rep(4L, 50), 1L + rpois(40, 2), rep(2L, 25), rep(1L, 30))) {
    nn <- length(b)
    l1 <- unname(structeval:::.wilson_ci(nn, nn, 0.95)[, 1])
    u0 <- unname(structeval:::.wilson_ci(0, nn, 0.95)[, 2])
    for (meth in c("cluster_bca", "cluster_delta")) {
      for (e in c(1L, 2L, 3L, 5L)) {
        for (seed in 1:3) {
          j <- sample.int(nn, e, replace = TRUE)
          a1 <- b
          for (jj in j) a1[jj] <- max(a1[jj] - 1L, 0L)
          hi <- structeval:::ratio_ci(a1, b, meth, B = 1000L, seed = seed)
          expect_lte(hi$lower, l1)
          a0 <- 0L * b
          for (jj in j) a0[jj] <- min(a0[jj] + 1L, b[jj])
          lo <- structeval:::ratio_ci(a0, b, meth, B = 1000L, seed = seed)
          expect_gte(lo$upper, u0)
          # The T1 p-value agrees with the decision off the boundary too.
          if (meth == "cluster_delta") {
            for (theta0 in c(0.85, 0.9, 0.95)) {
              t1 <- structeval:::.t1_core(a1, b, theta0, 0.025, meth, 1000L)
              expect_identical(t1$reject, t1$p_value < 0.025)
            }
          }
        }
      }
    }
  }
})
