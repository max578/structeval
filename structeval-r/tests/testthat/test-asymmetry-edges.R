test_that("cluster_t on a one-note field returns NA with a warning", {
  one <- data.frame(
    note = "1", field = "drugs", type = "list",
    cor = 2L, inc = 0L, mis = 3L, spu = 1L, tn = 0L
  )
  expect_warning(
    res <- test_asymmetry(one, "drugs", method = "cluster_t"),
    "fewer than two notes"
  )
  expect_true(is.na(res$p_value))
  expect_true(is.na(res$reject))
  expect_equal(res$n, 4)
  # The other methods are still defined on one note.
  expect_false(is.na(test_asymmetry(one, "drugs", method = "binomial")$p_value))
})
