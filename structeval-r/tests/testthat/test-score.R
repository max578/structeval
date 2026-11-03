long_table <- function(type, truth, pred, field = "f") {
  out <- data.frame(
    note = as.character(seq_along(truth)), field = field, type = type,
    stringsAsFactors = FALSE
  )
  out$truth <- truth
  out$pred <- pred
  out
}
counts_of <- function(x) as.matrix(x[, c("cor", "inc", "mis", "spu", "tn")])

test_that("binary counts cover TP, FN, FP, TN", {
  d <- long_table("binary", list("TRUE", "TRUE", "FALSE", "FALSE"),
    list("TRUE", "FALSE", "TRUE", "FALSE"))
  expect_equal(unname(counts_of(score(d))), rbind(
    c(1, 0, 0, 0, 0), c(0, 0, 1, 0, 0), c(0, 0, 0, 1, 0), c(0, 0, 0, 0, 1)
  ))
})

test_that("binary missing predictions follow the policy", {
  d <- long_table("binary", list("TRUE", "FALSE", "TRUE"),
    list(NA_character_, NA_character_, "-"))
  wrong <- score(d, missing = "wrong")
  expect_equal(unname(counts_of(wrong)), rbind(
    c(0, 0, 1, 0, 0), c(0, 0, 0, 1, 0), c(0, 0, 1, 0, 0)
  ))
  expect_warning(skipped <- score(d, missing = "skip"), "No scored notes")
  expect_equal(nrow(skipped), 0L)
})

test_that("scalar counts cover every branch including both empty", {
  d <- long_table("scalar",
    list("A", "A", "A", character(0), character(0)),
    list(" a ", "B", character(0), "B", character(0)))
  expect_equal(unname(counts_of(score(d))), rbind(
    c(1, 0, 0, 0, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0),
    c(0, 0, 0, 1, 0), c(0, 0, 0, 0, 1)
  ))
})

test_that("scalar missing predictions follow the policy", {
  d <- long_table("scalar", list("A", character(0)),
    list(NA_character_, NA_character_))
  expect_equal(unname(counts_of(score(d))), rbind(
    c(0, 0, 1, 0, 0), c(0, 0, 0, 0, 1)
  ))
  expect_warning(skipped <- score(d, missing = "skip"), "No scored notes")
  expect_equal(nrow(skipped), 0L)
})

test_that("list counts use set operations after normalisation", {
  d <- long_table("list",
    list(c("XELOX", "Docetaxel"), character(0), character(0), "A",
      c("a  b", "C")),
    list(c("xelox", "XELOX", "Cisplatin"), character(0), "B", character(0),
      c("A B", "-")))
  expect_equal(unname(counts_of(score(d))), rbind(
    c(1, 0, 1, 1, 0), c(0, 0, 0, 0, 1), c(0, 0, 0, 1, 0),
    c(0, 0, 1, 0, 0), c(1, 0, 1, 0, 0)
  ))
})

test_that("unlabelled rows are skipped, not scored", {
  d <- long_table("list", list(NA_character_, "A", NA),
    list("A", "A", "B"))
  s <- score(d)
  expect_equal(s$note, "2")
})

test_that("a field absent in every note gives NA metrics with a warning", {
  d <- rbind(
    long_table("scalar", list("A", "B"), list("A", "C"), field = "present"),
    long_table("list", list(character(0), character(0)),
      list(character(0), character(0)), field = "absent")
  )
  s <- score(d)
  expect_equal(sum(s$tn[s$field == "absent"]), 2)
  expect_warning(
    m <- metrics(s, field = "absent"),
    "'precision', 'recall', 'f1' undefined for field 'absent'"
  )
  # List fields report no specificity (0.1.1, base paper Table 1).
  expect_true(all(is.na(m$estimate)))
  expect_equal(m$metric, c("precision", "recall", "f1"))
  expect_error(metrics(s, field = "absent", metric = "specificity"),
    "not defined for list")
})

test_that("a field with no labelled row is reported, never dropped", {
  d <- rbind(
    long_table("scalar", list("A", character(0)), list("A", character(0)),
      field = "kept"),
    long_table("scalar", list(NA_character_), list("A"), field = "gone")
  )
  expect_warning(s <- score(d), "No scored notes for field\\(s\\): 'gone'")
  expect_warning(m <- metrics(s), "undefined for field 'gone'")
  expect_true("gone" %in% m$field)
  expect_true(all(is.na(m$estimate[m$field == "gone"])))
})

test_that("bad input is rejected with a named error", {
  d <- long_table("scalar", list(c("A", "B")), list("A"))
  expect_error(score(d), "more than one item")
  d <- long_table("binary", list("maybe"), list("TRUE"))
  expect_error(score(d), "TRUE or FALSE")
  expect_error(score(data.frame(x = 1)), "must be a data frame")
})

test_that("read_oncoshot reads the Oncoshot sample file unchanged", {
  csv <- test_path("..", "fixtures", "crosscheck.csv")
  long <- read_oncoshot(csv, types = c(
    "Has metastasis" = "binary", "Treatment Drugs" = "list"
  ))
  expect_equal(nrow(long), 120L)
  expect_true(is.list(long$truth))
  row2 <- long[long$note == "2" & long$field == "Treatment Drugs", ]
  expect_equal(row2$truth[[1]], c("Docetaxel", "Gemcitabine"))
  row5 <- long[long$note == "5" & long$field == "Has metastasis", ]
  expect_true(is.na(row5$truth[[1]]))
  expect_error(
    read_oncoshot(csv, types = c(Nope = "scalar")),
    "Columns not found"
  )
})

test_that("read_oncoshot accepts comma-separated lists and '-' cells", {
  csv <- tempfile(fileext = ".csv")
  writeLines(c(
    "Patient ID,Drugs,Res: Drugs",
    "1,\"A, B\",\"b, C\"",
    "2,-,",
    "3,,A"
  ), csv)
  long <- read_oncoshot(csv, types = c(Drugs = "list"))
  expect_equal(long$truth[[1]], c("A", "B"))
  expect_equal(long$truth[[2]], character(0))
  expect_true(is.na(long$pred[[2]]))
  s <- score(long)
  expect_equal(s$cor, c(1L, 0L))
  expect_equal(s$tn, c(0L, 1L))
})
