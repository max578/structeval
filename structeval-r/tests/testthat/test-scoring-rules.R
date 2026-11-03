# Scoring rules hardened in 0.1.1 (adversarial review SC-01 to SC-18, K1, K2).
# Each test names the finding it closes.

lt <- function(type, truth, pred, field = "f") {
  out <- data.frame(
    note = as.character(seq_along(truth)), field = field, type = type,
    stringsAsFactors = FALSE
  )
  out$truth <- truth
  out$pred <- pred
  out
}
cnt <- function(d, ...) {
  unname(as.matrix(score(d, ...)[, c("cor", "inc", "mis", "spu", "tn")]))
}
one <- function(cor = 0, inc = 0, mis = 0, spu = 0, tn = 0) {
  c(cor, inc, mis, spu, tn)
}
csv_file <- function(lines) {
  f <- tempfile(fileext = ".csv")
  writeLines(lines, f, useBytes = TRUE)
  f
}

# SC-01, SC-15: the reader -----------------------------------------------------

test_that("SC-01 a row with an extra field stops the reader, naming the line", {
  f <- csv_file(c("Patient ID,Drugs,Res: Drugs",
    "1,['A'],['B']", "2,['B'],['B'],extra", "3,['C'],['D']"))
  expect_error(read_oncoshot(f, c(Drugs = "list")), "line 3")
  late <- c("Patient ID,Drugs,Res: Drugs", paste0(1:6, ",['A'],['A']"),
    "7,['B'],['C'],['B']")
  expect_error(read_oncoshot(csv_file(late), c(Drugs = "list")), "line 8")
})

test_that("SC-15 short rows and unbalanced quotes stop the reader", {
  f <- csv_file(c("Patient ID,Drugs,Res: Drugs", "1,['A'],['A']", "2,['B']"))
  expect_error(read_oncoshot(f, c(Drugs = "list")), "fields")
  f <- csv_file(c("Patient ID,Drugs,Res: Drugs", "1,\"['A'],['A']",
    "2,['B'],['B']", "3,['C'],['C']"))
  expect_error(read_oncoshot(f, c(Drugs = "list")), "fields")
})

test_that("SC-15 a quoted cell spanning two lines is still read", {
  f <- csv_file(c("Patient ID,Hist,Res: Hist", "1,\"Adeno", "carcinoma\",x",
    "2,A,A"))
  long <- read_oncoshot(f, c(Hist = "scalar"))
  expect_equal(long$note, c("1", "2"))
})

test_that("SC-15 duplicated needed columns stop the reader", {
  f <- csv_file(c("Patient ID,Drugs,Res: Drugs,Res: Drugs",
    "1,['A'],['X'],['A']"))
  expect_error(read_oncoshot(f, c(Drugs = "list")), "more than once")
})

test_that("SC-15 invalid UTF-8 stops the reader", {
  f <- tempfile(fileext = ".csv")
  writeBin(c(charToRaw("Patient ID,Hist,Res: Hist\n1,Sj"), as.raw(0xf6),
    charToRaw("gren,x\n2,A,A\n")), f)
  expect_error(read_oncoshot(f, c(Hist = "scalar")), "UTF-8")
})

test_that("SC-15 an empty note identifier stops the reader", {
  f <- csv_file(c("Patient ID,Hist,Res: Hist", "1,A,A", " ,B,B"))
  expect_error(read_oncoshot(f, c(Hist = "scalar")), "empty")
})

# SC-17 ------------------------------------------------------------------------

test_that("SC-17 duplicated note and field rows stop score()", {
  d <- lt("scalar", list("A", "B"), list("A", "C"))
  d <- rbind(d, d[1, ])
  expect_error(score(d), "more than once")
})

# SC-02: binary vocabulary ------------------------------------------------------

test_that("SC-02 an unrecognised binary prediction stops and lists the values", {
  d <- lt("binary", list("TRUE", "TRUE", "FALSE", "FALSE"),
    list("Present", "Yes.", "Absent", "negative"))
  expect_error(score(d), "Present.*Yes\\..*Absent.*negative")
})

test_that("SC-02 recognised binary spellings, including 1.0 and 0.0, are flags", {
  d <- lt("binary", list("TRUE", "FALSE", "TRUE", "FALSE", "yes"),
    list("Yes", "no", "1.0", "0.0", " TRUE\u00a0"))
  expect_equal(cnt(d), rbind(one(cor = 1), one(tn = 1), one(cor = 1),
    one(tn = 1), one(cor = 1)))
})

test_that("SC-02 an empty binary prediction follows the missing policy", {
  d <- lt("binary", list("TRUE", "FALSE", "TRUE", "FALSE"),
    list("-", "None", "", "\u2013"))
  expect_equal(cnt(d), rbind(one(mis = 1), one(spu = 1), one(mis = 1),
    one(spu = 1)))
  expect_message(
    expect_warning(s <- score(d, missing = "skip"), "No scored notes"),
    "missing prediction"
  )
  expect_equal(nrow(s), 0L)
})

# SC-03, K2: numbers written as text ---------------------------------------------

test_that("K2 plain decimals compare as numbers", {
  d <- lt("scalar", list("4", "12", "4.50", "-0", "0.5", "100"),
    list("4.0", "12.00", "4.5", "0", "0.50", "100."))
  expect_equal(cnt(d)[, 1], c(1, 1, 1, 1, 1, 0))
})

test_that("K2 codes, exponents and special values stay text", {
  d <- lt("scalar", list("007", "02", "1e3", "inf", "1_000", ".5", "+4"),
    list("7", "2", "1000", "Infinity", "1000", "0.5", "4"))
  expect_equal(cnt(d)[, 2], rep(1, 7))
  d <- lt("scalar", list("007", "1E3"), list("007", "1e3"))
  expect_equal(cnt(d)[, 1], c(1, 1))
})

test_that("K2 a NaN label is text and is never dropped", {
  d <- lt("scalar", list("NaN", "NaN"), list("nan", "4"))
  expect_equal(cnt(d), rbind(one(cor = 1), one(inc = 1)))
})

test_that("K2 numeric list items compare as numbers", {
  d <- lt("list", list(c("4", "5.0", "007")), list(c("4.0", "5", "7")))
  expect_equal(cnt(d), rbind(one(cor = 2, mis = 1, spu = 1)))
})

# SC-04: null-like tokens ---------------------------------------------------------

test_that("SC-04 'null' and 'None' are empty values on both sides", {
  d <- lt("scalar", list("Adeno", "Adeno", "None", "null", character(0)),
    list("null", "NONE", "-", "Adeno", "None"))
  expect_equal(cnt(d), rbind(one(mis = 1), one(mis = 1), one(tn = 1),
    one(spu = 1), one(tn = 1)))
})

test_that("SC-04 N/A and NA (any case) are empty; unknown stays a value", {
  d <- lt("scalar", list("Adeno", "Adeno", "Adeno", "N/A", "na", "Adeno"),
    list("N/A", "NA", "unknown", "n/a", "B", "nA"))
  expect_equal(cnt(d), rbind(one(mis = 1), one(mis = 1), one(inc = 1),
    one(tn = 1), one(spu = 1), one(mis = 1)))
  b <- lt("binary", list("TRUE", "NA"), list("N/A", "TRUE"))
  expect_message(s <- score(b), "no information")
  expect_equal(cnt(b), rbind(one(mis = 1)))
})

test_that("SC-04 the literal NA is the same through read_oncoshot()", {
  f <- csv_file(c("Patient ID,Hist,Res: Hist", "1,Adeno,NA", "2,NA,n/a",
    "3,N/A,B"))
  via_reader <- read_oncoshot(f, c(Hist = "scalar"))
  direct <- lt("scalar", list("Adeno", "NA", "N/A"), list("NA", "n/a", "B"),
    field = "Hist")
  expect_equal(cnt(via_reader), cnt(direct))
  expect_equal(cnt(direct), rbind(one(mis = 1), one(tn = 1), one(spu = 1)))
})

test_that("SC-04 en and em dashes mean no information, like '-'", {
  d <- lt("scalar", list("\u2013", "-", "\u2014", "A"),
    list(character(0), "\u2013", "A", "\u2014"))
  expect_equal(cnt(d), rbind(one(tn = 1), one(tn = 1), one(spu = 1),
    one(mis = 1)))
})

test_that("SC-04 null_tokens = character(0) keeps null-like text as values", {
  d <- lt("scalar", list("Adeno"), list("None"))
  expect_equal(cnt(d, null_tokens = character(0)), rbind(one(inc = 1)))
})

# SC-05: bracket lists -------------------------------------------------------------

test_that("SC-05 malformed bracket lists stop with the cell named", {
  p <- structeval:::.parse_cell
  expect_error(p("['A', B]", "list"), "\\['A', B\\]")
  expect_error(p("['A', 'B'].", "list"), "bracket")
  expect_error(p("['A'; 'B']", "list"), "bracket")
})

test_that("SC-05 escaped quotes are unescaped; bare and blank items are handled", {
  p <- structeval:::.parse_cell
  expect_equal(p("['O\\'Brien']", "list"), "O'Brien")
  expect_equal(p("[\"O'Brien\", 'a\\\\b']", "list"), c("O'Brien", "a\\b"))
  expect_equal(p("[A, b]", "list"), c("A", "b"))
  expect_equal(p("['', 'A']", "list"), "A")
  expect_equal(p("['']", "list"), character(0))
  d <- lt("list", list(p("['O\\'Brien']", "list")),
    list(p("[\"O'Brien\"]", "list")))
  expect_equal(cnt(d), rbind(one(cor = 1)))
})

# SC-07, SC-16: blank and no-information references --------------------------------

test_that("SC-07 blank references are unlabelled for every type", {
  for (type in c("binary", "scalar", "list")) {
    d <- lt(type, list("", "   ", c(NA, NA), "\u00a0"),
      list("TRUE", "TRUE", "TRUE", "TRUE"))
    expect_warning(s <- score(d), "No scored notes", info = type)
    expect_equal(nrow(s), 0L, info = type)
    expect_equal(attr(s, "coverage")$unlabelled, 4L, info = type)
  }
})

test_that("SC-16 a '-' binary reference is skipped, counted and reported", {
  d <- lt("binary", list("TRUE", "-", "FALSE"), list("TRUE", "TRUE", "FALSE"))
  expect_message(s <- score(d), "no information")
  expect_equal(s$note, c("1", "3"))
  expect_equal(attr(s, "coverage")$no_information, 1L)
})

# SC-08: case folding ----------------------------------------------------------------

test_that("SC-08 non-ASCII items stop scoring outside a UTF-8 locale", {
  fold <- structeval:::.fold_case
  expect_error(fold("\u00c9tape II", utf8 = FALSE), "UTF-8")
  expect_equal(fold("Stage II", utf8 = FALSE), "stage ii")
  expect_equal(fold("\u00c9tape II", utf8 = TRUE), "\u00e9tape ii")
})

test_that("SC-08 case folding is the simple one-to-one lowercase mapping", {
  d <- lt("scalar", list("\u00c9tape II", "Stra\u00dfe", "\u039f\u0394\u039f\u03a3"),
    list("\u00e9tape ii", "STRASSE", "\u03bf\u03b4\u03bf\u03c3"))
  expect_equal(cnt(d)[, 1], c(1, 0, 1))
})

# SC-09, K1: no specificity for list fields ---------------------------------------------

test_that("K1 list fields report no specificity", {
  d <- lt("list", list(character(0), c("A", "B"), character(0)),
    list(character(0), c("A", "B"), "C"))
  m <- metrics(score(d), seed = 1, B = 50L)
  expect_equal(m$metric, c("precision", "recall", "f1"))
  expect_error(metrics(score(d), metric = "specificity"), "not defined for list")
  expect_error(
    test_threshold(score(d), "f", "specificity", theta0 = 0.9),
    "not defined for list"
  )
  expect_error(
    power_threshold(100, 0.9, 0.8, m = 2, type = "list",
      metric = "specificity", nsim = 10L, seed = 1),
    "not defined for list"
  )
})

# SC-11: whitespace ----------------------------------------------------------------------

test_that("SC-11 every Unicode space is trimmed and internal runs collapse", {
  d <- lt("scalar",
    list("Stage  IIA", "\u00a0Adeno", "\vAdeno", "Adeno\u2003carcinoma",
      "a\t b"),
    list("stage iia", "Adeno", "Adeno", "Adeno carcinoma", "A B"))
  expect_equal(cnt(d)[, 1], rep(1, 5))
})

# SC-12: the coverage record ---------------------------------------------------------------

test_that("SC-12 score() records skipped and missing rows per field", {
  d <- lt("scalar", list("A", "B", "C", "D", NA_character_),
    list("A", "B", NA_character_, NA_character_, "E"))
  s <- score(d)
  cov <- attr(s, "coverage")
  expect_equal(cov$rows, 5L)
  expect_equal(cov$unlabelled, 1L)
  expect_equal(cov$pred_missing, 2L)
  expect_equal(cov$scored, 4L)
  expect_message(k <- score(d, missing = "skip"), "2 row")
  expect_equal(attr(k, "coverage")$scored, 2L)
  m <- metrics(s, metric = "recall")
  expect_equal(attr(m, "coverage"), cov)
})

# SC-13, SC-14: documented conventions, held by tests -------------------------------------

test_that("SC-13 the missing policy differs by type, as documented", {
  sc <- lt("scalar", list(character(0)), list(NA_character_))
  expect_equal(cnt(sc), rbind(one(tn = 1)))
  bi <- lt("binary", list("FALSE"), list(NA_character_))
  expect_equal(cnt(bi), rbind(one(spu = 1)))
})

test_that("SC-14 F1 is 0 when the model returns nothing for a labelled field", {
  d <- lt("list", list("A", c("B", "C"), "D"),
    list(character(0), character(0), character(0)))
  expect_warning(m <- metrics(score(d), seed = 1, B = 50L), "'precision'")
  expect_true(is.na(m$estimate[m$metric == "precision"]))
  expect_equal(m$estimate[m$metric == "recall"], 0)
  expect_equal(m$estimate[m$metric == "f1"], 0)
})
