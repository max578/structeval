# Agreement with the Python reference implementation on two shared fixtures:
# crosscheck.csv (every edge case of the counting rules plus random rows) and
# crosscheck_adversarial.csv (inputs the two scored differently before 0.1.1).

fixtures <- list(
  crosscheck = c(
    "Has metastasis" = "binary", "First Primary Diagnosis" = "scalar",
    "First Primary Histology" = "scalar", "Treatment Drugs" = "list",
    "Radiotherapy Sites" = "list"
  ),
  crosscheck_adversarial = c(
    Met = "binary", Hist = "scalar", Size = "scalar", Drugs = "list"
  )
)

for (fx in names(fixtures)) {
  test_that(paste("counts, estimates and Wilson bounds match Python:", fx), {
    types <- fixtures[[fx]]
    fixture <- test_path("..", "fixtures", paste0(fx, ".csv"))
    reference <- test_path("..", "fixtures", paste0(fx, "_python.json"))
    skip_if_not(file.exists(reference), "Python reference output not built")
    skip_if_not_installed("jsonlite")
    py <- jsonlite::fromJSON(reference, simplifyVector = FALSE)$fields
    counts <- suppressMessages(score(read_oncoshot(fixture, types)))
    expect_setequal(names(py), names(types))

    for (f in names(types)) {
      ref <- py[[f]]
      rows <- counts[counts$field == f, ]

      # Per-note counts ---------------------------------------------------------
      py_names <- if (ref$type == "binary") {
        c(cor = "TP", spu = "FP", mis = "FN", tn = "TN")
      } else {
        c(cor = "cor", inc = "inc", mis = "mis", spu = "spu", tn = "tn")
      }
      expect_setequal(rows$note, names(ref$per_note))
      for (j in seq_along(py_names)) {
        py_vals <- vapply(ref$per_note[rows$note], function(x) {
          x[[py_names[[j]]]]
        }, numeric(1L))
        expect_equal(rows[[names(py_names)[j]]], unname(py_vals),
          tolerance = 0, info = paste(f, names(py_names)[j])
        )
      }

      # Point estimates and Wilson bounds ---------------------------------------
      method <- if (ref$type == "list") "naive_wilson" else "auto"
      wanted <- names(ref$estimates)
      # Every metric the R package reports for the type is compared.
      expect_setequal(wanted, structeval:::.default_metrics(ref$type))
      if (all(vapply(ref$estimates, is.null, logical(1L)))) {
        expect_warning(m <- metrics(counts, field = f, metric = wanted),
          "undefined")
      } else {
        m <- metrics(counts, field = f, metric = wanted, method = method,
          seed = 1)
      }
      for (metric in wanted) {
        r <- m[m$metric == metric, ]
        py_est <- ref$estimates[[metric]]
        if (is.null(py_est)) {
          expect_true(is.na(r$estimate), info = paste(f, metric))
        } else {
          expect_equal(r$estimate, py_est, tolerance = 1e-12,
            info = paste(f, metric))
        }
        py_ci <- ref$wilson[[metric]]
        if (!is.null(py_ci) && all(vapply(py_ci, is.null, logical(1L)))) {
          expect_true(all(is.na(c(r$lower, r$upper))), info = paste(f, metric))
        } else if (!is.null(py_ci)) {
          expect_match(r$method, "wilson")
          expect_equal(c(r$lower, r$upper), unlist(py_ci),
            tolerance = 1e-12, info = paste(f, metric, "wilson"))
        }
      }
    }
  })
}
