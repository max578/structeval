# Published Table 1 counts (Moldovan et al., medRxiv 2026.05.18.26353541) and
# the per-note counts they imply. Scalar fields are reconstructed exactly:
# every note contributes one of cor, inc, mis (field present) or spu, tn
# (field absent). List fields are spread one item per row, which keeps
# every pooled sum.

table1 <- data.frame(
  field = c(
    "First Primary Diagnosis", "First Primary Histology",
    "First Primary Overall Stage", "Treatment Drugs",
    "Biomarker Test Results"
  ),
  type = c("scalar", "scalar", "scalar", "list", "list"),
  cor = c(154, 124, 88, 352, 301), inc = c(3, 5, 0, 0, 0),
  mis = c(1, 1, 6, 40, 41), spu = c(2, 0, 8, 55, 16),
  tn = c(40, 70, 98, NA, NA),
  stringsAsFactors = FALSE
)

table1_published <- data.frame(
  field = rep(table1$field, c(4, 4, 4, 3, 3)),
  metric = c(
    rep(c("precision", "recall", "f1", "specificity"), 3),
    rep(c("precision", "recall", "f1"), 2)
  ),
  value = c(
    96.9, 97.5, 97.2, 95.2, 96.1, 95.4, 95.8, 100, 91.7, 93.6, 92.6, 92.5,
    86.5, 89.8, 88.1, 95.0, 88.0, 91.4
  ),
  stringsAsFactors = FALSE
)

expand_counts <- function(row) {
  outcome <- c("cor", "inc", "mis", "spu", "tn")
  times <- unlist(row[outcome])
  times[is.na(times)] <- 0
  lab <- rep(outcome, times)
  m <- vapply(outcome, function(o) as.integer(lab == o), integer(length(lab)))
  m <- matrix(m, ncol = 5L, dimnames = list(NULL, outcome))
  data.frame(
    note = as.character(seq_along(lab)), field = row$field, type = row$type,
    m, stringsAsFactors = FALSE
  )
}

table1_counts <- function() {
  do.call(rbind, lapply(seq_len(nrow(table1)), function(i) {
    expand_counts(table1[i, ])
  }))
}
