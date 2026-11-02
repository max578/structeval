# make_crosscheck.R -- writes crosscheck.csv, the shared fixture for the
# R and Python implementations (Oncoshot wide format, synthetic).
#
# Five fields: one binary, two scalar, two list. Rows 1-22 cover each edge
# case of the counting rules; rows 23-60 are drawn at random from small
# vocabularies. "Radiotherapy Sites" is absent in every note.
#
# Run from this directory: Rscript make_crosscheck.R

set.seed(20261029)

edge <- list(
  # Binary: TP, FN, FP, TN, unlabelled, missing prediction (both classes),
  # "-" as a prediction.
  c("TRUE", "TRUE"), c("TRUE", "FALSE"), c("FALSE", "TRUE"),
  c("FALSE", "FALSE"), c("", "TRUE"), c("TRUE", ""), c("FALSE", ""),
  c("TRUE", "-")
)

dx <- list(
  c("Colorectal Cancer", "Colorectal Cancer"),
  c("Colorectal Cancer", "colorectal cancer"),
  c("  Breast Cancer ", "Breast Cancer"),
  c("Breast Cancer", "Lung Cancer"),
  c("Lung Cancer", "-"),
  c("-", "Lung Cancer"),
  c("-", "-"),
  c("", "Lung Cancer"),
  c("Gastric Cancer", "")
)
hist <- list(
  c("Adenocarcinoma", "ADENOCARCINOMA"),
  c("Adenocarcinoma", "Squamous cell carcinoma"),
  c("-", "-"),
  c("Squamous cell carcinoma", "-"),
  c("-", "Adenocarcinoma"),
  c("", ""),
  c("Small cell carcinoma", "Small cell carcinoma")
)
drugs <- list(
  c("['XELOX']", "['XELOX']"),
  c("['Docetaxel', 'Gemcitabine']", "['gemcitabine', 'DOCETAXEL']"),
  c("['XELOX']", "['XELOX', 'xelox']"),
  c("['Docetaxel', 'Gemcitabine']", "['Docetaxel', 'XELOX']"),
  c("[]", "[]"),
  c("[]", "['Cisplatin']"),
  c("['Cisplatin']", "[]"),
  c("-", "['Cisplatin']"),
  c("", "['Cisplatin']"),
  c("['Oxaliplatin', 'Bevacizumab']", ""),
  c("[' Oxaliplatin ']", "['Oxaliplatin']"),
  c("['Pembrolizumab', 'Carboplatin', 'Pemetrexed']",
    "['Pembrolizumab', 'Carboplatin', 'Pemetrexed', 'Nivolumab']")
)

n_notes <- 60L
pad <- function(x, n) c(x, vector("list", n - length(x)))
edge <- pad(edge, n_notes)
dx <- pad(dx, n_notes)
hist <- pad(hist, n_notes)
drugs <- pad(drugs, n_notes)

# Random rows ----------------------------------------------------------------

dx_vocab <- c("Colorectal Cancer", "Breast Cancer", "Lung Cancer",
  "Gastric Cancer", "-")
hist_vocab <- c("Adenocarcinoma", "Squamous cell carcinoma",
  "Small cell carcinoma", "-")
drug_vocab <- c("XELOX", "Docetaxel", "Gemcitabine", "Cisplatin",
  "Oxaliplatin", "Bevacizumab", "Pembrolizumab", "Carboplatin")

draw_scalar <- function(vocab) {
  truth <- sample(vocab, 1L, prob = c(rep(0.8 / (length(vocab) - 1L),
    length(vocab) - 1L), 0.2))
  u <- runif(1L)
  pred <- if (u < 0.75) truth else sample(vocab, 1L)
  c(truth, pred)
}
as_list_cell <- function(x) {
  if (length(x) == 0L) "[]" else paste0("[", paste0("'", x, "'",
    collapse = ", "), "]")
}
draw_list <- function() {
  truth <- sample(drug_vocab, rpois(1L, 1.5))
  keep <- truth[runif(length(truth)) < 0.8]
  extra <- if (runif(1L) < 0.3) sample(drug_vocab, 1L) else character(0)
  c(as_list_cell(truth), as_list_cell(unique(c(keep, extra))))
}
for (i in seq_len(n_notes)) {
  if (is.null(edge[[i]])) {
    t <- sample(c("TRUE", "FALSE"), 1L)
    edge[[i]] <- c(t, if (runif(1L) < 0.85) t else setdiff(c("TRUE",
      "FALSE"), t))
  }
  if (is.null(dx[[i]])) dx[[i]] <- draw_scalar(dx_vocab)
  if (is.null(hist[[i]])) hist[[i]] <- draw_scalar(hist_vocab)
  if (is.null(drugs[[i]])) drugs[[i]] <- draw_list()
}
sites <- lapply(seq_len(n_notes), function(i) {
  switch(i %% 4L + 1L, c("-", "-"), c("[]", "[]"), c("-", "[]"), c("", "-"))
})

# Write ----------------------------------------------------------------------

col <- function(x, j) vapply(x, `[`, character(1L), j)
tab <- data.frame(
  "Patient ID" = as.character(seq_len(n_notes)),
  "Has metastasis" = col(edge, 1L), "Res: Has metastasis" = col(edge, 2L),
  "First Primary Diagnosis" = col(dx, 1L),
  "Res: First Primary Diagnosis" = col(dx, 2L),
  "First Primary Histology" = col(hist, 1L),
  "Res: First Primary Histology" = col(hist, 2L),
  "Treatment Drugs" = col(drugs, 1L), "Res: Treatment Drugs" = col(drugs, 2L),
  "Radiotherapy Sites" = col(sites, 1L),
  "Res: Radiotherapy Sites" = col(sites, 2L),
  check.names = FALSE
)
quote_cell <- function(x) {
  needs <- grepl("[,\"]", x) | grepl("^ | $", x)
  ifelse(needs, paste0("\"", gsub("\"", "\"\"", x), "\""), x)
}
lines <- c(
  paste(quote_cell(names(tab)), collapse = ","),
  apply(tab, 1L, function(r) paste(quote_cell(r), collapse = ","))
)
writeLines(lines, "crosscheck.csv", useBytes = TRUE)
