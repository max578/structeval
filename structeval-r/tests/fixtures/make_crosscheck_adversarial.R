# make_crosscheck_adversarial.R -- writes crosscheck_adversarial.csv, the
# second shared fixture for the R and Python implementations. Each row holds
# an input that the two implementations scored differently, or wrongly,
# before 0.1.1 (adversarial review SC-02 to SC-11, K2): binary spellings,
# null tokens (null, none, NA, N/A) and dashes, Unicode white space and case, numbers written as
# text, escaped and bare bracket lists. Every cell is quoted.
#
# Run from this directory: Rscript make_crosscheck_adversarial.R

rows <- list(
  c("1", "TRUE", "Yes", "Stage  IIA", "stage iia", "4", "4.0",
    "['O\\'Brien']", "[\"O'Brien\"]"),
  c("2", "FALSE", "no", " Adeno", "Adeno", "12", "12.00",
    "['A', 'B']", "[a, b]"),
  c("3", "TRUE", "1.0", "Étape II", "étape ii", "007", "7",
    "['4', '5.0']", "[4.0, 5]"),
  c("4", "FALSE", "0.0", "Straße", "STRASSE", "1e3", "1000",
    "['None']", "[]"),
  c("5", "yes", "None", "Adeno", "null", "inf", "Infinity",
    "['', 'A']", "['A']"),
  c("6", "no", "-", "Adeno", "None", "NaN", "NaN", "[]", "['null']"),
  c("7", "–", "TRUE", "None", "-", "NaN", "nan", "-", "–"),
  c("8", "null", "FALSE", "N/A", "n/a", "1_000", "1000", "['X']", ""),
  c("9", "TRUE", "", "NA", "Adeno", ".5", "0.5",
    "[' Oxaliplatin ']", "['oxaliplatin']"),
  c("10", "", "TRUE", "–", "", "-0", "0", "['A','A']", "['a']"),
  c("11", "FALSE", "n", "-", "—", "0.5", "0.50", "[\"x, y\"]",
    "['X, Y']"),
  c("12", "TRUE", "t", "ΟΔΟΣ", "οδοσ",
    "100", "100.", "[A]", "['B']"),
  c("13", "FALSE", "f", "Adeno carcinoma", "Adeno carcinoma", "+4", "4",
    "['a\\\\b']", "['A\\\\B']"),
  c("14", "TRUE", "TRUE", "\vAdeno", "Adeno", "02", "2",
    "['Cisplatin', 'XELOX']", "['cisplatin']"),
  c("15", "FALSE", "TRUE", "  ", "x", "4.50", "4.5", "", "['A']"),
  c("16", "TRUE", "FALSE", "unknown", "Adeno", "-", "4", "[]", "[]"),
  c("17", "N/A", "TRUE", "Adeno", "na", "NA", "-", "['NA', 'A']",
    "['a', 'n/a']"),
  c("18", "TRUE", "nA", "N/A", "Adeno", "n/a", "NA", "[N/A]", "['X']")
)
header <- c("Patient ID", "Met", "Res: Met", "Hist", "Res: Hist", "Size",
  "Res: Size", "Drugs", "Res: Drugs")
quote_all <- function(x) paste0("\"", gsub("\"", "\"\"", x), "\"")
lines <- vapply(c(list(header), rows), function(r) {
  paste(quote_all(enc2utf8(r)), collapse = ",")
}, character(1L))
writeLines(enc2utf8(lines), "crosscheck_adversarial.csv", useBytes = TRUE)
