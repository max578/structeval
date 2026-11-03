# Counts per note and field ------------------------------------------------

#' Count agreements per note and field
#'
#' @description
#' Compares reference and model items for every note and field and returns
#' five counts per row: correct (`cor`), incorrect (`inc`), missed (`mis`),
#' spurious (`spu`) and true negative (`tn`).
#'
#' @details
#' **Normalisation.** Every item passes through the same steps, in this
#' order, on both sides:
#'
#' 1. Whitespace: every Unicode space character (tab, line feed, vertical
#'    tab, form feed, carriage return, space, U+0085, U+00A0, U+1680,
#'    U+2000 to U+200A, U+2028, U+2029, U+202F, U+205F, U+3000) is
#'    trimmed at both ends, and internal runs collapse to one space.
#' 2. Case: the simple one-to-one Unicode lowercase mapping ([tolower()]),
#'    so the German sharp s does not match `"ss"`. This needs a UTF-8 locale when
#'    an item is not ASCII; outside one, `score()` stops rather than fold
#'    differently. Composed and decomposed forms of an accented letter are
#'    different strings (no Unicode normalisation).
#' 3. Numbers: an item that is a plain decimal, an optional minus sign,
#'    an integer part without leading zeros and an optional fraction
#'    (`-?(0|[1-9][0-9]*)(\.[0-9]+)?`), is rewritten in one canonical form
#'    (trailing fractional zeros dropped, `"-0"` read as `"0"`), so `"4"`,
#'    `"4.0"` and `"4.00"` match. Everything else stays text and is compared
#'    as text: codes with leading zeros (`"007"`, `"02"`), exponent forms
#'    (`"1e3"`), `"+4"`, `".5"`, `"inf"` and `"NaN"`. A `"NaN"` item is a
#'    value, never a missing one.
#' 4. Empty values: an item that is then empty, `"-"`, an en dash or an em
#'    dash, or one of `null_tokens` is dropped. The default tokens are
#'    `"null"`, `"none"`, `"na"` and `"n/a"`, matched after case folding, so
#'    `"NULL"`, `"None"`, `"NA"` and `"N/A"` are empty values. The base paper
#'    defines a missed value as one where "the model output is empty or
#'    null". `"unknown"` and `"NaN"` stay values. The rule is the same for a
#'    long table built by hand and for one from [read_oncoshot()], which
#'    keeps the text `"NA"` as written (only a blank cell becomes `NA`).
#' 5. Duplicate items count once.
#'
#' **Outcomes by field type.**
#'
#' * **binary**: the reference and the prediction are `TRUE` or `FALSE`
#'   (also `t`/`f`, `yes`/`no`, `y`/`n`, `1`/`0`, after normalisation, so
#'   `"1.0"` is `TRUE`). A true positive is `cor`, a false positive `spu`, a
#'   false negative `mis` and a true negative `tn`. Any other value stops
#'   with the distinct values listed, in the reference and in the
#'   prediction alike; recode them before scoring.
#' * **scalar**: zero or one item on each side. Both empty gives `tn = 1`;
#'   reference only gives `mis = 1`; prediction only gives `spu = 1`; both
#'   present gives `cor = 1` when equal and `inc = 1` otherwise.
#' * **list**: `cor` is the size of the intersection, `mis` the reference
#'   items not predicted, `spu` the predicted items not in the reference.
#'   `tn = 1` records a note in which both sides are empty (the
#'   "noncommittal" case of MUC-5). No list metric uses it: list fields
#'   have no specificity and no accuracy, as in the base paper.
#'
#' **Unlabelled and empty references.** A reference that is `NULL`, `NA`,
#' or only blank strings (after the whitespace rule) is unlabelled for every
#' type and the row is skipped. A reference that is empty after
#' normalisation (`character(0)`, `"-"`, a dash or a null token) holds no
#' information: for scalar and list fields it is scored (`tn` or `spu`);
#' for binary fields it has no class, so the row is skipped and reported
#' by a message.
#'
#' **Missing predictions.** A prediction that is `NA` or blank, and for
#' binary fields one that is empty after normalisation, is missing. Under
#' `missing = "wrong"` a missing binary prediction is scored as the wrong
#' class, while a missing scalar or list prediction is scored as an empty
#' answer: a missed value when the reference has one, and a true negative
#' when the reference is empty. Under `missing = "skip"` the row is dropped
#' and a message gives the count. A scalar or list prediction of `"-"` or a
#' null token is an empty answer under both policies.
#'
#' The attribute `coverage` records, per field, the rows read
#' (`rows`), skipped as unlabelled (`unlabelled`), skipped as binary
#' references with no information (`no_information`), with a missing
#' prediction (`pred_missing`), dropped under `missing = "skip"`
#' (`skipped_missing`) and scored (`scored`). [metrics()] carries it over.
#' A field with no scored row is reported by a warning and kept in the
#' `fields` attribute, so [metrics()] reports it as undefined rather than
#' dropping it. A note and field pair that appears twice stops `score()`.
#'
#' @param data Long table with columns `note`, `field`, `type`, `truth` and
#'   `pred`, as returned by [read_oncoshot()]. `truth` and `pred` are lists
#'   of character (or logical) vectors; `character(0)` means no information
#'   and `NA` means unlabelled.
#' @param missing Policy for a missing prediction, `"wrong"` or `"skip"`,
#'   as described in Details.
#' @param null_tokens Items that mean an empty value, compared after
#'   normalisation. `character(0)` keeps every such item as a value.
#'
#' @returns A data frame with columns `note`, `field`, `type`, `cor`, `inc`,
#'   `mis`, `spu` and `tn` (integers), one row per scored note and field. The
#'   attribute `fields` lists every field and its type, and the attribute
#'   `coverage` gives the row counts described in Details.
#'
#' @references
#' Chinchor, N., & Sundheim, B. (1993). MUC-5 evaluation metrics. In
#' *Proceedings of the Fifth Message Understanding Conference (MUC-5)*
#' (pp. 69-78). Association for Computational Linguistics.
#' \doi{10.3115/1072017.1072026}
#'
#' @seealso [metrics()] for estimates and intervals.
#'
#' @examples
#' long <- data.frame(
#'   note = c("1", "2", "3"),
#'   field = "Drugs",
#'   type = "list"
#' )
#' long$truth <- list(c("XELOX", "Bevacizumab"), character(0), "FOLFOX")
#' long$pred <- list(c("xelox"), character(0), c("FOLFOX", "Cetuximab"))
#' score(long)
#' attr(score(long), "coverage")
#' @export
score <- function(data, missing = c("wrong", "skip"),
                  null_tokens = c("null", "none", "na", "n/a")) {
  missing <- match.arg(missing)
  .check_long(data)
  if (!is.character(null_tokens) || anyNA(null_tokens)) {
    stop("`null_tokens` must be a character vector without NA.", call. = FALSE)
  }
  ctx <- list(
    missing = missing, utf8 = isTRUE(l10n_info()[["UTF-8"]]),
    null_tokens = character(0)
  )
  ctx$null_tokens <- .normalise(null_tokens, ctx)

  # Score each row ------------------------------------------------------------

  rows <- lapply(seq_len(nrow(data)), function(i) {
    .score_one(data$truth[[i]], data$pred[[i]], data$type[i], ctx)
  })
  status <- vapply(rows, `[[`, character(1L), "status")
  for (side in c("truth", "pred")) {
    bad <- status == paste0("bad_", side)
    if (any(bad)) {
      values <- unique(unlist(lapply(rows[bad], `[[`, "value")))
      stop(
        "Binary ", if (side == "truth") "reference" else "prediction",
        " values must be TRUE or FALSE (or t/f, yes/no, y/n, 1/0); got ",
        paste0("'", values, "'", collapse = ", "),
        ". Recode them before scoring.",
        call. = FALSE
      )
    }
  }
  kept <- status == "scored"
  mat <- matrix(
    c(integer(0), unlist(lapply(rows[kept], `[[`, "counts"),
      use.names = FALSE)),
    ncol = 5L, byrow = TRUE,
    dimnames = list(NULL, c("cor", "inc", "mis", "spu", "tn"))
  )
  out <- data.frame(
    note = as.character(data$note[kept]),
    field = as.character(data$field[kept]),
    type = as.character(data$type[kept]),
    stringsAsFactors = FALSE
  )
  out <- cbind(out, as.data.frame(mat))

  # Coverage per field ----------------------------------------------------------

  fields <- unique(data.frame(
    field = as.character(data$field), type = as.character(data$type),
    stringsAsFactors = FALSE
  ))
  rownames(fields) <- NULL
  pred_missing <- vapply(rows, `[[`, logical(1L), "pred_missing")
  by_field <- factor(as.character(data$field), levels = fields$field)
  tally <- function(x) as.integer(tapply(x, by_field, sum, default = 0L))
  coverage <- data.frame(
    field = fields$field, type = fields$type,
    rows = tally(rep(1L, nrow(data))),
    unlabelled = tally(status == "unlabelled"),
    no_information = tally(status == "no_information"),
    pred_missing = tally(pred_missing),
    skipped_missing = tally(status == "skipped_missing"),
    scored = tally(kept),
    stringsAsFactors = FALSE
  )
  .report_coverage(coverage)
  empty <- setdiff(fields$field, out$field)
  if (length(empty) > 0L) {
    warning(
      "No scored notes for field(s): ",
      paste0("'", empty, "'", collapse = ", "),
      " (every row unlabelled or skipped). Their metrics are undefined.",
      call. = FALSE
    )
  }
  attr(out, "fields") <- fields
  attr(out, "coverage") <- coverage
  out
}

#' Messages for rows that score() leaves out
#' @noRd
#' @keywords internal
.report_coverage <- function(coverage) {
  for (i in seq_len(nrow(coverage))) {
    f <- coverage$field[i]
    if (coverage$no_information[i] > 0L) {
      message(
        "Field '", f, "': ", coverage$no_information[i],
        " row(s) with a binary reference holding no information ('-' or ",
        "a null token) skipped."
      )
    }
    if (coverage$skipped_missing[i] > 0L) {
      message(
        "Field '", f, "': ", coverage$skipped_missing[i],
        " row(s) with a missing prediction skipped (missing = \"skip\")."
      )
    }
  }
  invisible(NULL)
}

#' Validate the long table
#' @noRd
#' @keywords internal
.check_long <- function(data) {
  cols <- c("note", "field", "type", "truth", "pred")
  absent <- setdiff(cols, names(data))
  if (!is.data.frame(data) || length(absent) > 0L) {
    stop(
      "`data` must be a data frame with columns ",
      paste0("'", cols, "'", collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.list(data$truth) || !is.list(data$pred)) {
    stop("`truth` and `pred` must be list columns.", call. = FALSE)
  }
  bad <- setdiff(unique(data$type), c("binary", "scalar", "list"))
  if (length(bad) > 0L) {
    stop(
      "Field types must be 'binary', 'scalar' or 'list'; got ",
      paste0("'", bad, "'", collapse = ", "), ".",
      call. = FALSE
    )
  }
  dup <- duplicated(data.frame(
    note = as.character(data$note), field = as.character(data$field)
  ))
  if (any(dup)) {
    stop(
      "Each note and field pair must appear once; repeated for note(s) ",
      paste0("'", unique(as.character(data$note[dup])), "'", collapse = ", "),
      " (a pair given more than once would be counted more than once).",
      call. = FALSE
    )
  }
  n_types <- tapply(data$type, data$field, function(x) length(unique(x)))
  if (any(n_types > 1L)) {
    stop(
      "Each field must have one type; mixed types for: ",
      paste0("'", names(n_types)[n_types > 1L], "'", collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Unicode space characters trimmed and collapsed by the normaliser
#' @noRd
#' @keywords internal
.space_class <- paste0(
  "[\t-\r \u0085\u00a0\u1680\u2000-\u200a",
  "\u2028\u2029\u202f\u205f\u3000]"
)

#' Characters that mean no information, like "-"
#' @noRd
#' @keywords internal
.dashes <- c("-", "\u2013", "\u2014")

#' Trim and collapse every Unicode space
#' @noRd
#' @keywords internal
.squish <- function(x) {
  x <- enc2utf8(as.character(x))
  x <- gsub(paste0(.space_class, "+"), " ", x, perl = TRUE)
  gsub("^ | $", "", x, perl = TRUE)
}

#' Lower case, refusing to fold non-ASCII text outside a UTF-8 locale
#' @noRd
#' @keywords internal
.fold_case <- function(x, utf8 = isTRUE(l10n_info()[["UTF-8"]])) {
  if (!utf8) {
    wide <- !is.na(x) & grepl("[^\\x01-\\x7F]", x, perl = TRUE)
    if (any(wide)) {
      stop(
        "Items with non-ASCII characters (e.g. '", x[wide][1L], "') need a ",
        "UTF-8 locale to be case-folded the same on every machine; the ",
        "current locale is ", Sys.getlocale("LC_CTYPE"), ". Set one, for ",
        "example Sys.setlocale(\"LC_CTYPE\", \"en_US.UTF-8\").",
        call. = FALSE
      )
    }
  }
  tolower(x)
}

#' Canonical form of plain decimals; every other item is returned unchanged
#' @noRd
#' @keywords internal
.canon_number <- function(x) {
  num <- grepl("^-?(0|[1-9][0-9]*)(\\.[0-9]+)?$", x, perl = TRUE)
  y <- x[num]
  y <- sub("(\\.[0-9]*[1-9])0+$", "\\1", y, perl = TRUE)
  y <- sub("\\.0+$", "", y, perl = TRUE)
  y[y == "-0"] <- "0"
  x[num] <- y
  x
}

#' Is a cell unlabelled (NULL, NA, or only blank strings)?
#' @noRd
#' @keywords internal
.is_unlabelled <- function(x) {
  if (is.null(x)) {
    return(TRUE)
  }
  if (length(x) == 0L) {
    return(FALSE)
  }
  x <- as.character(x)
  all(is.na(x) | !nzchar(.squish(x)))
}

#' Normalise items: whitespace, case, numbers, empty values, duplicates
#' @noRd
#' @keywords internal
.normalise <- function(x, ctx) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- .canon_number(.fold_case(.squish(x), ctx$utf8))
  unique(x[nzchar(x) & !x %in% c(.dashes, ctx$null_tokens)])
}

#' State of a binary cell: TRUE, FALSE, empty or unrecognised
#' @noRd
#' @keywords internal
.binary_state <- function(x, ctx) {
  x <- .normalise(x, ctx)
  if (length(x) == 0L) {
    return("empty")
  }
  if (length(x) == 1L && x %in% c("true", "t", "yes", "y", "1")) {
    return("true")
  }
  if (length(x) == 1L && x %in% c("false", "f", "no", "n", "0")) {
    return("false")
  }
  "other"
}

#' Counts and status for one note and field
#' @noRd
#' @keywords internal
.score_one <- function(truth, pred, type, ctx) {
  row <- function(status, counts = NULL, pred_missing = FALSE,
                  value = NULL) {
    list(
      status = status, counts = counts, pred_missing = pred_missing,
      value = value
    )
  }
  if (.is_unlabelled(truth)) {
    return(row("unlabelled"))
  }

  # Binary ----------------------------------------------------------------

  if (type == "binary") {
    t_state <- .binary_state(truth, ctx)
    if (t_state == "empty") {
      return(row("no_information"))
    }
    if (t_state == "other") {
      return(row("bad_truth", value = paste(truth, collapse = ", ")))
    }
    t_flag <- t_state == "true"
    p_state <- if (.is_unlabelled(pred)) "empty" else .binary_state(pred, ctx)
    if (p_state == "other") {
      return(row("bad_pred", value = paste(pred, collapse = ", ")))
    }
    if (p_state == "empty") {
      if (ctx$missing == "skip") {
        return(row("skipped_missing", pred_missing = TRUE))
      }
      p_flag <- !t_flag
    } else {
      p_flag <- p_state == "true"
    }
    return(row("scored", pred_missing = p_state == "empty", counts = c(
      cor = as.integer(t_flag && p_flag), inc = 0L,
      mis = as.integer(t_flag && !p_flag),
      spu = as.integer(!t_flag && p_flag),
      tn = as.integer(!t_flag && !p_flag)
    )))
  }

  # Scalar and list -------------------------------------------------------

  pred_missing <- .is_unlabelled(pred)
  if (pred_missing) {
    if (ctx$missing == "skip") {
      return(row("skipped_missing", pred_missing = TRUE))
    }
    pred <- character(0)
  }
  t_set <- .normalise(truth, ctx)
  p_set <- .normalise(pred, ctx)
  if (type == "scalar") {
    if (length(t_set) > 1L || length(p_set) > 1L) {
      stop(
        "A scalar field holds more than one item: '",
        paste(c(t_set, p_set), collapse = "', '"), "'.",
        call. = FALSE
      )
    }
    both <- length(t_set) == 1L && length(p_set) == 1L
    return(row("scored", pred_missing = pred_missing, counts = c(
      cor = as.integer(both && t_set == p_set),
      inc = as.integer(both && t_set != p_set),
      mis = as.integer(length(t_set) == 1L && length(p_set) == 0L),
      spu = as.integer(length(t_set) == 0L && length(p_set) == 1L),
      tn = as.integer(length(t_set) == 0L && length(p_set) == 0L)
    )))
  }
  row("scored", pred_missing = pred_missing, counts = c(
    cor = length(intersect(t_set, p_set)), inc = 0L,
    mis = length(setdiff(t_set, p_set)),
    spu = length(setdiff(p_set, t_set)),
    tn = as.integer(length(t_set) == 0L && length(p_set) == 0L)
  ))
}
