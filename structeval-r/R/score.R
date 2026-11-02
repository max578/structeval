# Counts per note and field ------------------------------------------------

#' Count agreements per note and field
#'
#' @description
#' Compares reference and model items for every note and field and returns
#' five counts per row: correct (`cor`), incorrect (`inc`), missed (`mis`),
#' spurious (`spu`) and true negative (`tn`).
#'
#' @details
#' Items are normalised before comparison: surrounding whitespace is trimmed,
#' case is folded, internal runs of whitespace collapse to one space, and
#' empty items and `"-"` items are dropped. Duplicate items count once.
#'
#' * **binary**: the reference and the prediction are `TRUE` or `FALSE`
#'   (also `t`/`f`, `yes`/`no`, `1`/`0`). A true positive is `cor`, a false
#'   positive `spu`, a false negative `mis` and a true negative `tn`.
#' * **scalar**: zero or one item on each side. Both empty gives `tn = 1`;
#'   reference only gives `mis = 1`; prediction only gives `spu = 1`; both
#'   present gives `cor = 1` when equal and `inc = 1` otherwise.
#' * **list**: `cor` is the size of the intersection, `mis` the reference
#'   items not predicted, `spu` the predicted items not in the reference;
#'   both empty gives `tn = 1`.
#'
#' A row whose reference is `NA` is unlabelled and skipped. A field with no
#' labelled row is reported by a warning and kept in the `fields` attribute,
#' so [metrics()] reports it as undefined rather than dropping it.
#'
#' @param data Long table with columns `note`, `field`, `type`, `truth` and
#'   `pred`, as returned by [read_oncoshot()]. `truth` and `pred` are lists
#'   of character (or logical) vectors; `character(0)` means no information
#'   and `NA` means unlabelled.
#' @param missing Policy for a missing prediction (`NA`, or a binary value
#'   that is not `TRUE`/`FALSE`). `"wrong"` scores a binary prediction as the
#'   wrong class and a scalar or list prediction as empty. `"skip"` drops the
#'   row.
#'
#' @returns A data frame with columns `note`, `field`, `type`, `cor`, `inc`,
#'   `mis`, `spu` and `tn` (integers), one row per scored note and field. The
#'   attribute `fields` lists every field and its type.
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
#' @export
score <- function(data, missing = c("wrong", "skip")) {
  missing <- match.arg(missing)
  .check_long(data)

  # Score each row ------------------------------------------------------------

  rows <- lapply(seq_len(nrow(data)), function(i) {
    .score_one(data$truth[[i]], data$pred[[i]], data$type[i], missing)
  })
  kept <- !vapply(rows, is.null, logical(1L))
  mat <- matrix(
    c(integer(0), unlist(rows[kept], use.names = FALSE)),
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

  # Fields with no labelled row -----------------------------------------------

  fields <- unique(data.frame(
    field = as.character(data$field), type = as.character(data$type),
    stringsAsFactors = FALSE
  ))
  rownames(fields) <- NULL
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
  out
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

#' Is a cell unlabelled?
#' @noRd
#' @keywords internal
.is_unlabelled <- function(x) {
  is.null(x) || (length(x) == 1L && is.na(x))
}

#' Normalise items: trim, case-fold, collapse whitespace, drop empties
#' @noRd
#' @keywords internal
.normalise <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- tolower(gsub("[[:space:]]+", " ", trimws(x)))
  unique(x[nzchar(x) & x != "-"])
}

#' Map a binary cell to TRUE, FALSE or NA
#' @noRd
#' @keywords internal
.as_flag <- function(x) {
  x <- .normalise(x)
  if (length(x) != 1L) {
    return(NA)
  }
  if (x %in% c("true", "t", "yes", "y", "1")) {
    return(TRUE)
  }
  if (x %in% c("false", "f", "no", "n", "0")) {
    return(FALSE)
  }
  NA
}

#' Counts for one note and field; NULL when the row is skipped
#' @noRd
#' @keywords internal
.score_one <- function(truth, pred, type, missing) {
  if (.is_unlabelled(truth)) {
    return(NULL)
  }

  # Binary ----------------------------------------------------------------

  if (type == "binary") {
    t_flag <- .as_flag(truth)
    if (is.na(t_flag)) {
      if (length(.normalise(truth)) == 0L) {
        return(NULL)
      }
      stop(
        "Binary reference values must be TRUE or FALSE; got '",
        paste(truth, collapse = ", "), "'.",
        call. = FALSE
      )
    }
    p_flag <- if (.is_unlabelled(pred)) NA else .as_flag(pred)
    if (is.na(p_flag)) {
      if (missing == "skip") {
        return(NULL)
      }
      p_flag <- !t_flag
    }
    return(c(
      cor = as.integer(t_flag && p_flag), inc = 0L,
      mis = as.integer(t_flag && !p_flag),
      spu = as.integer(!t_flag && p_flag),
      tn = as.integer(!t_flag && !p_flag)
    ))
  }

  # Scalar and list -------------------------------------------------------

  if (.is_unlabelled(pred)) {
    if (missing == "skip") {
      return(NULL)
    }
    pred <- character(0)
  }
  t_set <- .normalise(truth)
  p_set <- .normalise(pred)
  if (type == "scalar") {
    if (length(t_set) > 1L || length(p_set) > 1L) {
      stop(
        "A scalar field holds more than one item: '",
        paste(c(t_set, p_set), collapse = "', '"), "'.",
        call. = FALSE
      )
    }
    both <- length(t_set) == 1L && length(p_set) == 1L
    return(c(
      cor = as.integer(both && t_set == p_set),
      inc = as.integer(both && t_set != p_set),
      mis = as.integer(length(t_set) == 1L && length(p_set) == 0L),
      spu = as.integer(length(t_set) == 0L && length(p_set) == 1L),
      tn = as.integer(length(t_set) == 0L && length(p_set) == 0L)
    ))
  }
  c(
    cor = length(intersect(t_set, p_set)), inc = 0L,
    mis = length(setdiff(t_set, p_set)),
    spu = length(setdiff(p_set, t_set)),
    tn = as.integer(length(t_set) == 0L && length(p_set) == 0L)
  )
}
