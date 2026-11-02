# Reader for the Oncoshot wide CSV ------------------------------------------

#' Read an Oncoshot wide CSV into the long table
#'
#' @description
#' Reads the wide evaluation file used by the Oncoshot pipeline: one row per
#' note, a column `<field>` with the reference value and a column
#' `Res: <field>` with the model output, for each field. The result is the
#' long table that [score()] consumes, with one row per note and field.
#'
#' @details
#' Cell conventions follow the source format. A blank cell is unlabelled and
#' becomes `NA` (the row is skipped by [score()]). The value `"-"` means the
#' note holds no information and becomes `character(0)`. List cells may be
#' comma-separated (`"a, b"`) or written as a bracketed list of quoted items
#' (`"['a', 'b']"`); `"[]"` is the empty list. Scalar and binary cells hold
#' one value. Normalisation (trimming, case folding, collapsing internal
#' whitespace) is applied by [score()], so the long table keeps the values as
#' written.
#'
#' @param file Path to a UTF-8 CSV file. A byte-order mark is accepted.
#' @param types Named character vector giving the structure of each field to
#'   read, one of `"binary"`, `"scalar"` or `"list"`. The names are the field
#'   names as they appear in the reference columns. Fields not named here are
#'   ignored.
#' @param id Name of the note identifier column.
#' @param pred_prefix Prefix of the prediction columns.
#'
#' @returns A data frame with columns `note` (character), `field`
#'   (character), `type` (character), and the list columns `truth` and `pred`,
#'   each element a character vector of items or `NA` when unlabelled.
#'
#' @seealso [score()] to count agreements in the long table.
#'
#' @examples
#' csv <- tempfile(fileext = ".csv")
#' writeLines(c(
#'   "Patient ID,Stage,Res: Stage,Drugs,Res: Drugs",
#'   "1,II,II,\"['XELOX']\",\"['XELOX', 'Bevacizumab']\"",
#'   "2,-,III,[],[]",
#'   "3,,IV,\"['FOLFOX']\",-"
#' ), csv)
#' long <- read_oncoshot(csv, types = c(Stage = "scalar", Drugs = "list"))
#' long
#' score(long)
#' @export
read_oncoshot <- function(file, types, id = "Patient ID",
                          pred_prefix = "Res: ") {
  .check_types(types)
  raw <- utils::read.csv(
    file,
    check.names = FALSE, colClasses = "character", na.strings = character(0),
    fileEncoding = "UTF-8-BOM", strip.white = FALSE
  )

  # Columns ---------------------------------------------------------------

  fields <- names(types)
  needed <- c(id, fields, paste0(pred_prefix, fields))
  absent <- setdiff(needed, names(raw))
  if (length(absent) > 0L) {
    stop(
      "Columns not found in `file`: ",
      paste0("'", absent, "'", collapse = ", "), ".",
      call. = FALSE
    )
  }
  notes <- trimws(raw[[id]])
  if (anyDuplicated(notes) > 0L) {
    stop("Note identifiers in column '", id, "' are not unique.", call. = FALSE)
  }

  # Long table --------------------------------------------------------------

  per_field <- lapply(fields, function(f) {
    out <- data.frame(
      note = notes, field = f, type = unname(types[f]),
      stringsAsFactors = FALSE
    )
    out$truth <- lapply(raw[[f]], .parse_cell, type = types[[f]])
    out$pred <- lapply(raw[[paste0(pred_prefix, f)]], .parse_cell,
      type = types[[f]]
    )
    out
  })
  out <- do.call(rbind, per_field)
  rownames(out) <- NULL
  out
}

#' Parse one wide-format cell into items
#' @noRd
#' @keywords internal
.parse_cell <- function(x, type) {
  x <- trimws(x)
  if (is.na(x) || !nzchar(x)) {
    return(NA_character_)
  }
  if (x == "-") {
    return(character(0))
  }
  if (startsWith(x, "[") && endsWith(x, "]")) {
    return(.parse_bracket_list(x))
  }
  if (type == "list") {
    return(trimws(strsplit(x, ",", fixed = TRUE)[[1L]]))
  }
  x
}

#' Parse a bracketed list of quoted items
#' @noRd
#' @keywords internal
.parse_bracket_list <- function(x) {
  inner <- trimws(substr(x, 2L, nchar(x) - 1L))
  if (!nzchar(inner)) {
    return(character(0))
  }
  pattern <- "'(?:[^'\\\\]|\\\\.)*'|\"(?:[^\"\\\\]|\\\\.)*\""
  tokens <- regmatches(inner, gregexpr(pattern, inner, perl = TRUE))[[1L]]
  if (length(tokens) == 0L) {
    return(trimws(strsplit(inner, ",", fixed = TRUE)[[1L]]))
  }
  substr(tokens, 2L, nchar(tokens) - 1L)
}

#' Validate a named field-type vector
#' @noRd
#' @keywords internal
.check_types <- function(types) {
  ok <- is.character(types) && length(types) > 0L && !is.null(names(types)) &&
    all(nzchar(names(types))) && !anyDuplicated(names(types))
  if (!ok) {
    stop(
      "`types` must be a named character vector with one unique name per ",
      "field.",
      call. = FALSE
    )
  }
  bad <- setdiff(unique(types), c("binary", "scalar", "list"))
  if (length(bad) > 0L) {
    stop(
      "Field types must be 'binary', 'scalar' or 'list'; got ",
      paste0("'", bad, "'", collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
