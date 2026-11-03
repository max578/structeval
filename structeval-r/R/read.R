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
#' Cell conventions follow the source format. A blank cell (empty or only
#' white space) is unlabelled and becomes `NA` (the row is skipped by
#' [score()]). The value `"-"`, or an en or em dash, means the note holds no
#' information and becomes `character(0)`. List cells may be
#' comma-separated (`"a, b"`) or written as a bracketed list (`"['a', 'b']"`,
#' the syntax of a Python list of strings); `"[]"` is the empty list. In a
#' bracketed list either every item is quoted, with `'` or `"` and the
#' escapes `\'`, `\"` and `\\`, or no item is (`"[a, b]"`). Blank items are
#' dropped. Anything else, such as a mix of quoted and bare items or text
#' after the closing bracket, stops with the cell named. A bracketed cell
#' in a scalar or binary column is read the same way. Normalisation
#' (white space, case, numbers, null tokens) is applied by [score()], so the
#' long table keeps the values as written.
#'
#' The file is checked before it is parsed, and each check stops with a
#' named error: every line must be valid UTF-8; every record must have as
#' many fields as the header (an unquoted comma in a cell would otherwise
#' shift the columns); a warning from the CSV parser, such as an
#' unterminated quote, stops the read; the identifier, reference and
#' prediction columns must each appear once; note identifiers must be
#' non-empty and unique.
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
  .check_csv(file)
  raw <- withCallingHandlers(
    utils::read.csv(
      file,
      check.names = FALSE, colClasses = "character",
      na.strings = character(0), fileEncoding = "UTF-8-BOM",
      strip.white = FALSE
    ),
    warning = function(w) {
      stop("`file` could not be read cleanly: ", conditionMessage(w),
        call. = FALSE
      )
    }
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
  repeated <- intersect(needed, names(raw)[duplicated(names(raw))])
  if (length(repeated) > 0L) {
    stop(
      "Columns appear more than once in `file`: ",
      paste0("'", repeated, "'", collapse = ", "), ".",
      call. = FALSE
    )
  }
  notes <- .trim_space(raw[[id]])
  if (any(!nzchar(notes))) {
    stop(
      "Note identifiers in column '", id, "' are empty in data row(s) ",
      paste(which(!nzchar(notes)), collapse = ", "), ".",
      call. = FALSE
    )
  }
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

#' Check a CSV file before parsing: UTF-8 and one field count
#' @noRd
#' @keywords internal
.check_csv <- function(file) {
  lines <- readLines(file, warn = FALSE)
  bad <- which(!validUTF8(lines))
  if (length(bad) > 0L) {
    stop(
      "`file` is not valid UTF-8 on line(s) ",
      paste(utils::head(bad, 10L), collapse = ", "), ".",
      call. = FALSE
    )
  }
  # quote = "\"" as in read.csv(); the default would treat ' as a quote.
  n <- utils::count.fields(
    file, sep = ",", quote = "\"", comment.char = "",
    blank.lines.skip = FALSE
  )
  bad <- which(!is.na(n) & n != 0L & n != n[1L])
  if (length(bad) > 0L) {
    stop(
      "Every record must have ", n[1L], " fields, as the header has; ",
      paste0("line ", utils::head(bad, 10L), " has ", n[utils::head(bad, 10L)],
        collapse = ", "
      ),
      ". Quote cells that contain a comma.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Trim every Unicode space at both ends
#' @noRd
#' @keywords internal
.trim_space <- function(x) {
  x <- enc2utf8(as.character(x))
  gsub(paste0("^", .space_class, "+|", .space_class, "+$"), "", x,
    perl = TRUE
  )
}

#' Parse one wide-format cell into items
#' @noRd
#' @keywords internal
.parse_cell <- function(x, type) {
  x <- .trim_space(x)
  if (is.na(x) || !nzchar(x)) {
    return(NA_character_)
  }
  if (x %in% .dashes) {
    return(character(0))
  }
  if (startsWith(x, "[")) {
    return(.parse_bracket_list(x))
  }
  if (type == "list") {
    items <- .trim_space(strsplit(x, ",", fixed = TRUE)[[1L]])
    return(items[nzchar(items)])
  }
  x
}

#' Parse a bracketed list: all items quoted, or none
#' @noRd
#' @keywords internal
.parse_bracket_list <- function(x) {
  fail <- function(why) {
    stop("Malformed bracket list '", x, "': ", why, ".", call. = FALSE)
  }
  if (!endsWith(x, "]")) {
    fail("text after the closing bracket, or no closing bracket")
  }
  inner <- .trim_space(substr(x, 2L, nchar(x) - 1L))
  if (!nzchar(inner)) {
    return(character(0))
  }
  pattern <- "'(?:[^'\\\\]|\\\\.)*'|\"(?:[^\"\\\\]|\\\\.)*\""
  hits <- gregexpr(pattern, inner, perl = TRUE)
  tokens <- regmatches(inner, hits)[[1L]]
  if (length(tokens) == 0L) {
    if (grepl("['\"\\[\\]\\\\]", inner, perl = TRUE)) {
      fail("unbalanced quotes, brackets or backslashes in a bare list")
    }
    items <- .trim_space(strsplit(inner, ",", fixed = TRUE)[[1L]])
    return(items[nzchar(items)])
  }
  shape <- inner
  regmatches(shape, hits) <- list(rep("Q", length(tokens)))
  shape <- gsub(paste0(.space_class, "+"), "", shape, perl = TRUE)
  if (!grepl("^Q(,Q)*$", shape, perl = TRUE)) {
    fail("every item must be quoted and separated by one comma")
  }
  body <- substr(tokens, 2L, nchar(tokens) - 1L)
  escapes <- unlist(regmatches(body, gregexpr("\\\\.", body, perl = TRUE)))
  if (any(!escapes %in% c("\\'", "\\\"", "\\\\"))) {
    fail("only the escapes \\', \\\" and \\\\ are supported")
  }
  items <- .trim_space(gsub("\\\\(.)", "\\1", body, perl = TRUE))
  items[nzchar(items)]
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
