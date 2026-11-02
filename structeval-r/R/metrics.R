# Metrics with routed intervals --------------------------------------------

.metric_names <- c(
  "precision", "recall", "f1", "specificity", "accuracy", "mis_share"
)

#' Note-level numerators and denominators of a metric
#' @noRd
#' @keywords internal
.metric_ab <- function(cor, inc, mis, spu, tn, metric) {
  switch(metric,
    precision = list(a = cor, b = cor + spu + inc),
    recall = list(a = cor, b = cor + mis + inc),
    f1 = list(a = 2 * cor, b = 2 * cor + mis + spu + 2 * inc),
    specificity = list(a = tn, b = tn + as.numeric(cor + mis + inc == 0 & spu > 0)),
    accuracy = list(a = cor + tn, b = cor + tn + inc + mis + spu),
    mis_share = list(a = mis, b = mis + spu),
    stop("Unknown metric '", metric, "'.", call. = FALSE)
  )
}

#' Default metrics for a field type
#' @noRd
#' @keywords internal
.default_metrics <- function(type) {
  if (type == "list") {
    c("precision", "recall", "f1", "specificity")
  } else {
    c("precision", "recall", "f1", "specificity", "accuracy")
  }
}

#' Field table of a counts frame (keeps fields with no scored rows)
#' @noRd
#' @keywords internal
.field_table <- function(counts) {
  fields <- attr(counts, "fields")
  if (is.null(fields)) {
    fields <- unique(counts[, c("field", "type")])
  }
  rownames(fields) <- NULL
  fields
}

#' Validate a counts frame and select fields
#' @noRd
#' @keywords internal
.check_counts <- function(counts, field = NULL) {
  cols <- c("note", "field", "type", "cor", "inc", "mis", "spu", "tn")
  absent <- setdiff(cols, names(counts))
  if (!is.data.frame(counts) || length(absent) > 0L) {
    stop(
      "`counts` must be a data frame from score() or simulate_counts() with ",
      "columns ", paste0("'", cols, "'", collapse = ", "), ".",
      call. = FALSE
    )
  }
  fields <- .field_table(counts)
  if (!is.null(field)) {
    unknown <- setdiff(field, fields$field)
    if (length(unknown) > 0L) {
      stop(
        "Field(s) not in `counts`: ",
        paste0("'", unknown, "'", collapse = ", "), ".",
        call. = FALSE
      )
    }
    fields <- fields[match(field, fields$field), , drop = FALSE]
  }
  fields
}

#' Numerators and denominators of a metric for one field
#' @noRd
#' @keywords internal
.field_ab <- function(counts, field, metric) {
  rows <- counts[counts$field == field, , drop = FALSE]
  ab <- .metric_ab(rows$cor, rows$inc, rows$mis, rows$spu, rows$tn, metric)
  ab$note <- rows$note
  ab
}

#' Point estimates and intervals for each field
#'
#' @description
#' Pools counts over notes (micro averaging) and reports each metric as a
#' ratio of note-level sums, with an interval matched to the structure of
#' the field.
#'
#' @details
#' Each metric is \eqn{R = \sum_i a_i / \sum_i b_i} over notes:
#'
#' | metric | \eqn{a_i} | \eqn{b_i} |
#' |---|---|---|
#' | precision | cor | cor + spu + inc |
#' | recall | cor | cor + mis + inc |
#' | f1 | 2 cor | 2 cor + mis + spu + 2 inc |
#' | specificity | tn | tn + (spurious note with an empty reference) |
#' | accuracy | cor + tn | cor + tn + inc + mis + spu |
#' | mis_share | mis | mis + spu |
#'
#' F1 computed this way equals \eqn{2PR/(P+R)}. Accuracy applies to binary
#' and scalar fields. Specificity is a note-level proportion for every type.
#'
#' With `method = "auto"`, a metric whose note-level denominators are all 0
#' or 1 (binary and scalar precision, recall, specificity and accuracy) gets
#' the Wilson score interval. Every other metric (list precision, recall and
#' F1, scalar F1) gets the note-level cluster bootstrap with BCa percentiles.
#'
#' **Boundary rule.** For `cluster_bca` and `cluster_delta`, when the
#' estimate is 0 or 1, the bootstrap distribution has no spread, or the
#' acceleration is undefined, the interval is the Wilson interval on
#' `round(R * n_notes)` successes in `n_notes` trials (the notes are the
#' trials) and `boundary` is `TRUE`. No zero-width interval is returned. For
#' `percentile`, `boundary` flags the same condition but the percentile
#' interval is returned unchanged, so its failure stays visible.
#'
#' A metric with \eqn{\sum b_i = 0} is undefined: the estimate and interval
#' are `NA`, `n` is 0, and a warning names the field and the metric.
#'
#' @param counts Counts from [score()] or [simulate_counts()].
#' @param field Fields to report. `NULL` reports every field.
#' @param metric Metrics to report. `NULL` reports the default set for each
#'   field type: precision, recall, F1 and specificity, plus accuracy for
#'   binary and scalar fields. `"mis_share"` is also available.
#' @param method Interval method: `"auto"` (the router described in
#'   Details), `"wilson"`, `"naive_wilson"` (Wilson on the pooled item counts,
#'   which treats items as independent), `"cluster_bca"`, `"cluster_delta"`
#'   (linearised ratio variance on the logit scale) or `"percentile"` (note
#'   level percentile bootstrap).
#' @param level Two-sided confidence level.
#' @param B Bootstrap replicates.
#' @param seed Optional integer seed for the bootstrap. The caller's random
#'   number state is restored afterwards.
#'
#' @returns A data frame with one row per field and metric and columns
#'   `field`, `type`, `metric`, `estimate`, `lower`, `upper`, `method`, `k`
#'   (\eqn{\sum a_i}), `n` (\eqn{\sum b_i}), `n_notes` (notes with
#'   \eqn{b_i > 0}) and `boundary`.
#'
#' @references
#' Wilson, E. B. (1927). Probable inference, the law of succession, and
#' statistical inference. *Journal of the American Statistical Association*,
#' 22(158), 209-212.
#'
#' Efron, B. (1987). Better bootstrap confidence intervals. *Journal of the
#' American Statistical Association*, 82(397), 171-185.
#'
#' @seealso [test_threshold()] to certify a metric against a threshold.
#'
#' @examples
#' counts <- simulate_counts(n = 100, theta = 0.9, m = 3, rho = 0.2,
#'   type = "list", seed = 1)
#' metrics(counts, seed = 1)
#' metrics(counts, metric = "f1", method = "naive_wilson")
#' @export
metrics <- function(counts, field = NULL, metric = NULL, method = "auto",
                    level = 0.95, B = 2000L, seed = NULL) {
  method <- match.arg(method, c("auto", .ci_methods))
  fields <- .check_counts(counts, field)
  if (!is.null(metric)) {
    unknown <- setdiff(metric, .metric_names)
    if (length(unknown) > 0L) {
      stop(
        "Unknown metric(s): ", paste0("'", unknown, "'", collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  # One row per field and metric ----------------------------------------------

  rows <- .with_seed(seed, {
    out <- list()
    for (i in seq_len(nrow(fields))) {
      f <- fields$field[i]
      type <- fields$type[i]
      wanted <- if (is.null(metric)) .default_metrics(type) else metric
      if (type == "list" && "accuracy" %in% wanted) {
        stop(
          "Accuracy is not defined for list field '", f, "'.",
          call. = FALSE
        )
      }
      undefined <- character(0)
      for (m in wanted) {
        ab <- .field_ab(counts, f, m)
        used <- .route_method(ab$b, method)
        res <- ratio_ci(ab$a, ab$b, used, level = level, B = B)
        if (res$n_notes == 0L) {
          undefined <- c(undefined, m)
        }
        out[[length(out) + 1L]] <- data.frame(
          field = f, type = type, metric = m, estimate = res$estimate,
          lower = res$lower, upper = res$upper,
          method = if (res$n_notes == 0L) "none" else used,
          k = res$k, n = res$n, n_notes = res$n_notes,
          boundary = res$boundary, stringsAsFactors = FALSE
        )
      }
      if (length(undefined) > 0L) {
        warning(
          "Metric(s) ", paste0("'", undefined, "'", collapse = ", "),
          " undefined for field '", f, "' (zero denominator); reported as NA.",
          call. = FALSE
        )
      }
    }
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
