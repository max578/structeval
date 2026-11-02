# Error decomposition plot ---------------------------------------------------

#' Plot the error decomposition of each field
#'
#' @description
#' Stacked bars of correct, incorrect, missed and spurious items per field,
#' as counts or as shares of all items. Requires the ggplot2 package.
#'
#' @param counts Counts from [score()] or [simulate_counts()].
#' @param field Fields to show. `NULL` shows every field.
#' @param proportion If `TRUE`, bars show shares of each field's items;
#'   otherwise counts.
#'
#' @returns A ggplot object.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' counts <- simulate_counts(n = 100, theta = 0.9, m = 3, rho = 0.2,
#'   spu_ratio = 0.6, type = "list", seed = 1)
#' plot_errors(counts)
#' @export
plot_errors <- function(counts, field = NULL, proportion = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_errors() needs the ggplot2 package.", call. = FALSE)
  }
  fields <- .check_counts(counts, field)

  # Long table of totals ------------------------------------------------------

  outcomes <- c("cor", "inc", "mis", "spu")
  labels <- c("Correct", "Incorrect", "Missed", "Spurious")
  rows <- counts[counts$field %in% fields$field, , drop = FALSE]
  totals <- vapply(outcomes, function(o) {
    tapply(rows[[o]], factor(rows$field, levels = fields$field), sum)
  }, numeric(nrow(fields)))
  totals <- matrix(totals, nrow = nrow(fields), dimnames = list(NULL, outcomes))
  totals[is.na(totals)] <- 0
  value <- as.vector(totals)
  if (proportion) {
    value <- as.vector(totals / pmax(rowSums(totals), 1))
  }
  df <- data.frame(
    field = factor(rep(fields$field, 4L), levels = rev(fields$field)),
    outcome = factor(rep(labels, each = nrow(fields)), levels = rev(labels)),
    value = value
  )

  # Plot ----------------------------------------------------------------------

  ggplot2::ggplot(
    df, ggplot2::aes(x = .data$value, y = .data$field, fill = .data$outcome)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(
      values = c(
        Correct = "#1b9e77", Incorrect = "#7570b3", Missed = "#d95f02",
        Spurious = "#e7298a"
      ),
      breaks = labels
    ) +
    ggplot2::labs(
      x = if (proportion) "Share of items" else "Items", y = NULL,
      fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}
