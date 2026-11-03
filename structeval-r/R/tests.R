# Hypothesis tests T1 to T3 ------------------------------------------------

#' Test row in the common output shape
#' @noRd
#' @keywords internal
.test_row <- function(estimate = NA_real_, lower = NA_real_, upper = NA_real_,
                      p_value = NA_real_, reject = NA, method = NA_character_,
                      n = 0, n_notes = 0L, boundary = FALSE) {
  data.frame(
    estimate = estimate, lower = lower, upper = upper, p_value = p_value,
    reject = reject, method = method, n = n, n_notes = n_notes,
    boundary = boundary, stringsAsFactors = FALSE
  )
}

#' Two-sided sign-flip p-value for the sum of x
#' @noRd
#' @keywords internal
.signflip_p <- function(x, n_perm, signs = NULL) {
  if (is.null(signs)) {
    signs <- .sign_matrix(length(x), n_perm)
  } else {
    signs <- signs[seq_along(x), , drop = FALSE]
  }
  stat <- abs(sum(x))
  perm <- abs(drop(crossprod(signs, x)))
  tol <- sqrt(.Machine$double.eps) * max(1, stat)
  (1 + sum(perm >= stat - tol)) / (1 + ncol(signs))
}

#' Draw a sign matrix for the sign-flip tests
#' @noRd
#' @keywords internal
.sign_matrix <- function(n, n_perm) {
  matrix(2L * stats::rbinom(n * n_perm, 1L, 0.5) - 1L, nrow = n)
}

# T1: threshold test -------------------------------------------------------

#' Threshold test on note-level vectors
#' @noRd
#' @keywords internal
.t1_core <- function(a, b, theta0, alpha, method, B) {
  used <- .route_method(b, method)
  res <- ratio_ci(a, b, used, level = 1 - 2 * alpha, B = B)
  if (res$n_notes == 0L) {
    return(.test_row(method = "none"))
  }
  r <- res$estimate
  if (used %in% c("wilson", "naive_wilson")) {
    p <- .score_p_upper(res$k, res$n, theta0)
    label <- used
  } else if (is.finite(res$se_logit) && res$se_logit > 0) {
    z <- (stats::qlogis(r) - stats::qlogis(theta0)) / res$se_logit
    p <- stats::pnorm(z, lower.tail = FALSE)
    label <- if (used == "cluster_delta") {
      used
    } else {
      paste0(used, " (p: cluster_delta)")
    }
    if (used %in% c("cluster_bca", "cluster_delta")) {
      # The boundary cap: never below the perfect record's p on these notes.
      p_cap <- .score_p_upper(res$n_notes, res$n_notes, theta0)
      if (p_cap > p) {
        p <- p_cap
        label <- paste0(used, " (p: boundary cap)")
      }
    }
  } else {
    p <- .score_p_upper(round(r * res$n_notes), res$n_notes, theta0)
    label <- paste0(used, " (p: wilson on notes)")
  }
  .test_row(
    estimate = r, lower = res$lower, upper = res$upper, p_value = p,
    reject = res$lower > theta0, method = label, n = res$n,
    n_notes = res$n_notes, boundary = res$boundary
  )
}

#' Test whether a metric exceeds a threshold (T1)
#'
#' @description
#' Tests H0: \eqn{R \le \theta_0} against \eqn{R > \theta_0} for each field.
#' The decision is that the lower bound of the two-sided
#' \eqn{1 - 2\alpha} interval from the routed method exceeds
#' \eqn{\theta_0}.
#'
#' @details
#' For the Wilson methods the p-value is the one-sided score test
#' \eqn{z = (R - \theta_0)/\sqrt{\theta_0(1-\theta_0)/n}}, which gives the
#' same decision as the Wilson bound. For the cluster methods it is the
#' cluster-robust logit test
#' \eqn{z = (\mathrm{logit} R - \mathrm{logit}\,\theta_0)/SE_{logit}}; for
#' `cluster_bca` and `percentile` this p-value comes from the
#' `cluster_delta` variance and the `method` column says so. At the
#' boundary the score test on notes as trials is used. Under `cluster_bca`
#' or `cluster_delta` the p-value is never below that score test for a
#' perfect record on the same notes, in line with the boundary cap of
#' [metrics()], and the `method` column then reads `(p: boundary cap)`.
#'
#' @inheritParams metrics
#' @param field One or more fields.
#' @param metric Metric to test, recycled over `field`.
#' @param theta0 Threshold, recycled over `field`.
#' @param alpha One-sided significance level. The default 0.025 matches the
#'   lower bound of a two-sided 95% interval.
#' @param adjust `"none"` or `"holm"`. With `"holm"` a column `p_adjusted`
#'   holds Holm-adjusted p-values across the rows.
#'
#' @returns A data frame with one row per field and columns `field`,
#'   `metric`, `theta0`, `estimate`, `lower`, `upper`, `p_value`, `reject`,
#'   `method`, `n`, `n_notes` and `boundary`, plus `p_adjusted` when
#'   requested.
#'
#' @seealso [n_threshold()] and [power_threshold()] to plan the sample.
#'
#' @examples
#' counts <- simulate_counts(n = 200, theta = 0.95, type = "proportion",
#'   seed = 2)
#' test_threshold(counts, "field", "recall", theta0 = 0.90)
#' @export
test_threshold <- function(counts, field, metric, theta0, alpha = 0.025,
                           method = "auto", B = 2000L, seed = NULL,
                           adjust = c("none", "holm")) {
  method <- match.arg(method, c("auto", .ci_methods))
  adjust <- match.arg(adjust)
  .check_counts(counts, field)
  .check_unit(theta0, "theta0", open = TRUE)
  metric <- rep_len(metric, length(field))
  theta0 <- rep_len(theta0, length(field))

  rows <- .with_seed(seed, lapply(seq_along(field), function(i) {
    ab <- .field_ab(counts, field[i], metric[i])
    row <- .t1_core(ab$a, ab$b, theta0[i], alpha, method, B)
    if (row$n_notes == 0L) {
      warning(
        "Metric '", metric[i], "' is undefined for field '", field[i],
        "' (zero denominator); reported as NA.",
        call. = FALSE
      )
    }
    cbind(
      data.frame(
        field = field[i], metric = metric[i], theta0 = theta0[i],
        stringsAsFactors = FALSE
      ),
      row
    )
  }))
  out <- do.call(rbind, rows)
  if (adjust == "holm") {
    out$p_adjusted <- stats::p.adjust(out$p_value, method = "holm")
  }
  rownames(out) <- NULL
  out
}

# T2: paired comparison ------------------------------------------------------

#' Paired comparison on matched note-level vectors
#' @noRd
#' @keywords internal
.t2_core <- function(aa, ba, ab, bb, method, alpha, B, n_perm, signs = NULL) {
  keep <- ba + bb > 0
  aa <- aa[keep]
  ba <- ba[keep]
  ab <- ab[keep]
  bb <- bb[keep]
  n_notes <- length(aa)
  if (n_notes < 2L || sum(ba) == 0 || sum(bb) == 0) {
    return(.test_row(method = method, n_notes = n_notes))
  }
  ra <- sum(aa) / sum(ba)
  rb <- sum(ab) / sum(bb)
  delta <- ra - rb

  # Cluster-robust variance of the difference ---------------------------------

  u <- (aa - ra * ba) / sum(ba) - (ab - rb * bb) / sum(bb)
  se <- sqrt(n_notes / (n_notes - 1) * sum((u - mean(u))^2))
  z_crit <- stats::qnorm(1 - alpha / 2)
  degenerate <- !(se > 0)
  lower <- delta - z_crit * se
  upper <- delta + z_crit * se
  p_delta <- if (degenerate) NA_real_ else 2 * stats::pnorm(-abs(delta / se))

  # Method-specific p-value and interval --------------------------------------

  p <- p_delta
  label <- method
  boundary <- degenerate
  if (method == "cluster_bca" && !degenerate) {
    w <- stats::rmultinom(B, n_notes, rep(1 / n_notes, n_notes))
    rs <- drop(crossprod(w, aa)) / drop(crossprod(w, ba)) -
      drop(crossprod(w, ab)) / drop(crossprod(w, bb))
    jack <- (sum(aa) - aa) / (sum(ba) - ba) - (sum(ab) - ab) / (sum(bb) - bb)
    ci <- .bca_interval(rs, delta, jack, 1 - alpha)
    if (is.null(ci)) {
      boundary <- TRUE
      label <- "cluster_bca (interval: cluster_delta)"
    } else {
      lower <- ci[1L]
      upper <- ci[2L]
      label <- "cluster_bca (p: cluster_delta)"
    }
  } else if (method == "signflip") {
    # Influence values centred at the pooled ratio sum exactly to delta.
    r0 <- (sum(aa) + sum(ab)) / (sum(ba) + sum(bb))
    v <- (aa - r0 * ba) / sum(ba) - (ab - r0 * bb) / sum(bb)
    p <- .signflip_p(v, n_perm, signs)
    label <- "signflip (interval: cluster_delta)"
  } else if (method == "mcnemar") {
    if (any(!ba %in% c(0, 1)) || any(ba != bb)) {
      stop(
        "`mcnemar` needs the same note-level 0/1 denominator for both ",
        "systems (binary or scalar recall, accuracy, or specificity).",
        call. = FALSE
      )
    }
    n10 <- sum(aa == 1 & ab == 0)
    n01 <- sum(aa == 0 & ab == 1)
    p <- if (n10 + n01 == 0L) 1 else stats::binom.test(n10, n10 + n01)$p.value
    label <- "mcnemar (interval: cluster_delta)"
  }
  reject <- if (method == "cluster_bca" && !boundary) {
    lower > 0 || upper < 0
  } else {
    p < alpha
  }
  out <- .test_row(
    estimate = delta, lower = lower, upper = upper, p_value = p,
    reject = reject, method = label, n = sum(ba), n_notes = n_notes,
    boundary = boundary
  )
  out$estimate_a <- ra
  out$estimate_b <- rb
  out$n_b <- sum(bb)
  out
}

#' Compare two systems scored on the same notes (T2)
#'
#' @description
#' Tests H0: \eqn{\Delta = R_A - R_B = 0} for one metric, with the notes as
#' matched pairs.
#'
#' @details
#' Notes are matched by `note`. With influence values
#' \eqn{u_i = (a^A_i - R_A b^A_i)/\sum b^A - (a^B_i - R_B b^B_i)/\sum b^B},
#' the cluster-robust variance is
#' \eqn{[n/(n-1)] \sum (u_i - \bar u)^2}.
#'
#' * `cluster_delta`: Wald test and interval from that variance.
#' * `cluster_bca`: paired note resampling with BCa percentiles; the
#'   decision is whether the interval excludes 0 and the p-value is the
#'   `cluster_delta` one.
#' * `signflip`: random sign flips of the influence values centred at the
#'   pooled ratio \eqn{R_0}, whose sum equals \eqn{\Delta} exactly
#'   (two-sided). The uncentred \eqn{u_i} sum to zero and carry no test.
#' * `mcnemar`: exact binomial test on discordant notes; only for note-level
#'   0/1 metrics with the same denominator in both systems.
#'
#' @inheritParams metrics
#' @param counts_a,counts_b Counts for system A and system B.
#' @param field One field.
#' @param metric One metric.
#' @param method Test method, see Details.
#' @param alpha Two-sided significance level.
#' @param n_perm Sign flips for `signflip`.
#'
#' @returns A one-row data frame with `field`, `metric`, `estimate`
#'   (\eqn{\Delta}), `lower`, `upper`, `p_value`, `reject`, `method`, `n`
#'   (\eqn{\sum b^A}), `n_notes`, `boundary`, `estimate_a`, `estimate_b` and
#'   `n_b` (\eqn{\sum b^B}).
#'
#' @seealso [n_paired()] and [power_paired()].
#'
#' @examples
#' pair <- simulate_counts(n = 200, theta = 0.9, m = 2, rho = 0.2,
#'   type = "list", delta_logit = 0.5, kappa = 0.5, seed = 3)
#' test_paired(pair$a, pair$b, "field", "f1")
#' @export
test_paired <- function(counts_a, counts_b, field, metric,
                        method = c(
                          "cluster_delta", "cluster_bca", "signflip",
                          "mcnemar"
                        ),
                        alpha = 0.05, B = 2000L, n_perm = 10000L,
                        seed = NULL) {
  method <- match.arg(method)
  .check_counts(counts_a, field)
  .check_counts(counts_b, field)

  # Match notes ---------------------------------------------------------------

  ab_a <- .field_ab(counts_a, field, metric)
  ab_b <- .field_ab(counts_b, field, metric)
  notes <- intersect(ab_a$note, ab_b$note)
  unmatched <- length(union(ab_a$note, ab_b$note)) - length(notes)
  if (unmatched > 0L) {
    warning(
      unmatched, " note(s) scored for only one system were dropped.",
      call. = FALSE
    )
  }
  ia <- match(notes, ab_a$note)
  ib <- match(notes, ab_b$note)

  res <- .with_seed(seed, .t2_core(
    ab_a$a[ia], ab_a$b[ia], ab_b$a[ib], ab_b$b[ib],
    method, alpha, B, n_perm
  ))
  if (is.na(res$p_value)) {
    warning(
      "The paired test is undefined for field '", field, "' (too few ",
      "notes or zero note-level variance of the difference).",
      call. = FALSE
    )
  }
  out <- cbind(
    data.frame(field = field, metric = metric, stringsAsFactors = FALSE),
    res
  )
  rownames(out) <- NULL
  out
}

# T3: omission versus commission ---------------------------------------------

#' Asymmetry test on note-level mis and spu
#' @noRd
#' @keywords internal
.t3_core <- function(mis, spu, method, alpha, n_perm, level, signs = NULL) {
  d <- mis - spu
  share <- ratio_ci(mis, mis + spu, "cluster_delta", level = level)
  total <- sum(mis) + sum(spu)
  if (total == 0) {
    return(.test_row(method = method, n_notes = share$n_notes))
  }
  p <- switch(method,
    signflip = .signflip_p(d, n_perm, signs),
    cluster_t = {
      if (length(d) < 2L) {
        NA_real_
      } else if (stats::sd(d) > 0) {
        stats::t.test(d)$p.value
      } else if (all(d == 0)) {
        1
      } else {
        NA_real_
      }
    },
    binomial = stats::binom.test(sum(mis), total)$p.value
  )
  .test_row(
    estimate = share$estimate, lower = share$lower, upper = share$upper,
    p_value = p, reject = p < alpha, method = method, n = total,
    n_notes = share$n_notes, boundary = share$boundary
  )
}

#' Test whether omissions and commissions are balanced (T3)
#'
#' @description
#' Tests H0: \eqn{E[d_i] = 0} for the per-note difference
#' \eqn{d_i = mis_i - spu_i}, that is, whether the system misses reference
#' items as often as it adds spurious ones.
#'
#' @details
#' * `signflip`: random sign flips of \eqn{d_i} over every scored note of
#'   the field, two-sided. Exact under within-note symmetry of the two error
#'   labels.
#' * `cluster_t`: one-sample t test on \eqn{d_i}; `NA` with a warning when
#'   the field has fewer than two notes or constant non-zero differences.
#' * `binomial`: exact binomial test of \eqn{\sum mis} out of
#'   \eqn{\sum mis + \sum spu} against 0.5. Items are treated as independent;
#'   this is the naive comparison method.
#'
#' A list substitution (a wrong item in place of a right one) adds one `mis`
#' and one `spu` to the same note. It leaves \eqn{d_i}, and so the sign-flip
#' test, unchanged, while it dilutes the binomial share toward 0.5.
#'
#' The `estimate` is the MIS share \eqn{\sum mis / (\sum mis + \sum spu)}
#' with its `cluster_delta` interval at level `level`.
#'
#' @inheritParams test_paired
#' @param counts Counts from [score()] or [simulate_counts()].
#' @param level Confidence level for the MIS share interval.
#'
#' @returns A one-row data frame with `field`, `estimate`, `lower`, `upper`,
#'   `p_value`, `reject`, `method`, `n` (\eqn{\sum mis + \sum spu}),
#'   `n_notes` and `boundary`.
#'
#' @seealso [power_asymmetry()].
#'
#' @examples
#' counts <- simulate_counts(n = 150, theta = 0.9, m = 4, rho = 0.2,
#'   spu_ratio = 0.5, type = "list", seed = 4)
#' test_asymmetry(counts, "field", seed = 4)
#' test_asymmetry(counts, "field", method = "binomial")
#' @export
test_asymmetry <- function(counts, field,
                           method = c("signflip", "cluster_t", "binomial"),
                           alpha = 0.05, n_perm = 10000L, level = 0.95,
                           seed = NULL) {
  method <- match.arg(method)
  .check_counts(counts, field)
  rows <- counts[counts$field == field, , drop = FALSE]
  res <- .with_seed(seed, .t3_core(
    rows$mis, rows$spu, method, alpha, n_perm, level
  ))
  if (res$n == 0) {
    warning(
      "No missed or spurious items for field '", field,
      "'; the asymmetry test is undefined.",
      call. = FALSE
    )
  } else if (is.na(res$p_value) && nrow(rows) < 2L) {
    warning(
      "Field '", field, "' has fewer than two notes; the t test is ",
      "undefined.",
      call. = FALSE
    )
  } else if (is.na(res$p_value)) {
    warning(
      "The per-note differences are constant and non-zero for field '",
      field, "'; the t test is undefined.",
      call. = FALSE
    )
  }
  out <- cbind(
    data.frame(field = field, stringsAsFactors = FALSE),
    res
  )
  rownames(out) <- NULL
  out
}

#' Check a probability argument
#' @noRd
#' @keywords internal
.check_unit <- function(x, name, open = TRUE) {
  ok <- is.numeric(x) && length(x) > 0L && all(is.finite(x)) &&
    if (open) all(x > 0 & x < 1) else all(x >= 0 & x <= 1)
  if (!ok) {
    stop(
      "`", name, "` must be a number in ", if (open) "(0, 1)" else "[0, 1]",
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
