# Interval engine ----------------------------------------------------------
#
# One function, ratio_ci(), computes every interval in the package. A metric
# is a ratio of note-level sums R = sum(a) / sum(b); notes with b = 0 carry no
# information about R and are dropped before any computation.

.ci_methods <- c(
  "wilson", "naive_wilson", "cluster_bca", "cluster_delta", "percentile"
)

# Evaluate expr under a fixed seed and restore the caller's RNG state.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(expr)
  }
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_seed) {
    old <- get(".Random.seed", envir = env, inherits = FALSE)
  }
  on.exit(
    if (had_seed) {
      assign(".Random.seed", old, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  )
  set.seed(seed)
  expr
}

# Wilson score interval for k successes in n trials (vectorised).
.wilson_ci <- function(k, n, level) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- k / n
  z2 <- z * z
  den <- 1 + z2 / n
  centre <- (p + z2 / (2 * n)) / den
  half <- (z / den) * sqrt(p * (1 - p) / n + z2 / (4 * n * n))
  lower <- pmax(0, centre - half)
  upper <- pmin(1, centre + half)
  empty <- rep_len(n <= 0, length(lower))
  lower[empty] <- NA_real_
  upper[empty] <- NA_real_
  cbind(lower = lower, upper = upper)
}

# One-sided score-test p-value for H0: p <= theta0 (the inversion of Wilson).
.score_p_upper <- function(k, n, theta0) {
  z <- (k / n - theta0) / sqrt(theta0 * (1 - theta0) / n)
  stats::pnorm(z, lower.tail = FALSE)
}

# Cluster-robust standard error of logit(R) by the linearised ratio variance.
.delta_se_logit <- function(a, b, r) {
  n_notes <- length(b)
  if (n_notes < 2L || r <= 0 || r >= 1) {
    return(NA_real_)
  }
  v <- n_notes / (n_notes - 1) * sum((a - r * b)^2) / sum(b)^2
  sqrt(v) / (r * (1 - r))
}

# Note-level bootstrap replicates of R via multinomial note weights.
.boot_ratio <- function(a, b, B) {
  n_notes <- length(b)
  w <- stats::rmultinom(B, n_notes, rep(1 / n_notes, n_notes))
  drop(crossprod(w, a)) / drop(crossprod(w, b))
}

# BCa percentiles (Efron 1987); NULL when the adjustment is undefined.
.bca_interval <- function(rs, est, jack, level) {
  rs <- rs[is.finite(rs)]
  if (length(rs) < 2L || stats::sd(rs) == 0 || any(!is.finite(jack))) {
    return(NULL)
  }
  z0 <- stats::qnorm(mean(rs < est) + 0.5 * mean(rs == est))
  if (!is.finite(z0)) {
    return(NULL)
  }
  d <- mean(jack) - jack
  den <- 6 * sum(d^2)^1.5
  if (!(den > 0)) {
    return(NULL)
  }
  acc <- sum(d^3) / den
  alpha <- 1 - level
  s <- z0 + stats::qnorm(c(alpha / 2, 1 - alpha / 2))
  den2 <- 1 - acc * s
  if (any(den2 <= 0)) {
    return(NULL)
  }
  q <- stats::pnorm(z0 + s / den2)
  stats::quantile(rs, q, type = 7, names = FALSE)
}

# Boundary cap for the cluster methods on n_notes notes: the lower bound never
# exceeds the perfect record's (Wilson, n_notes successes in n_notes trials)
# and the upper bound is never below the empty record's (0 successes), so a
# sample with fewer correct items never certifies more than a perfect one.
.cap_ci <- function(ci, n_notes, level) {
  w <- .wilson_ci(c(n_notes, 0), n_notes, level)
  c(min(ci[[1L]], w[1L, "lower"]), max(ci[[2L]], w[2L, "upper"]))
}

# Choose the interval method for method = "auto" (spec section 4).
.route_method <- function(b, method) {
  if (method != "auto") {
    return(method)
  }
  b <- b[b > 0]
  if (length(b) > 0L && all(b == 1)) "wilson" else "cluster_bca"
}

#' Interval for a ratio of note-level sums
#'
#' The single interval engine behind [metrics()], [test_threshold()] and the
#' power functions. The metric is \eqn{R = \sum_i a_i / \sum_i b_i} over notes
#' with \eqn{b_i > 0}.
#'
#' @param a,b Numeric vectors of note-level numerators and denominators.
#' @param method One of `"wilson"`, `"naive_wilson"`, `"cluster_bca"`,
#'   `"cluster_delta"` or `"percentile"`.
#' @param level Two-sided confidence level.
#' @param B Bootstrap replicates for the bootstrap methods.
#' @param seed Optional integer seed; the caller's RNG state is restored.
#' @return A list with `estimate`, `lower`, `upper`, `method`, `k`, `n`,
#'   `n_notes`, `boundary` and `se_logit`.
#' @noRd
ratio_ci <- function(a, b, method = "wilson", level = 0.95, B = 2000L,
                     seed = NULL) {
  method <- match.arg(method, .ci_methods)
  keep <- b > 0
  a <- as.numeric(a[keep])
  b <- as.numeric(b[keep])
  n_notes <- length(b)
  out <- list(
    estimate = NA_real_, lower = NA_real_, upper = NA_real_,
    method = method, k = sum(a), n = sum(b), n_notes = n_notes,
    boundary = FALSE, se_logit = NA_real_
  )
  if (n_notes == 0L) {
    return(out)
  }
  r <- out$k / out$n
  out$estimate <- r
  out$se_logit <- .delta_se_logit(a, b, r)

  ci <- switch(method,
    wilson = {
      if (any(b != 1)) {
        stop(
          "`wilson` needs every note-level denominator in {0, 1}; ",
          "use `cluster_bca` or `cluster_delta`.",
          call. = FALSE
        )
      }
      .wilson_ci(out$k, out$n, level)
    },
    naive_wilson = .wilson_ci(out$k, out$n, level),
    cluster_delta = {
      if (is.finite(out$se_logit) && out$se_logit > 0) {
        z <- stats::qnorm(1 - (1 - level) / 2)
        stats::plogis(stats::qlogis(r) + c(-1, 1) * z * out$se_logit)
      }
    },
    percentile = {
      rs <- .with_seed(seed, .boot_ratio(a, b, B))
      out$boundary <- r %in% c(0, 1) || stats::sd(rs) == 0
      alpha <- 1 - level
      stats::quantile(rs, c(alpha / 2, 1 - alpha / 2), names = FALSE)
    },
    cluster_bca = {
      if (!(r %in% c(0, 1)) && n_notes >= 2L) {
        rs <- .with_seed(seed, .boot_ratio(a, b, B))
        jack <- (out$k - a) / (out$n - b)
        .bca_interval(rs, r, jack, level)
      }
    }
  )

  if (is.null(ci)) {
    # Boundary rule: the clusters are the trials.
    out$boundary <- TRUE
    ci <- .wilson_ci(round(r * n_notes), n_notes, level)
  } else if (method %in% c("cluster_bca", "cluster_delta")) {
    ci <- .cap_ci(ci, n_notes, level)
  }
  out$lower <- ci[[1L]]
  out$upper <- ci[[2L]]
  out
}
