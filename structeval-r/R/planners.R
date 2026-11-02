# Closed-form planners -------------------------------------------------------

#' Sample size to certify a metric above a threshold
#'
#' @description
#' The number of notes needed so that the one-sided test of H0:
#' \eqn{R \le \theta_0} has the requested power when the true value is
#' \eqn{\theta_1}.
#'
#' @details
#' The number of independent units is
#' \deqn{n_{units} = \left\lceil \frac{[z_{1-\alpha}\sqrt{\theta_0(1-\theta_0)}
#' + z_{power}\sqrt{\theta_1(1-\theta_1)}]^2}{(\theta_1-\theta_0)^2}
#' \right\rceil.}
#' The design effect is \eqn{DEFF = 1 + (m-1)\rho}. For `m = 1` the unit is
#' the note where the field is present, so `n_present = n_units`. For
#' `m > 1` the unit is the item: `n_items = ceiling(n_units * DEFF)` and
#' `n_present = ceiling(n_items / m)`. In both cases
#' `n_notes = ceiling(n_present / prevalence)`. All arguments are
#' vectorised.
#'
#' @param theta0 Threshold to certify against.
#' @param theta1 Assumed true value, greater than `theta0`.
#' @param alpha One-sided significance level.
#' @param power Target power.
#' @param m Mean items per note where the field is present.
#' @param rho Within-note correlation of item outcomes (intraclass
#'   correlation).
#' @param prevalence Proportion of notes where the field is present.
#'
#' @returns A data frame with the inputs and `deff`, `n_units`, `n_items`,
#'   `n_present` and `n_notes`.
#'
#' @seealso [power_threshold()] for the simulated power at a given n.
#'
#' @examples
#' n_threshold(0.90, 0.95)
#' n_threshold(0.90, 0.95, m = 4, rho = 0.2, prevalence = 0.5)
#' @export
n_threshold <- function(theta0, theta1, alpha = 0.025, power = 0.80, m = 1,
                        rho = 0, prevalence = 1) {
  .check_unit(theta0, "theta0")
  .check_unit(theta1, "theta1")
  .check_unit(alpha, "alpha")
  .check_unit(power, "power")
  .check_unit(prevalence, "prevalence", open = FALSE)
  if (any(theta1 <= theta0)) {
    stop("`theta1` must be greater than `theta0`.", call. = FALSE)
  }
  if (any(m < 1) || any(rho < 0 | rho >= 1) || any(prevalence <= 0)) {
    stop(
      "`m` must be at least 1, `rho` in [0, 1) and `prevalence` positive.",
      call. = FALSE
    )
  }
  num <- stats::qnorm(1 - alpha) * sqrt(theta0 * (1 - theta0)) +
    stats::qnorm(power) * sqrt(theta1 * (1 - theta1))
  n_units <- ceiling(num^2 / (theta1 - theta0)^2)
  deff <- 1 + (m - 1) * rho
  n_items <- ifelse(m > 1, ceiling(n_units * deff), n_units)
  n_present <- ifelse(m > 1, ceiling(n_items / m), n_units)
  data.frame(
    theta0 = theta0, theta1 = theta1, alpha = alpha, power = power, m = m,
    rho = rho, prevalence = prevalence, deff = deff,
    n_units = as.integer(n_units), n_items = as.integer(n_items),
    n_present = as.integer(n_present),
    n_notes = as.integer(ceiling(n_present / prevalence))
  )
}

#' Sample size for a paired comparison of note-level metrics
#'
#' @description
#' The number of notes for McNemar's test to detect a difference in a
#' note-level proportion between two systems scored on the same notes
#' (Connor 1987).
#'
#' @details
#' \deqn{n = \left\lceil \frac{[z_{1-\alpha}\sqrt{\psi} +
#' z_{power}\sqrt{\psi - \delta^2}]^2}{\delta^2} \right\rceil,}
#' where \eqn{\psi} is the proportion of discordant notes and \eqn{\delta}
#' the difference in proportions. The result is divided by `prevalence`.
#' Item-level paired designs are planned with [power_paired()].
#'
#' @param delta Difference in proportions to detect.
#' @param discordance Expected proportion of discordant notes; must exceed
#'   `delta^2`.
#' @param alpha One-sided significance level (0.025 matches a two-sided test
#'   at 0.05).
#' @param power Target power.
#' @param prevalence Proportion of notes where the field is present.
#'
#' @returns A data frame with the inputs and `n_present` and `n_notes`.
#'
#' @references
#' Connor, R. J. (1987). Sample size for testing differences in proportions
#' for the paired-sample design. *Biometrics*, 43(1), 207-211.
#'
#' @examples
#' n_paired(delta = 0.05, discordance = 0.10)
#' @export
n_paired <- function(delta, discordance, alpha = 0.025, power = 0.80,
                     prevalence = 1) {
  .check_unit(abs(delta), "delta")
  .check_unit(discordance, "discordance")
  .check_unit(alpha, "alpha")
  .check_unit(power, "power")
  .check_unit(prevalence, "prevalence", open = FALSE)
  if (any(discordance <= delta^2) || any(prevalence <= 0)) {
    stop(
      "`discordance` must exceed `delta^2` and `prevalence` must be ",
      "positive.",
      call. = FALSE
    )
  }
  num <- stats::qnorm(1 - alpha) * sqrt(discordance) +
    stats::qnorm(power) * sqrt(discordance - delta^2)
  n_present <- ceiling(num^2 / delta^2)
  data.frame(
    delta = delta, discordance = discordance, alpha = alpha, power = power,
    prevalence = prevalence, n_present = as.integer(n_present),
    n_notes = as.integer(ceiling(n_present / prevalence))
  )
}

#' Smallest count that certifies a threshold at a given n
#'
#' @description
#' The smallest number of successes \eqn{k} out of `n` whose Wilson lower
#' bound (two-sided level \eqn{1 - 2\alpha}) exceeds `theta0`.
#'
#' @param n Number of trials (note-level units), vectorised.
#' @param theta0 Threshold.
#' @param alpha One-sided significance level.
#'
#' @returns An integer vector of counts, `NA` where even `n` out of `n` does
#'   not certify.
#'
#' @examples
#' # At n = 50 only a perfect 50/50 certifies 0.90.
#' min_certifiable(50, 0.90)
#' min_certifiable(c(50, 100, 200), 0.90)
#' @export
min_certifiable <- function(n, theta0, alpha = 0.025) {
  .check_unit(theta0, "theta0")
  .check_unit(alpha, "alpha")
  if (!is.numeric(n) || any(n < 1) || any(n != round(n))) {
    stop("`n` must be positive whole numbers.", call. = FALSE)
  }
  vapply(n, function(nn) {
    k <- seq.int(0L, as.integer(nn))
    ok <- which(.wilson_ci(k, nn, 1 - 2 * alpha)[, "lower"] > theta0)
    if (length(ok) == 0L) NA_integer_ else k[min(ok)]
  }, integer(1L))
}
