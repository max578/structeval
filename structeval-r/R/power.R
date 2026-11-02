# Monte Carlo power ----------------------------------------------------------

#' Note-level a and b of a metric from a simulated count matrix
#' @noRd
#' @keywords internal
.mat_ab <- function(mat, metric) {
  .metric_ab(
    mat[, "cor"], 0L, mat[, "mis"], mat[, "spu"], mat[, "tn"], metric
  )
}

#' Rate and Monte Carlo standard error; an undefined test counts as no
#' rejection and is tallied
#' @noRd
#' @keywords internal
.rate <- function(x) {
  undefined <- sum(is.na(x))
  x[is.na(x)] <- 0
  r <- mean(x)
  c(rate = r, se = sqrt(r * (1 - r) / length(x)), undefined = undefined)
}

#' Simulated power of the threshold test
#'
#' @description
#' Draws `nsim` data sets from [simulate_counts()] at the true value
#' `theta1` and reports how often [test_threshold()] rejects H0:
#' \eqn{R \le \theta_0}. With `theta1 = theta0` the rate is the size of the
#' test.
#'
#' @inheritParams n_threshold
#' @param n Number of notes per data set.
#' @param theta0 Threshold under test.
#' @param theta1 True value used to generate the data.
#' @param m Mean items per present note. With `m = 1` and the default `type`
#'   the data are note-level proportions.
#' @param spu_ratio Expected spurious items per missed item.
#' @param type `"proportion"` or `"list"`; defaults to `"proportion"` when
#'   `m = 1`.
#' @param metric Metric under test; defaults to recall for proportions and
#'   F1 for lists.
#' @param method Interval method passed to the test (`"auto"` routes).
#' @param nsim Number of simulated data sets.
#' @param B Bootstrap replicates inside each data set.
#' @param seed Optional integer seed.
#'
#' @returns A one-row data frame with the design, `method` (the routed
#'   method used most often), `nsim`, `B`, `power`, `power_se` (Monte Carlo
#'   standard error), `n_undefined` (data sets where the test was undefined,
#'   counted as no rejection), `boundary_rate` and `seconds`.
#'
#' @examples
#' power_threshold(n = 100, theta0 = 0.85, theta1 = 0.95, nsim = 200,
#'   seed = 1)
#' @export
power_threshold <- function(n, theta0, theta1, m = 1, rho = 0,
                            prevalence = 1, spu_ratio = 1,
                            type = if (m == 1) "proportion" else "list",
                            metric = NULL, method = "auto", alpha = 0.025,
                            nsim = 2000L, B = 1000L, seed = NULL) {
  type <- match.arg(type, c("proportion", "list"))
  method <- match.arg(method, c("auto", .ci_methods))
  .check_sim(n, theta1, m, rho, prevalence, spu_ratio, type)
  .check_unit(theta0, "theta0")
  if (is.null(metric)) {
    metric <- if (type == "list") "f1" else "recall"
  }
  start <- proc.time()[["elapsed"]]

  # Simulation loop -----------------------------------------------------------

  res <- .with_seed(seed, vapply(seq_len(nsim), function(i) {
    mat <- .sim_core(n, theta1, m, rho, prevalence, spu_ratio, type)
    ab <- .mat_ab(mat, metric)
    row <- .t1_core(ab$a, ab$b, theta0, alpha, method, B)
    c(
      reject = as.numeric(row$reject), boundary = as.numeric(row$boundary),
      wilson = as.numeric(startsWith(row$method, "wilson"))
    )
  }, numeric(3L)))
  rate <- .rate(res["reject", ])
  routed <- if (method != "auto") {
    method
  } else if (mean(res["wilson", ]) >= 0.5) {
    "wilson"
  } else {
    "cluster_bca"
  }
  data.frame(
    n = n, theta0 = theta0, theta1 = theta1, m = m, rho = rho,
    prevalence = prevalence, type = type, metric = metric, method = routed,
    nsim = nsim, B = B, power = rate[["rate"]], power_se = rate[["se"]],
    n_undefined = as.integer(rate[["undefined"]]),
    boundary_rate = mean(res["boundary", ]),
    seconds = proc.time()[["elapsed"]] - start, stringsAsFactors = FALSE
  )
}

#' Simulated power of the asymmetry test
#'
#' @description
#' Draws `nsim` list data sets with spurious-to-missed ratio `spu_ratio` and
#' reports how often [test_asymmetry()] rejects balance. With
#' `spu_ratio = 1` the rate is the size of the test.
#'
#' @details
#' For `signflip`, one sign matrix is drawn per call and shared by every
#' simulated data set. Each p-value is valid on its own; sharing the matrix
#' makes the loop fast.
#'
#' @inheritParams power_threshold
#' @param theta Target success probability.
#' @param method `"signflip"`, `"cluster_t"` or `"binomial"`.
#' @param alpha Two-sided significance level.
#' @param n_perm Sign flips for `signflip`.
#'
#' @returns A one-row data frame with the design, `method`, `nsim`,
#'   `reject`, `reject_se`, `n_undefined` (counted as no rejection) and
#'   `seconds`.
#'
#' @examples
#' power_asymmetry(n = 100, theta = 0.9, m = 4, rho = 0.2, spu_ratio = 0.5,
#'   nsim = 100, n_perm = 1000, seed = 1)
#' @export
power_asymmetry <- function(n, theta, m = 2, rho = 0, prevalence = 1,
                            spu_ratio = 1,
                            method = c("signflip", "cluster_t", "binomial"),
                            alpha = 0.05, nsim = 2000L, n_perm = 10000L,
                            seed = NULL) {
  method <- match.arg(method)
  .check_sim(n, theta, m, rho, prevalence, spu_ratio, "list")
  start <- proc.time()[["elapsed"]]
  res <- .with_seed(seed, {
    signs <- if (method == "signflip") .sign_matrix(n, n_perm) else NULL
    vapply(seq_len(nsim), function(i) {
      mat <- .sim_core(n, theta, m, rho, prevalence, spu_ratio, "list")
      row <- .t3_core(
        mat[, "mis"], mat[, "spu"], method, alpha, n_perm, 0.95, signs
      )
      as.numeric(row$reject)
    }, numeric(1L))
  })
  rate <- .rate(res)
  data.frame(
    n = n, theta = theta, m = m, rho = rho, prevalence = prevalence,
    spu_ratio = spu_ratio, method = method, nsim = nsim,
    reject = rate[["rate"]], reject_se = rate[["se"]],
    n_undefined = as.integer(rate[["undefined"]]),
    seconds = proc.time()[["elapsed"]] - start, stringsAsFactors = FALSE
  )
}

#' Simulated power of the paired comparison
#'
#' @description
#' Draws `nsim` pairs of systems from [simulate_counts()] with a logit shift
#' `delta_logit` and item concordance `kappa`, and reports how often
#' [test_paired()] rejects \eqn{\Delta = 0}. With `delta_logit = 0` the rate
#' is the size of the test.
#'
#' @inheritParams power_threshold
#' @param theta_a Target success probability of system A.
#' @param delta_logit Logit-scale shift of system B relative to A.
#' @param kappa Item-level concordance between the systems.
#' @param method `"cluster_delta"`, `"cluster_bca"`, `"signflip"` or
#'   `"mcnemar"`.
#' @param alpha Two-sided significance level.
#' @param n_perm Sign flips for `signflip`.
#'
#' @returns A one-row data frame with the design, `metric`, `method`,
#'   `nsim`, `reject`, `reject_se`, `n_undefined` (counted as no rejection), `mean_delta`
#'   (mean observed difference)
#'   and `seconds`.
#'
#' @examples
#' power_paired(n = 100, delta_logit = 0, kappa = 0.5, m = 2, rho = 0.2,
#'   nsim = 100, seed = 1)
#' @export
power_paired <- function(n, theta_a = 0.9, delta_logit, kappa = 0.5, m = 2,
                         rho = 0.2, prevalence = 1, spu_ratio = 1,
                         type = if (m == 1) "proportion" else "list",
                         metric = NULL,
                         method = c(
                           "cluster_delta", "cluster_bca", "signflip",
                           "mcnemar"
                         ),
                         alpha = 0.05, nsim = 2000L, B = 1000L,
                         n_perm = 2000L, seed = NULL) {
  method <- match.arg(method)
  type <- match.arg(type, c("proportion", "list"))
  .check_sim(n, theta_a, m, rho, prevalence, spu_ratio, type)
  .check_unit(kappa, "kappa", open = FALSE)
  if (is.null(metric)) {
    metric <- if (type == "list") "f1" else "recall"
  }
  start <- proc.time()[["elapsed"]]
  res <- .with_seed(seed, {
    signs <- if (method == "signflip") .sign_matrix(n, n_perm) else NULL
    vapply(seq_len(nsim), function(i) {
      sim <- .sim_core(
        n, theta_a, m, rho, prevalence, spu_ratio, type, delta_logit, kappa
      )
      ab_a <- .mat_ab(sim$a, metric)
      ab_b <- .mat_ab(sim$b, metric)
      row <- .t2_core(
        ab_a$a, ab_a$b, ab_b$a, ab_b$b, method, alpha, B, n_perm, signs
      )
      c(as.numeric(row$reject), row$estimate)
    }, numeric(2L))
  })
  rate <- .rate(res[1L, ])
  data.frame(
    n = n, theta_a = theta_a, delta_logit = delta_logit, kappa = kappa,
    m = m, rho = rho, prevalence = prevalence, type = type, metric = metric,
    method = method, nsim = nsim, reject = rate[["rate"]],
    reject_se = rate[["se"]], n_undefined = as.integer(rate[["undefined"]]),
    mean_delta = mean(res[2L, ], na.rm = TRUE),
    seconds = proc.time()[["elapsed"]] - start, stringsAsFactors = FALSE
  )
}
