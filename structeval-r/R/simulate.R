# Generative model of clustered extraction errors ----------------------------

#' Draw note-level counts as a matrix (internal, used in the power loops)
#' @noRd
#' @keywords internal
.sim_core <- function(n, theta, m, rho, prevalence, spu_ratio, type,
                      delta_logit = NULL, kappa = 0) {
  present <- stats::runif(n) < prevalence
  n_present <- sum(present)

  # Note-level success probabilities ------------------------------------------

  p <- if (rho == 0) {
    rep(theta, n_present)
  } else {
    stats::rbeta(
      n_present, theta * (1 - rho) / rho, (1 - theta) * (1 - rho) / rho
    )
  }
  m_i <- if (type == "list") {
    1L + stats::rpois(n_present, m - 1)
  } else {
    rep(1L, n_present)
  }
  spread <- function(cor, mis, spu) {
    out <- matrix(0L, n, 4L, dimnames = list(NULL, c("cor", "mis", "spu", "tn")))
    out[present, "cor"] <- cor
    out[present, "mis"] <- mis
    out[present, "spu"] <- spu
    out[!present, "tn"] <- 1L
    out
  }

  # One system ----------------------------------------------------------------

  if (is.null(delta_logit)) {
    cor <- stats::rbinom(n_present, m_i, p)
    spu <- if (type == "list") {
      stats::rpois(n_present, spu_ratio * m_i * (1 - p))
    } else {
      integer(n_present)
    }
    return(spread(cor, m_i - cor, spu))
  }

  # Two systems with item-level concordance kappa -----------------------------

  p_b <- stats::plogis(stats::qlogis(p) + delta_logit)
  idx <- rep.int(seq_len(n_present), m_i)
  u_a <- stats::runif(length(idx))
  shared <- stats::runif(length(idx)) < kappa
  u_b <- ifelse(shared, u_a, stats::runif(length(idx)))
  cor_a <- tabulate(idx[u_a < p[idx]], n_present)
  cor_b <- tabulate(idx[u_b < p_b[idx]], n_present)
  spu_a <- spu_b <- integer(n_present)
  if (type == "list") {
    spu_a <- stats::rpois(n_present, spu_ratio * m_i * (1 - p))
    spu_b <- stats::rpois(n_present, spu_ratio * m_i * (1 - p_b))
  }
  list(
    a = spread(cor_a, m_i - cor_a, spu_a),
    b = spread(cor_b, m_i - cor_b, spu_b)
  )
}

#' Matrix of counts to the counts data frame
#' @noRd
#' @keywords internal
.as_counts <- function(mat, field, type) {
  out <- data.frame(
    note = as.character(seq_len(nrow(mat))), field = field,
    type = type, cor = mat[, "cor"], inc = 0L, mis = mat[, "mis"],
    spu = mat[, "spu"], tn = mat[, "tn"], stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

#' Check generative-model arguments
#' @noRd
#' @keywords internal
.check_sim <- function(n, theta, m, rho, prevalence, spu_ratio, type) {
  .check_unit(theta, "theta")
  .check_unit(prevalence, "prevalence", open = FALSE)
  if (length(n) != 1L || n < 2 || n != round(n)) {
    stop("`n` must be a whole number of at least 2.", call. = FALSE)
  }
  if (rho < 0 || rho >= 1) {
    stop("`rho` must be in [0, 1).", call. = FALSE)
  }
  if (m < 1 || spu_ratio < 0) {
    stop("`m` must be at least 1 and `spu_ratio` non-negative.", call. = FALSE)
  }
  if (type == "proportion" && m != 1) {
    stop("`type = \"proportion\"` has one item per note; set `m = 1`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Simulate note-level counts from the generative model
#'
#' @description
#' Draws counts for `n` notes from a model of clustered extraction errors
#' with a target accuracy `theta`, within-note correlation `rho` and an
#' expected ratio of spurious to missed items `spu_ratio`.
#'
#' @details
#' A note holds the field with probability `prevalence`; a note without it
#' contributes `tn = 1` and nothing else. A note with the field has success
#' probability \eqn{p_i = \theta} when \eqn{\rho = 0}, otherwise
#' \eqn{p_i \sim Beta(\theta(1-\rho)/\rho, (1-\theta)(1-\rho)/\rho)}, whose
#' intraclass correlation is \eqn{\rho}.
#'
#' * `"proportion"`: one trial per note, `cor` with probability \eqn{p_i},
#'   else `mis`. Returned with type `"scalar"`.
#' * `"list"`: \eqn{m_i = 1 + Poisson(m - 1)} reference items,
#'   \eqn{cor_i \sim Binomial(m_i, p_i)}, \eqn{mis_i = m_i - cor_i} and
#'   \eqn{spu_i \sim Poisson(r\, m_i (1 - p_i))} where \eqn{r} is
#'   `spu_ratio`. With \eqn{r = 1}, precision, recall and F1 all equal \eqn{\theta} as ratios of
#'   expectations.
#'
#' With `delta_logit` set, two systems are drawn on the same notes. System A
#' uses \eqn{p_i}; system B uses \eqn{expit(logit(p_i) + \delta)}. Each
#' reference item draws a uniform \eqn{U_{ij}} for A; B reuses it with
#' probability `kappa` and draws its own otherwise. The difference in a
#' metric implied by `delta_logit` can be computed by one draw at large `n`.
#'
#' @param n Number of notes.
#' @param theta Target success probability (F1 = precision = recall when
#'   `spu_ratio = 1`).
#' @param m Mean reference items per note where the field is present.
#' @param rho Within-note correlation, in \[0, 1).
#' @param prevalence Probability that a note holds the field.
#' @param spu_ratio Expected spurious items per missed item.
#' @param type `"proportion"` or `"list"`.
#' @param delta_logit Optional logit-scale shift of system B relative to A.
#' @param kappa Item-level concordance between the two systems, in \[0, 1].
#' @param seed Optional integer seed. The caller's random number state is
#'   restored afterwards.
#' @param field Field name written into the output.
#'
#' @returns A counts data frame in the shape returned by [score()], or with
#'   `delta_logit` a list with elements `a` and `b`.
#'
#' @examples
#' counts <- simulate_counts(n = 100, theta = 0.9, m = 4, rho = 0.2,
#'   type = "list", seed = 1)
#' head(counts)
#' metrics(counts, metric = "f1", seed = 1)
#' @export
simulate_counts <- function(n, theta, m = 1, rho = 0, prevalence = 1,
                            spu_ratio = 1, type = c("proportion", "list"),
                            delta_logit = NULL, kappa = 0, seed = NULL,
                            field = "field") {
  type <- match.arg(type)
  .check_sim(n, theta, m, rho, prevalence, spu_ratio, type)
  .check_unit(kappa, "kappa", open = FALSE)
  out_type <- if (type == "list") "list" else "scalar"
  sim <- .with_seed(seed, .sim_core(
    n, theta, m, rho, prevalence, spu_ratio, type, delta_logit, kappa
  ))
  if (is.null(delta_logit)) {
    return(.as_counts(sim, field, out_type))
  }
  list(
    a = .as_counts(sim$a, field, out_type),
    b = .as_counts(sim$b, field, out_type)
  )
}
