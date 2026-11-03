# Cell definitions and per-cell evaluators for the size-and-power study.
# Every interval and test reuses the structeval internals; nothing is
# re-implemented except the sharing of one bootstrap matrix between the
# percentile and BCa intervals inside a Grid A replicate.

wilson_ci <- structeval:::.wilson_ci
delta_se_logit <- structeval:::.delta_se_logit
bca_interval <- structeval:::.bca_interval
cap_ci <- structeval:::.cap_ci
sim_core <- structeval:::.sim_core
mat_ab <- structeval:::.mat_ab
t2_core <- structeval:::.t2_core
t3_core <- structeval:::.t3_core
sign_matrix <- structeval:::.sign_matrix

SEED_BASE <- 20261029L
N_LARGE <- 1e6

# Grid definitions ---------------------------------------------------------

grid_a_cells <- function() {
  prop <- expand.grid(
    n = c(20, 50, 100, 200, 400), theta = c(0.85, 0.90, 0.95, 0.98),
    m = 1, rho = 0, prevalence = 1, type = "proportion",
    stringsAsFactors = FALSE
  )
  lst <- expand.grid(
    n = c(50, 100, 200, 400), theta = c(0.85, 0.90, 0.95, 0.98),
    m = c(2, 4), rho = c(0, 0.2, 0.5), prevalence = 1, type = "list",
    stringsAsFactors = FALSE
  )
  check <- data.frame(
    n = 100, theta = 0.90, m = 4, rho = 0.2, prevalence = 0.5,
    type = "list", stringsAsFactors = FALSE
  )
  out <- rbind(prop, lst, check)
  out$grid <- "a"
  out
}

grid_b_cells <- function() {
  out <- expand.grid(
    n = c(100, 200, 400), theta = c(0.85, 0.90, 0.95), m = c(2, 4),
    rho = c(0, 0.2, 0.5), spu_ratio = c(1, 0.6, 0.4)
  )
  out$grid <- "b"
  out
}

grid_c_cells <- function() {
  out <- expand.grid(
    type = c("list", "proportion"), n = c(100, 200, 400),
    delta = c(0, 0.03, 0.05), kappa = c(0.5, 0.9),
    stringsAsFactors = FALSE
  )
  out$theta_a <- 0.90
  # Proportion cells: one item per note, rho irrelevant (key contract m1.rho0).
  out$m <- ifelse(out$type == "list", 2, 1)
  out$rho <- ifelse(out$type == "list", 0.2, 0)
  out$grid <- "c"
  out
}

planner_scenarios <- function() {
  data.frame(
    scenario = c("P1", "P2", "P3", "P4", "L1", "L2", "L3", "L4", "L5", "L6",
      "L7", "L8"),
    type = rep(c("proportion", "list"), c(4, 8)),
    theta0 = c(0.90, 0.90, 0.85, 0.85, 0.90, 0.90, 0.90, 0.90, 0.85, 0.85,
      0.85, 0.90),
    theta1 = c(0.95, 0.98, 0.90, 0.93, 0.95, 0.95, 0.95, 0.95, 0.93, 0.93,
      0.90, 0.98),
    m = c(1, 1, 1, 1, 2, 2, 4, 4, 2, 4, 4, 2),
    rho = c(0, 0, 0, 0, 0, 0.2, 0.2, 0.5, 0.5, 0, 0.2, 0.2),
    stringsAsFactors = FALSE
  )
}

# Grid A: one data set per replicate, every method on it -------------------

# Intervals of every method for one replicate. Returns a 3 x k matrix
# (lower, upper, boundary) with one column per method.
eval_intervals <- function(a, b, type, B, level = 0.95) {
  keep <- b > 0
  a <- as.numeric(a[keep])
  b <- as.numeric(b[keep])
  nn <- length(b)
  methods <- if (type == "list") {
    c("naive_wilson", "percentile", "cluster_bca", "cluster_delta")
  } else {
    c("wilson", "percentile", "cluster_bca")
  }
  out <- matrix(NA_real_, 3L, length(methods),
    dimnames = list(c("lower", "upper", "boundary"), methods)
  )
  if (nn == 0L) {
    return(out)
  }
  k <- sum(a)
  n <- sum(b)
  r <- k / n
  fallback <- wilson_ci(round(r * nn), nn, level)

  # Wilson on items (exact Wilson for proportions, naive for lists).
  out[1:2, 1] <- wilson_ci(k, n, level)
  out[3, 1] <- 0

  # One bootstrap matrix for percentile and BCa.
  w <- stats::rmultinom(B, nn, rep(1 / nn, nn))
  rs <- drop(crossprod(w, a)) / drop(crossprod(w, b))
  out[1:2, "percentile"] <- stats::quantile(
    rs, c((1 - level) / 2, 1 - (1 - level) / 2), names = FALSE
  )
  out[3, "percentile"] <- as.numeric(r %in% c(0, 1) || stats::sd(rs) == 0)

  ci <- NULL
  if (!(r %in% c(0, 1)) && nn >= 2L) {
    jack <- (k - a) / (n - b)
    ci <- bca_interval(rs, r, jack, level)
  }
  if (is.null(ci)) {
    out[1:2, "cluster_bca"] <- fallback
    out[3, "cluster_bca"] <- 1
  } else {
    out[1:2, "cluster_bca"] <- cap_ci(ci, nn, level)
    out[3, "cluster_bca"] <- 0
  }

  if (type == "list") {
    se <- delta_se_logit(a, b, r)
    if (is.finite(se) && se > 0) {
      z <- stats::qnorm(1 - (1 - level) / 2)
      out[1:2, "cluster_delta"] <- cap_ci(stats::plogis(stats::qlogis(r) +
        c(-1, 1) * z * se), nn, level)
      out[3, "cluster_delta"] <- 0
    } else {
      out[1:2, "cluster_delta"] <- fallback
      out[3, "cluster_delta"] <- 1
    }
  }
  out
}

ab_of <- function(mat, type) {
  if (type == "list") {
    list(a = 2 * mat[, "cor"], b = 2 * mat[, "cor"] + mat[, "mis"] +
      mat[, "spu"])
  } else {
    list(a = mat[, "cor"], b = mat[, "cor"] + mat[, "mis"])
  }
}

true_metric <- function(theta, m, rho, prevalence, type, spu_ratio = 1) {
  big <- sim_core(N_LARGE, theta, m, rho, prevalence, spu_ratio, type)
  ab <- ab_of(big, type)
  sum(ab$a) / sum(ab$b)
}

mc_se <- function(p, reps) sqrt(p * (1 - p) / reps)

run_cell_a <- function(cell, idx, reps, B) {
  t0 <- proc.time()[["elapsed"]]
  set.seed(SEED_BASE + idx)
  th <- true_metric(cell$theta, cell$m, cell$rho, cell$prevalence, cell$type)
  res <- vapply(seq_len(reps), function(i) {
    mat <- sim_core(cell$n, cell$theta, cell$m, cell$rho, cell$prevalence,
      1, cell$type)
    ab <- ab_of(mat, cell$type)
    as.vector(eval_intervals(ab$a, ab$b, cell$type, B))
  }, numeric(if (cell$type == "list") 12L else 9L))
  methods <- if (cell$type == "list") {
    c("naive_wilson", "percentile", "cluster_bca", "cluster_delta")
  } else {
    c("wilson", "percentile", "cluster_bca")
  }
  secs <- proc.time()[["elapsed"]] - t0
  rows <- lapply(seq_along(methods), function(j) {
    lo <- res[3 * (j - 1) + 1, ]
    up <- res[3 * (j - 1) + 2, ]
    bd <- res[3 * (j - 1) + 3, ]
    undefined <- sum(is.na(lo))
    cover <- ifelse(is.na(lo), 0, lo <= th & th <= up)
    size <- ifelse(is.na(lo), 0, lo > th)
    power <- ifelse(is.na(lo), 0, lo > th - 0.05)
    width <- up - lo
    data.frame(
      type = cell$type, n = cell$n, theta = cell$theta, m = cell$m,
      rho = cell$rho, prevalence = cell$prevalence,
      metric = if (cell$type == "list") "f1" else "recall",
      method = methods[j], reps = reps, B = B,
      coverage = mean(cover), coverage_se = mc_se(mean(cover), reps),
      width = mean(width, na.rm = TRUE),
      boundary_rate = mean(bd, na.rm = TRUE),
      size = mean(size), size_se = mc_se(mean(size), reps),
      power = mean(power), power_se = mc_se(mean(power), reps),
      theta_true_metric = th, seconds = secs,
      width_se = stats::sd(width, na.rm = TRUE) / sqrt(reps),
      n_undefined = undefined, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (cell$type == "proportion") {
    dup <- out[out$method == "wilson", ]
    dup$method <- "naive_wilson"
    out <- rbind(out, dup)
  }
  out
}

# Grid B: T3 asymmetry ------------------------------------------------------

run_cell_b <- function(cell, idx, reps, n_perm) {
  t0 <- proc.time()[["elapsed"]]
  set.seed(SEED_BASE + idx)
  signs <- sign_matrix(cell$n, n_perm)
  methods <- c("binomial", "cluster_t", "signflip")
  res <- vapply(seq_len(reps), function(i) {
    mat <- sim_core(cell$n, cell$theta, cell$m, cell$rho, 1,
      cell$spu_ratio, "list")
    vapply(methods, function(meth) {
      as.numeric(t3_core(mat[, "mis"], mat[, "spu"], meth, 0.05, n_perm,
        0.95, signs)$reject)
    }, numeric(1L))
  }, numeric(3L))
  secs <- proc.time()[["elapsed"]] - t0
  do.call(rbind, lapply(seq_along(methods), function(j) {
    x <- res[j, ]
    undefined <- sum(is.na(x))
    x[is.na(x)] <- 0
    data.frame(
      n = cell$n, theta = cell$theta, m = cell$m, rho = cell$rho,
      spu_ratio = cell$spu_ratio, method = methods[j], reps = reps,
      reject = mean(x), reject_se = mc_se(mean(x), reps), seconds = secs,
      n_undefined = undefined, stringsAsFactors = FALSE
    )
  }))
}

# Grid C: T2 paired ---------------------------------------------------------

metric_c <- function(type) if (type == "list") "f1" else "recall"

# True metric-scale difference R_B - R_A from one draw at n = 1e6.
true_delta <- function(dl, type, m, rho, kappa, seed) {
  set.seed(seed)
  s <- sim_core(N_LARGE, 0.90, m, rho, 1, 1, type, dl, kappa)
  ma <- mat_ab(s$a, metric_c(type))
  mb <- mat_ab(s$b, metric_c(type))
  sum(mb$a) / sum(mb$b) - sum(ma$a) / sum(ma$b)
}

calibrate_delta_logit <- function(target, type, m, rho, seed) {
  if (target == 0) {
    return(0)
  }
  stats::uniroot(
    function(dl) true_delta(dl, type, m, rho, 1, seed) - target,
    c(0, 3), tol = 1e-6
  )$root
}

run_cell_c <- function(cell, idx, reps, B, n_perm) {
  t0 <- proc.time()[["elapsed"]]
  methods <- c("cluster_delta", "signflip", "cluster_bca")
  if (cell$type == "proportion") methods <- c(methods, "mcnemar")
  realised <- true_delta(cell$delta_logit, cell$type, cell$m, cell$rho,
    cell$kappa, SEED_BASE + idx)
  set.seed(SEED_BASE + idx)
  signs <- sign_matrix(cell$n, n_perm)
  met <- metric_c(cell$type)
  res <- vapply(seq_len(reps), function(i) {
    s <- sim_core(cell$n, cell$theta_a, cell$m, cell$rho, 1, 1, cell$type,
      cell$delta_logit, cell$kappa)
    aa <- mat_ab(s$a, met)
    bb <- mat_ab(s$b, met)
    vapply(methods, function(meth) {
      as.numeric(t2_core(aa$a, aa$b, bb$a, bb$b, meth, 0.05, B, n_perm,
        signs)$reject)
    }, numeric(1L))
  }, numeric(length(methods)))
  res <- matrix(res, nrow = length(methods))
  secs <- proc.time()[["elapsed"]] - t0
  do.call(rbind, lapply(seq_along(methods), function(j) {
    x <- res[j, ]
    undefined <- sum(is.na(x))
    x[is.na(x)] <- 0
    data.frame(
      type = cell$type, n = cell$n, theta_a = cell$theta_a,
      delta = cell$delta, kappa = cell$kappa, m = cell$m, rho = cell$rho,
      method = methods[j], reps = reps, reject = mean(x),
      reject_se = mc_se(mean(x), reps), seconds = secs,
      n_undefined = undefined, delta_logit = cell$delta_logit,
      delta_true = realised, stringsAsFactors = FALSE
    )
  }))
}

# Planner validation --------------------------------------------------------

run_cell_p <- function(cell, idx, reps, B) {
  t0 <- proc.time()[["elapsed"]]
  out <- structeval::power_threshold(
    n = cell$n_run, theta0 = cell$theta0, theta1 = cell$theta1, m = cell$m,
    rho = cell$rho, type = cell$type, nsim = reps, B = B,
    seed = SEED_BASE + idx
  )
  data.frame(
    scenario = cell$scenario, type = cell$type, theta0 = cell$theta0,
    theta1 = cell$theta1, m = cell$m, rho = cell$rho,
    n_planned = cell$n_planned, n_run = cell$n_run, power = out$power,
    power_se = out$power_se, factor = cell$factor,
    n_undefined = out$n_undefined, method = out$method,
    seconds = proc.time()[["elapsed"]] - t0, stringsAsFactors = FALSE
  )
}
