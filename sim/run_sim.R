# =============================================================================
# Size-and-power simulation study for the structeval methods paper
#
# Author      : Max Moldovan
# Version date: 29 October 2026
# Usage       : Rscript sim/run_sim.R [reps]   (run from the project root;
#               reps defaults to 2000; finished cells are cached in
#               sim/results/cells/ and skipped on a re-run)
#
# Data        : simulated from structeval::simulate_counts() (spec section 7);
#               deployment counts from the base paper Table 1 (spec section 9)
# Requirements: structeval 0.1.0, jsonlite, parallel; at most 5 cores
#
# Outputs     : sim/results/{grid_a,grid_b,grid_c,planner_validation}.csv,
#               sim/results/grid_c_calibration.csv, sim/results/numbers.json,
#               sim/log.txt, sim/sessionInfo.txt
#
# Structure   : 1. Speed-up check (5 cores against 1)
#               2. Deployment example and planner headline numbers
#               3. Grid C calibration of the logit shift
#               4. All grid cells in parallel
#               5. Write the CSVs
#               6. Write numbers.json
# =============================================================================

require(structeval)
require(jsonlite)
require(parallel)

rm(list = ls())

ptm1 <- proc.time()   # measure execution time

args_vec <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args_vec) > 0L) as.integer(args_vec[1]) else 2000L
B <- 1000L          # bootstrap replicates inside each simulated data set
n_perm_b <- 2000L   # sign flips inside Grid B (deployment example uses 10,000)
n_perm_c <- 2000L   # sign flips inside Grid C
cores <- 5L

dir_sim <- "sim/"
dir_res <- paste0(dir_sim, "results/")
dir_cells <- paste0(dir_res, "cells_reps", reps, "/")
dir.create(dir_cells, showWarnings = FALSE, recursive = TRUE)
file_log <- paste0(dir_sim, "log.txt")

source(paste0(dir_sim, "cells.R"))

log_line <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", ...)
  cat(line, "\n", file = file_log, append = TRUE, sep = "")
  print(line)
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
log_line("START reps=", reps, " B=", B, " cores=", cores, " seed_base=",
  SEED_BASE)

# =============================================================================
# 1. Speed-up check (5 cores against 1)
# =============================================================================
print("==== 1. Speed-up check ====")
ptm_local <- proc.time()

# Ten identical Grid A list replicates batches, timed serially and forked.
job_speed <- function(i1) {
  set.seed(i1)
  for (i8 in seq_len(20L)) {
    mat <- sim_core(200, 0.9, 4, 0.2, 1, 1, "list")
    ab <- ab_of(mat, "list")
    eval_intervals(ab$a, ab$b, "list", B)
  }
  i1
}
time_1 <- system.time(lapply(seq_len(10L), job_speed))[["elapsed"]]
time_5 <- system.time(mclapply(seq_len(10L), job_speed,
  mc.cores = cores))[["elapsed"]]
speedup <- time_1 / time_5
log_line(sprintf("SPEEDUP 1 core %.2f s, %d cores %.2f s, ratio %.2f",
  time_1, cores, time_5, speedup))


time_section <- proc.time() - ptm_local
print("Section: 1. Speed-up check")
print(time_section)

# =============================================================================
# 2. Deployment example and planner headline numbers
# =============================================================================
print("==== 2. Deployment example ====")
ptm_local <- proc.time()

nums <- list()

# per-note table: each scalar note is one of cor/inc/mis/spu/tn
tab_t1 <- data.frame(
  key = c("diagnosis", "histology", "stage", "drugs", "biomarkers"),
  field = c("First Primary Diagnosis", "First Primary Histology",
    "First Primary Overall Stage", "Treatment Drugs",
    "Biomarker Test Results"),
  type = c("scalar", "scalar", "scalar", "list", "list"),
  cor = c(154, 124, 88, 352, 301), inc = c(3, 5, 0, 0, 0),
  mis = c(1, 1, 6, 40, 41), spu = c(2, 0, 8, 55, 16),
  tn = c(40, 70, 98, 0, 0), stringsAsFactors = FALSE
)
expand_field <- function(row) {
  outcome_vec <- c("cor", "inc", "mis", "spu", "tn")
  lab <- rep(outcome_vec, unlist(row[outcome_vec]))
  m <- vapply(outcome_vec, function(o) as.integer(lab == o),
    integer(length(lab)))
  data.frame(note = paste0(row$key, "_", seq_along(lab)), field = row$field,
    type = row$type, m, stringsAsFactors = FALSE)
}
counts <- do.call(rbind, lapply(seq_len(nrow(tab_t1)), function(i1) {
  expand_field(tab_t1[i1, ])
}))

published <- list(
  diagnosis = list(precision = c(94.1, 99.4), recall = c(94.8, 99.4),
    f1 = c(94.7, 99.1), specificity = c(87.8, 100)),
  histology = list(precision = c(92.4, 99.2), recall = c(91.4, 98.5),
    f1 = c(92.0, 98.8), specificity = c(100, 100)),
  stage = list(precision = c(85.7, 96.8), recall = c(88.2, 97.2),
    f1 = c(88.4, 96.2), specificity = c(87.0, 97.1)),
  drugs = list(precision = c(80.7, 91.6), recall = c(83.2, 95.4),
    f1 = c(83.4, 92.1)),
  biomarkers = list(precision = c(90.3, 98.4), recall = c(82.4, 92.8),
    f1 = c(87.3, 94.7))
)

for (i1 in seq_len(nrow(tab_t1))) {
  key_local <- tab_t1$key[i1]
  field_local <- tab_t1$field[i1]
  scalar_local <- tab_t1$type[i1] == "scalar"
  metric_vec <- names(published[[key_local]])
  for (metric_local in metric_vec) {
    pre <- paste0("deployment.", key_local, ".", metric_local, ".")
    pub <- published[[key_local]][[metric_local]] / 100
    nums[[paste0(pre, "published_lower")]] <- pub[1]
    nums[[paste0(pre, "published_upper")]] <- pub[2]
    if (scalar_local) {
      # Wilson for note-level proportions, cluster_bca for scalar F1
      if (metric_local == "f1") {
        res <- metrics(counts, field_local, "f1", B = 5000L, seed = SEED_BASE)
        tt <- test_threshold(counts, field_local, "f1", theta0 = 0.90,
          B = 5000L, seed = SEED_BASE)
        stopifnot(res$method == "cluster_bca", abs(tt$lower - res$lower) < 1e-12)
        nums[[paste0(pre, "cluster_bca_lower")]] <- res$lower
        nums[[paste0(pre, "cluster_bca_upper")]] <- res$upper
        nums[[paste0(pre, "boundary")]] <- res$boundary
      } else {
        res <- metrics(counts, field_local, metric_local)
        tt <- test_threshold(counts, field_local, metric_local, theta0 = 0.90)
        stopifnot(res$method == "wilson")
        nums[[paste0(pre, "wilson_lower")]] <- res$lower
        nums[[paste0(pre, "wilson_upper")]] <- res$upper
      }
      nums[[paste0(pre, "estimate")]] <- res$estimate
      nums[[paste0(pre, "t1_reject_theta90")]] <- isTRUE(tt$reject)
      nums[[paste0(pre, "t1_p_theta90")]] <- tt$p_value
    } else {
      # list: naive Wilson on pooled items; T1 decided by the published bound
      res <- metrics(counts, field_local, metric_local, method = "naive_wilson")
      nums[[paste0(pre, "estimate")]] <- res$estimate
      nums[[paste0(pre, "wilson_lower")]] <- res$lower
      nums[[paste0(pre, "wilson_upper")]] <- res$upper
      nums[[paste0(pre, "naive_wilson_lower")]] <- res$lower
      nums[[paste0(pre, "naive_wilson_upper")]] <- res$upper
      nums[[paste0(pre, "width_ratio")]] <- (res$upper - res$lower) /
        (pub[2] - pub[1])
      nums[[paste0(pre, "t1_reject_theta90")]] <- pub[1] > 0.90
    }
  }
  if (!scalar_local) {
    t3 <- test_asymmetry(counts, field_local, method = "binomial")
    pre <- paste0("deployment.", key_local, ".mis_share.")
    nums[[paste0(pre, "t3_binomial_p")]] <- t3$p_value
    nums[[paste0(pre, "t3_share_mis")]] <- t3$estimate
    nums[[paste0(pre, "estimate")]] <- t3$estimate
    nums[[paste0(pre, "t3_reject")]] <- isTRUE(t3$reject)
  }
}

# certifiable values and Wilson bounds quoted in the text
mc_vec <- min_certifiable(c(50, 100, 200), 0.90)
nums[["min_certifiable.n50.theta90"]] <- mc_vec[1]
nums[["min_certifiable.n100.theta90"]] <- mc_vec[2]
nums[["min_certifiable.n200.theta90"]] <- mc_vec[3]
nums[["wilson_lower.k49.n50"]] <- wilson_ci(49, 50, 0.95)[, "lower"]
nums[["wilson_lower.k194.n200"]] <- wilson_ci(194, 200, 0.95)[, "lower"]

# planner headline keys
for (pair_local in list(c(0.90, 0.95), c(0.90, 0.97), c(0.90, 0.98),
  c(0.85, 0.92))) {
  nums[[sprintf("planner.n_units.theta%d.theta%d", round(100 * pair_local[1]),
    round(100 * pair_local[2]))]] <- n_threshold(pair_local[1],
    pair_local[2])$n_units
}
nums[["planner.deff.m4.rho2"]] <- n_threshold(0.90, 0.95, m = 4,
  rho = 0.2)$deff
for (m_local in 1:4) {
  for (rho_local in c(0, 0.2, 0.5)) {
    nt <- n_threshold(0.90, 0.95, m = m_local, rho = rho_local)
    nums[[sprintf("planner.n_notes.list.theta90.theta95.m%d.rho%d", m_local,
      round(10 * rho_local))]] <- nt$n_notes
  }
}
stopifnot(nums[["planner.n_units.theta90.theta95"]] == 239,
  nums[["min_certifiable.n50.theta90"]] == 50)
log_line("DEPLOYMENT and planner keys computed (", length(nums), " keys)")


time_section <- proc.time() - ptm_local
print("Section: 2. Deployment example")
print(time_section)

# =============================================================================
# 3. Grid C calibration of the logit shift
# =============================================================================
print("==== 3. Grid C calibration ====")
ptm_local <- proc.time()

# The logit shift that gives the target metric-scale difference at n = 1e6;
# the B marginal does not depend on kappa, so one shift per type and delta,
# found with kappa = 1 (shared uniforms) to cut the noise of the difference.
tab_cal <- expand.grid(type = c("list", "proportion"), delta = c(0, 0.03, 0.05),
  stringsAsFactors = FALSE)
tab_cal$m <- ifelse(tab_cal$type == "list", 2, 1)
tab_cal$rho <- ifelse(tab_cal$type == "list", 0.2, 0)
tab_cal$delta_logit <- unlist(mclapply(seq_len(nrow(tab_cal)), function(i1) {
  calibrate_delta_logit(tab_cal$delta[i1], tab_cal$type[i1], tab_cal$m[i1],
    tab_cal$rho[i1], SEED_BASE)
}, mc.cores = cores))
write.csv(tab_cal, paste0(dir_res, "grid_c_calibration.csv"),
  row.names = FALSE)
for (i1 in seq_len(nrow(tab_cal))) {
  log_line(sprintf("CALIBRATE type=%s delta=%.2f delta_logit=%.6f",
    tab_cal$type[i1], tab_cal$delta[i1], tab_cal$delta_logit[i1]))
}


time_section <- proc.time() - ptm_local
print("Section: 3. Grid C calibration")
print(time_section)

# =============================================================================
# 4. All grid cells in parallel
# =============================================================================
print("==== 4. Grid cells ====")
ptm_local <- proc.time()

cells_a <- grid_a_cells()
cells_b <- grid_b_cells()
cells_c <- grid_c_cells()
cells_c$delta_logit <- tab_cal$delta_logit[match(
  paste(cells_c$type, cells_c$delta), paste(tab_cal$type, tab_cal$delta))]
scen <- planner_scenarios()
cells_p <- do.call(rbind, lapply(seq_len(nrow(scen)), function(i1) {
  nt <- n_threshold(scen$theta0[i1], scen$theta1[i1], m = scen$m[i1],
    rho = scen$rho[i1])
  data.frame(scen[rep(i1, 3), ], n_planned = nt$n_notes,
    factor = c(0.8, 1, 1.2), n_run = round(c(0.8, 1, 1.2) * nt$n_notes))
}))

# one job list with a global cell index (seed = SEED_BASE + index)
jobs <- c(
  lapply(seq_len(nrow(cells_a)), function(i1) list(g = "a", cell = cells_a[i1, ])),
  lapply(seq_len(nrow(cells_b)), function(i1) list(g = "b", cell = cells_b[i1, ])),
  lapply(seq_len(nrow(cells_c)), function(i1) list(g = "c", cell = cells_c[i1, ])),
  lapply(seq_len(nrow(cells_p)), function(i1) list(g = "p", cell = cells_p[i1, ]))
)
for (i1 in seq_along(jobs)) jobs[[i1]]$idx <- i1

# heaviest first: bootstrap cells by n, then the rest
cost_vec <- vapply(jobs, function(j) {
  c(a = 1, b = 0.15, c = 1.6, p = 0.8)[[j$g]] * j$cell$n_run %||% j$cell$n
}, numeric(1L))
order_vec <- order(-cost_vec)

run_job <- function(j) {
  file_local <- paste0(dir_cells, j$g, "_", sprintf("%03d", j$idx), ".rds")
  if (file.exists(file_local)) {
    return(readRDS(file_local))
  }
  out <- switch(j$g,
    a = run_cell_a(j$cell, j$idx, reps, B),
    b = run_cell_b(j$cell, j$idx, reps, n_perm_b),
    c = run_cell_c(j$cell, j$idx, reps, B, n_perm_c),
    p = run_cell_p(j$cell, j$idx, reps, B)
  )
  saveRDS(out, file_local)
  desc <- paste(names(j$cell), unlist(lapply(j$cell, as.character)),
    sep = "=", collapse = " ")
  log_line(sprintf("CELL %s%03d %.1f s %s", j$g, j$idx, out$seconds[1], desc))
  out
}

log_line("GRIDS begin: ", length(jobs), " cells (A ", nrow(cells_a), ", B ",
  nrow(cells_b), ", C ", nrow(cells_c), ", planner ", nrow(cells_p), ")")
t_grid <- proc.time()[["elapsed"]]
res_list <- mclapply(jobs[order_vec], run_job, mc.cores = cores,
  mc.preschedule = FALSE)
res_list[order_vec] <- res_list
failed <- vapply(res_list, function(x) inherits(x, "try-error") ||
  !is.data.frame(x), logical(1L))
if (any(failed)) {
  log_line("FAILED cells: ", paste(which(failed), collapse = ", "))
  print(res_list[failed][[1]])
  stop("Some cells failed; re-run to resume from the cache.")
}
grid_seconds <- proc.time()[["elapsed"]] - t_grid
log_line(sprintf("GRIDS end: %.1f s wall", grid_seconds))


time_section <- proc.time() - ptm_local
print("Section: 4. Grid cells")
print(time_section)

# =============================================================================
# 5. Write the CSVs
# =============================================================================
print("==== 5. Write the CSVs ====")
ptm_local <- proc.time()

g_vec <- vapply(jobs, function(j) j$g, character(1L))
tab_a <- do.call(rbind, res_list[g_vec == "a"])
tab_b <- do.call(rbind, res_list[g_vec == "b"])
tab_c <- do.call(rbind, res_list[g_vec == "c"])
tab_p <- do.call(rbind, res_list[g_vec == "p"])

cols_a <- c("type", "n", "theta", "m", "rho", "prevalence", "metric", "method",
  "reps", "B", "coverage", "coverage_se", "width", "boundary_rate", "size",
  "size_se", "power", "power_se", "theta_true_metric", "seconds")
cols_b <- c("n", "theta", "m", "rho", "spu_ratio", "method", "reps", "reject",
  "reject_se", "seconds")
cols_c <- c("type", "n", "theta_a", "delta", "kappa", "m", "rho", "method",
  "reps", "reject", "reject_se", "seconds")
cols_p <- c("scenario", "type", "theta0", "theta1", "m", "rho", "n_planned",
  "n_run", "power", "power_se")
write.csv(tab_a[, cols_a], paste0(dir_res, "grid_a.csv"), row.names = FALSE)
write.csv(tab_b[, cols_b], paste0(dir_res, "grid_b.csv"), row.names = FALSE)
write.csv(tab_c[, cols_c], paste0(dir_res, "grid_c.csv"), row.names = FALSE)
write.csv(tab_p[, cols_p], paste0(dir_res, "planner_validation.csv"),
  row.names = FALSE)

# undefined tests and realised true deltas go to the log, not the contracts
log_line("UNDEFINED A=", sum(tab_a$n_undefined), " B=", sum(tab_b$n_undefined),
  " C=", sum(tab_c$n_undefined), " planner=", sum(tab_p$n_undefined))
log_line(sprintf("THETA_TRUE max |theta_true_metric - theta| = %.2e",
  max(abs(tab_a$theta_true_metric - tab_a$theta))))
tab_cd <- unique(tab_c[, c("type", "n", "delta", "kappa", "delta_logit",
  "delta_true")])
log_line(sprintf("DELTA_TRUE range by target: %s", paste(tapply(
  tab_cd$delta_true, paste(tab_cd$type, tab_cd$delta), function(x) {
    sprintf("[%.4f, %.4f]", min(x), max(x))
  }), names(tapply(tab_cd$delta_true, paste(tab_cd$type, tab_cd$delta),
  length)), sep = " ", collapse = "; ")))


time_section <- proc.time() - ptm_local
print("Section: 5. Write the CSVs")
print(time_section)

# =============================================================================
# 6. Write numbers.json
# =============================================================================
print("==== 6. Write numbers.json ====")
ptm_local <- proc.time()

k100 <- function(x) sprintf("%d", as.integer(round(100 * x)))
k10 <- function(x) sprintf("%d", as.integer(round(10 * x)))
type_key <- function(t) ifelse(t == "list", "list", "prop")

# grid A
for (i1 in seq_len(nrow(tab_a))) {
  r <- tab_a[i1, ]
  stem <- paste0(type_key(r$type), ".n", r$n, ".theta", k100(r$theta),
    ".m", r$m, ".rho", k10(r$rho), if (r$prevalence != 1) {
      paste0(".prev", k10(r$prevalence))
    })
  boundary_se <- mc_se(r$boundary_rate, r$reps)
  q_vec <- c(coverage = r$coverage, coverage_se = r$coverage_se,
    width = r$width, width_se = r$width_se, boundary_rate = r$boundary_rate,
    boundary_rate_se = boundary_se, size = r$size, size_se = r$size_se,
    power = r$power, power_se = r$power_se)
  for (q in names(q_vec)) {
    nums[[paste0("a.", q, ".", stem, ".", r$method)]] <- q_vec[[q]]
  }
  nums[[paste0("a.theta_true_metric.", stem)]] <- r$theta_true_metric
}

# grid B
for (i1 in seq_len(nrow(tab_b))) {
  r <- tab_b[i1, ]
  stem <- paste0("list.n", r$n, ".theta", k100(r$theta), ".m", r$m, ".rho",
    k10(r$rho), ".r", k10(r$spu_ratio), ".", r$method)
  nums[[paste0("b.reject.", stem)]] <- r$reject
  nums[[paste0("b.reject_se.", stem)]] <- r$reject_se
}

# grid C
for (i1 in seq_len(nrow(tab_c))) {
  r <- tab_c[i1, ]
  stem <- paste0(type_key(r$type), ".n", r$n, ".theta", k100(r$theta_a),
    ".m", r$m, ".rho", k10(r$rho), ".delta", k100(r$delta), ".kappa",
    k10(r$kappa), ".", r$method)
  nums[[paste0("c.reject.", stem)]] <- r$reject
  nums[[paste0("c.reject_se.", stem)]] <- r$reject_se
}

# planner validation
for (i1 in seq_len(nrow(tab_p))) {
  r <- tab_p[i1, ]
  tag <- c("0.8" = "power_at_0p8n", "1" = "power_at_planned",
    "1.2" = "power_at_1p2n")[[as.character(r$factor)]]
  nums[[paste0("planner_validation.", tag, ".", r$scenario)]] <- r$power
  nums[[paste0("planner_validation.", tag, "_se.", r$scenario)]] <- r$power_se
  if (r$factor == 1) {
    nums[[paste0("planner_validation.n_planned.", r$scenario)]] <- r$n_planned
  }
}
at_planned <- tab_p[tab_p$factor == 1, ]
nums[["planner_validation.min_power"]] <- min(at_planned$power)
nums[["planner_validation.max_power"]] <- max(at_planned$power)

total_seconds <- (proc.time() - ptm1)[["elapsed"]]
nums[["compute.total_seconds"]] <- total_seconds
nums[["compute.grid_seconds"]] <- grid_seconds
nums[["compute.cores"]] <- cores
nums[["compute.speedup_5v1"]] <- speedup
nums[["compute.reps"]] <- reps
nums[["compute.B"]] <- B
nums[["compute.seed_base"]] <- SEED_BASE
nums[["compute.n_perm_b"]] <- n_perm_b
nums[["compute.n_perm_c"]] <- n_perm_c
nums[["compute.n_large"]] <- N_LARGE
nums[["compute.deployment_B"]] <- 5000L

writeLines(toJSON(nums, auto_unbox = TRUE, digits = NA, pretty = TRUE),
  paste0(dir_res, "numbers.json"))
writeLines(capture.output(sessionInfo()), paste0(dir_sim, "sessionInfo.txt"))
log_line(sprintf("END total %.1f s, %d keys", total_seconds, length(nums)))


time_section <- proc.time() - ptm_local
print("Section: 6. Write numbers.json")
print(time_section)

time_total <- proc.time() - ptm1
print("Total time")
print(time_total)
