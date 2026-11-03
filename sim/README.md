# Size-and-power simulation study

Run from the project root:

```
Rscript sim/run_sim.R 2000              # full study, about 14 min on 5 cores
Rscript sim/check_keys.R                # every manuscript key present in numbers.json
Rscript paper/figures/fig2_sizepower.R  # Figure 2 from the CSVs
```

Finished cells are cached in `sim/results/cells_reps<reps>/` and skipped on a
re-run. Delete that folder to recompute from scratch.

Run of 3 November 2026 with structeval 0.1.1: reps 2,000, B 1,000, 5 cores, 1,000 s in total
(grids 990 s). Measured speed-up of a 5-core batch over 1 core: 2.26 (`sim/log.txt`, line 734),
lower than the 3.65 of earlier runs because a package check ran at the same time.
Lines 1 to 366 of `sim/log.txt` hold the run of 29 October 2026 with structeval 0.1.0, and
lines 367 to 732 an earlier run of 3 November with a development version of 0.1.1.

## Design (spec v0.2, sections 5 to 9)

| Grid | Cells | Methods | Notes |
|---|---|---|---|
| A | 117 = 20 proportion (n 20 to 400) + 96 list + 1 prevalence check | wilson (= naive_wilson), percentile, cluster_bca for proportions; naive_wilson, percentile, cluster_bca, cluster_delta for lists | recall for proportions, F1 for lists |
| B | 162 | binomial, cluster_t, signflip | sign flips n_perm = 2,000 |
| C | 36 | cluster_delta, signflip, cluster_bca, mcnemar (proportions) | n_perm = 2,000, B = 1,000 |
| Planner | 12 scenarios x 3 n | routed T1 (wilson or cluster_bca) | `structeval::power_threshold()` |

Seed of each cell is 20261029 + its global index (A first, then B, C, planner).

**Engine.** Grid A draws one data set per replicate with `structeval:::.sim_core()` and
evaluates every method on it. Wilson, the cluster_delta standard error and the BCa
adjustment come from `structeval:::.wilson_ci()`, `.delta_se_logit()` and
`.bca_interval()`. One multinomial weight matrix (B = 1,000) per replicate gives
both the percentile and the BCa interval. The boundary rule and the boundary cap of
`ratio_ci()` are reproduced exactly. Size (θ0 = true value) and power (θ0 = true value - 0.05)
come from the same 95% interval. Grids B and C call `structeval:::.t3_core()` and
`.t2_core()` on one draw per replicate, with one sign matrix per cell. An
undefined test counts as no rejection.

**True values.** Grid A coverage is judged against `theta_true_metric`, one draw at
n = 10^6 per cell. Its largest distance from θ was 7.3e-4. Under r = 1 the
large-n F1 equals θ exactly, so this distance is draw noise. Grid C targets a
metric-scale Δ. The logit shift for each type and Δ was found by root search on
one n = 10^6 draw with κ = 1 (`results/grid_c_calibration.csv`). The realised
true Δ at n = 10^6 per cell lay within 0.0004 of its target (`log.txt`).

## Choices to confirm

1. **Grid C proportion cells run at m = 1, ρ = 0.** This follows the key contract
   (`prop ... m1.rho0`). Spec section 8 says ρ = 0.2 for Grid C.
2. **Prevalence check cell** (list, n 100, θ 0.90, m 4, ρ 0.2, prevalence 0.5).
   Its keys carry `.prev5` before the method.
3. **Planner scenarios.** P1 to P4 are proportions (0.90/0.95, 0.90/0.98,
   0.85/0.90, 0.85/0.93). L1 to L8 are lists: 0.90/0.95 at (m, ρ) = (2, 0), (2, 0.2),
   (4, 0.2), (4, 0.5); 0.85/0.93 at (2, 0.5), (4, 0); 0.85/0.90 at (4, 0.2);
   0.90/0.98 at (2, 0.2).
4. **Deployment list fields.** `t1_reject_theta90` is the published note-level lower
   bound compared with 0.90, as spec section 9 says. `wilson_lower/upper` for list
   fields hold the naive Wilson bounds (also written as `naive_wilson_*`).
   `width_ratio` is naive Wilson width / published width, per the key contract
   and the manuscript sentence. The vignette prints the inverse (published / naive).
5. **Scalar F1** uses `metrics()` with cluster_bca, B = 5,000, seed 20261029. T1 uses
   the same interval (asserted equal in the script).
6. **numbers.json** stores numbers at full precision and T1/T3 verdicts as JSON
   booleans. Extra keys: `*_se` for every rate, `a.theta_true_metric.*`,
   `planner_validation.power_at_0p8n.*` and `power_at_1p2n.*`,
   `planner.n_notes.list.theta90.theta95.m<m>.rho<r>` for the 12 points of panel (b).

## Findings that bear on the manuscript

- **T3 binomial is conservative here, not liberal.** Across the 54 Grid B null cells
  its size runs from 0.010 to 0.050 (median 0.024). In the deployment-like cell it is
  0.035, 0.023 and 0.018 at ρ = 0, 0.2 and 0.5. The manuscript (Results, Treatment
  Drugs paragraph) says the binomial test "rejected too often (Grid B)". Under this
  generative model that sentence is contradicted. Within a note, mis is binomial and
  spu Poisson, both driven by the same p_i, so the share is under-dispersed against
  the binomial.
- **The cluster intervals are mildly liberal one-sided.** In the deployment-like cell
  (list, n 100, θ 0.90, m 4), T1 size is 0.030 / 0.040 / 0.033 (cluster_bca) and
  0.027 / 0.039 / 0.033 (cluster_delta) at ρ = 0 / 0.2 / 0.5, nominal 0.025. Naive
  Wilson is 0.025 / 0.126 / 0.208.
- **The planner is liberal at ρ = 0.5.** L4 gives 0.719 and L5 gives 0.725 at the
  planned n. The DEFF of 1 + (m - 1)ρ understates the variance of the F1 ratio there.
  At ρ = 0 the list planner is very conservative (L1 0.984, L6 0.979).
- **Paired cluster_bca is liberal.** Grid C size runs from 0.053 to 0.077.
  Sign-flip and McNemar on proportions are very conservative (median size 0.018),
  because there are few discordant notes.
