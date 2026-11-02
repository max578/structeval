# Size-and-power simulation study

Run from the repository root:

```
Rscript sim/run_sim.R 2000              # full study, about 14 min on 5 cores
Rscript sim/figures/fig2_sizepower.R    # Figure 2 from the result tables
```

Finished cells are cached in `sim/results/cells_reps<reps>/` and skipped on a
re-run. Delete that folder to recompute from scratch.

Run of 29 October 2026: reps 2,000, B 1,000, 5 cores, 830 s in total (grids 822 s).
Measured speed-up of a 5-core batch over 1 core: 3.46 (`sim/log.txt`, line 2).
`sim/run_console.txt` holds the console output of that run and
`sim/sessionInfo.txt` the R session.

## Design

| Grid | Cells | Methods | Notes |
|---|---|---|---|
| A | 117 = 20 proportion (n 20 to 400) + 96 list + 1 prevalence check | wilson (= naive_wilson), percentile and cluster_bca for proportions. naive_wilson, percentile, cluster_bca and cluster_delta for lists | recall for proportions, F1 for lists |
| B | 162 | binomial, cluster_t, signflip | sign flips n_perm = 2,000 |
| C | 36 | cluster_delta, signflip, cluster_bca, mcnemar (proportions) | n_perm = 2,000, B = 1,000 |
| Planner | 12 scenarios x 3 n | routed T1 (wilson or cluster_bca) | `structeval::power_threshold()` |

Seed of each cell is 20261029 + its global index (A first, then B, C, planner).

**Engine.** Grid A draws one data set per replicate with `structeval:::.sim_core()` and
evaluates every method on it. Wilson, the cluster_delta standard error and the BCa
adjustment come from `structeval:::.wilson_ci()`, `.delta_se_logit()` and
`.bca_interval()`. One multinomial weight matrix (B = 1,000) per replicate gives
both the percentile and the BCa interval. The boundary rule of `ratio_ci()` is
reproduced exactly. Size (θ0 = true value) and power (θ0 = true value - 0.05)
come from the same 95% interval. Grids B and C call `structeval:::.t3_core()` and
`.t2_core()` on one draw per replicate, with one sign matrix per cell. An
undefined test counts as no rejection.

**True values.** Grid A coverage is judged against `theta_true_metric`, one draw at
n = 10^6 per cell. Its largest distance from θ was 7.3e-4. Under r = 1 the
large-n F1 equals θ exactly, so this distance is draw noise. Grid C targets a
metric-scale Δ. The logit shift for each type and Δ was found by root search on
one n = 10^6 draw with κ = 1 (`results/grid_c_calibration.csv`). The realised
true Δ at n = 10^6 per cell lay within 0.0004 of its target (`log.txt`).

## Result files

| File | Contents |
|---|---|
| `results/grid_a.csv` | Coverage, size and power of every interval method, Grid A |
| `results/grid_b.csv` | Size and power of the asymmetry tests, Grid B |
| `results/grid_c.csv` | Size and power of the paired tests, Grid C |
| `results/grid_c_calibration.csv` | Logit shifts that realise each target Δ |
| `results/planner_validation.csv` | Realised power of the planners at 0.8, 1.0 and 1.2 times the planned n |
| `results/numbers.json` | Every reported number at full precision, keyed |
| `results/cells_reps2000/` | One cached result object per cell |

## Results in brief

- The binomial asymmetry test is conservative under this generative model. Across
  the 54 Grid B null cells its size runs from 0.010 to 0.050 (median 0.024). Within
  a note, omissions are binomial and fabrications Poisson, both driven by the same
  p_i, so their share is under-dispersed against the binomial.
- The cluster intervals are mildly liberal one-sided. In the deployment-like cell
  (list, n 100, θ 0.90, m 4), the certification test's size is 0.030 / 0.040 / 0.033
  (cluster_bca) and 0.027 / 0.039 / 0.033 (cluster_delta) at ρ = 0 / 0.2 / 0.5,
  against a nominal 0.025. Naive Wilson gives 0.025 / 0.126 / 0.208.
- The list-field planner is liberal at ρ = 0.5 (realised power 0.719 and 0.725 at
  the planned n in scenarios L4 and L5) and very conservative at ρ = 0 (0.984 and
  0.979 in L1 and L6). The design effect 1 + (m - 1)ρ understates the variance of
  the F1 ratio at high ρ.
- The paired cluster_bca test is liberal (Grid C size 0.053 to 0.077). Sign-flip and
  McNemar on proportions are very conservative (median size 0.018) because few notes
  are discordant.
