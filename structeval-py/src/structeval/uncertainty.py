"""Uncertainty quantification: closed-form Wilson intervals where the metric is
a genuine binomial proportion, BCa bootstrap everywhere else.

Design choices that make this "optimal yet simple":

* **Wilson score intervals** (exact ``z`` via :class:`statistics.NormalDist`,
  no SciPy) for binary precision / recall / accuracy / specificity. These have
  correct small-sample coverage and need no resampling.
* **BCa bootstrap** (bias-corrected and accelerated) for F-scores and all
  set-field metrics -- ratios that are not clean binomials. The reported point
  estimate is the **observed** statistic, never the bootstrap mean.
* **Vectorised**: every replicate is a row-index resample followed by a single
  column sum, so there is no per-draw group-by. Orders of magnitude faster than
  recomputing a full metric pipeline per draw.
* **Cluster bootstrap** (``clusters=``) resamples groups, not rows, for
  repeated-measures data (e.g. several documents per patient).
* Undefined replicates (zero denominator) are **counted and reported**, never
  silently dropped.
"""
from __future__ import annotations

from statistics import NormalDist

import numpy as np

_NORM = NormalDist()


def z_for(ci: float) -> float:
    """Two-sided ``z`` multiplier for a confidence level (e.g. 0.95 -> 1.96)."""
    return _NORM.inv_cdf(1.0 - (1.0 - ci) / 2.0)


def wilson(k: float, n: float, z: float) -> tuple[float, float]:
    """Wilson score interval for ``k`` successes in ``n`` Bernoulli trials."""
    if n <= 0:
        return (np.nan, np.nan)
    phat = k / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (phat + z2 / (2 * n)) / denom
    half = (z / denom) * np.sqrt(phat * (1 - phat) / n + z2 / (4 * n * n))
    return (max(0.0, center - half), min(1.0, center + half))


def _bca_bounds(boot: np.ndarray, theta: float, jack: np.ndarray,
                ci: float) -> tuple[float, float]:
    """BCa lower/upper percentiles; falls back to the plain percentile method
    if the bias/acceleration adjustment is numerically degenerate."""
    alpha = 1.0 - ci
    try:
        prop = float(np.mean(boot < theta))
        prop = min(max(prop, 1e-6), 1.0 - 1e-6)
        z0 = _NORM.inv_cdf(prop)

        jack = jack[~np.isnan(jack)]
        if jack.size >= 2:
            d = jack.mean() - jack
            den = 6.0 * (np.sum(d ** 2) ** 1.5)
            a = float(np.sum(d ** 3) / den) if den != 0 else 0.0
        else:
            a = 0.0

        def adj(q: float) -> float:
            zq = _NORM.inv_cdf(q)
            num = z0 + zq
            return _NORM.cdf(z0 + num / (1.0 - a * num))

        qlo, qhi = adj(alpha / 2.0), adj(1.0 - alpha / 2.0)
        if not (0.0 < qlo < qhi < 1.0):
            raise ValueError("degenerate BCa quantiles")
        return (float(np.percentile(boot, 100 * qlo)),
                float(np.percentile(boot, 100 * qhi)))
    except (ValueError, ZeroDivisionError, FloatingPointError):
        return (float(np.percentile(boot, 100 * (alpha / 2.0))),
                float(np.percentile(boot, 100 * (1.0 - alpha / 2.0))))


def _cluster_groups(clusters: np.ndarray) -> list[np.ndarray]:
    """Row indices grouped by cluster id, in first-appearance order."""
    order, seen = [], {}
    for i, c in enumerate(clusters):
        seen.setdefault(c, []).append(i)
        if c not in order:
            order.append(c)
    return [np.asarray(seen[c]) for c in order]


def bootstrap_ci(counts: np.ndarray, metric_fn, ci: float, n_boot: int,
                 rng: np.random.Generator,
                 clusters: np.ndarray | None = None) -> dict:
    """BCa bootstrap CI for ``metric_fn(column_totals)`` over per-case counts.

    Parameters
    ----------
    counts : (n_cases, k) array of per-case integer counts.
    metric_fn : callable mapping a length-``k`` totals vector to a scalar metric.
    """
    n = counts.shape[0]
    theta_hat = float(metric_fn(counts.sum(0)))

    boot = np.empty(n_boot, dtype=float)
    if clusters is None:
        for b in range(n_boot):
            idx = rng.integers(0, n, n)
            boot[b] = metric_fn(counts[idx].sum(0))
    else:
        groups = _cluster_groups(clusters)
        g = len(groups)
        for b in range(n_boot):
            chosen = rng.integers(0, g, g)
            idx = np.concatenate([groups[c] for c in chosen])
            boot[b] = metric_fn(counts[idx].sum(0))

    n_undef = int(np.isnan(boot).sum())
    boot = boot[~np.isnan(boot)]
    if np.isnan(theta_hat) or boot.size < 2:
        return {"estimate": theta_hat, "ci_low": np.nan, "ci_high": np.nan,
                "method": "bootstrap", "n_undefined_boot": n_undef}

    total = counts.sum(0)
    if clusters is None:
        jack = np.array([metric_fn(total - counts[i]) for i in range(n)])
    else:
        jack = np.array([metric_fn(total - counts[grp].sum(0))
                         for grp in _cluster_groups(clusters)])

    lo, hi = _bca_bounds(boot, theta_hat, jack, ci)
    return {"estimate": theta_hat, "ci_low": lo, "ci_high": hi,
            "method": "bca", "n_undefined_boot": n_undef}


def paired_difference_ci(counts_a: np.ndarray, counts_b: np.ndarray, metric_fn,
                         ci: float, n_boot: int, rng: np.random.Generator,
                         clusters: np.ndarray | None = None) -> dict:
    """BCa CI and two-sided bootstrap p-value for the *difference* in a metric
    between two models scored on the SAME cases (a paired comparison).

    The two count matrices are resampled with a shared index, so within-case
    correlation between the models is preserved -- the paired design that makes
    a difference detectable when the marginal intervals overlap.
    """
    n = counts_a.shape[0]
    diff = lambda ia: float(metric_fn(counts_a[ia].sum(0))
                            - metric_fn(counts_b[ia].sum(0)))
    theta = float(metric_fn(counts_a.sum(0)) - metric_fn(counts_b.sum(0)))

    boot = np.empty(n_boot, dtype=float)
    groups = _cluster_groups(clusters) if clusters is not None else None
    for b in range(n_boot):
        if groups is None:
            idx = rng.integers(0, n, n)
        else:
            chosen = rng.integers(0, len(groups), len(groups))
            idx = np.concatenate([groups[c] for c in chosen])
        boot[b] = diff(idx)

    boot = boot[~np.isnan(boot)]
    if np.isnan(theta) or boot.size < 2:
        return {"estimate": theta, "ci_low": np.nan, "ci_high": np.nan,
                "p_value": np.nan, "method": "paired-bootstrap"}

    all_idx = np.arange(n)
    if groups is None:
        jack = np.array([diff(np.delete(all_idx, i)) for i in range(n)])
    else:
        jack = np.array([diff(np.concatenate([groups[k] for k in range(len(groups)) if k != j]))
                         for j in range(len(groups))])
    lo, hi = _bca_bounds(boot, theta, jack, ci)
    # two-sided achieved significance level: how often replicates cross zero
    p = 2.0 * min(float(np.mean(boot <= 0)), float(np.mean(boot >= 0)))
    return {"estimate": theta, "ci_low": lo, "ci_high": hi,
            "p_value": min(1.0, p), "method": "paired-bca"}


# proportions that have a clean (k, n) binomial form -> Wilson, not bootstrap
def binary_wilson_cis(tp, fp, fn, tn, ci: float) -> dict:
    """Wilson CIs for the binomial binary metrics, keyed by metric name."""
    z = z_for(ci)
    specs = {
        "precision": (tp, tp + fp),
        "recall": (tp, tp + fn),
        "accuracy": (tp + tn, tp + fp + fn + tn),
        "specificity": (tn, tn + fp),
    }
    out = {}
    for name, (k, nn) in specs.items():
        lo, hi = wilson(k, nn, z)
        out[name] = {"ci_low": lo, "ci_high": hi, "method": "wilson"}
    return out
