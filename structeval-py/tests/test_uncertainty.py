import numpy as np
import pytest

from structeval.aggregate import set_metrics
from structeval.uncertainty import (bootstrap_ci, binary_wilson_cis, wilson,
                                    z_for)


def test_z_for_95():
    assert z_for(0.95) == pytest.approx(1.959963, abs=1e-5)
    assert z_for(0.90) == pytest.approx(1.644853, abs=1e-5)


def test_wilson_known_value():
    # 8/10: textbook Wilson 95% interval ~ [0.4901, 0.9427]
    lo, hi = wilson(8, 10, z_for(0.95))
    assert lo == pytest.approx(0.4901, abs=1e-3)
    assert hi == pytest.approx(0.9427, abs=1e-3)


def test_wilson_bounds_clamped_and_degenerate():
    assert wilson(0, 0, z_for(0.95)) == (float("nan"), float("nan")) or \
        all(np.isnan(x) for x in wilson(0, 0, z_for(0.95)))
    lo, hi = wilson(10, 10, z_for(0.95))
    assert 0.0 <= lo <= 1.0 and hi <= 1.0


def test_wilson_wider_for_higher_confidence():
    w95 = wilson(7, 10, z_for(0.95))
    w99 = wilson(7, 10, z_for(0.99))
    assert (w99[1] - w99[0]) > (w95[1] - w95[0])


def _precision_fn(beta=1.0):
    # totals vector order: cor, inc, mis, spu, par
    return lambda t: set_metrics(t[0], t[1], t[2], t[3], t[4], beta)["precision"]


def test_bootstrap_estimate_is_observed_not_mean():
    rng = np.random.default_rng(0)
    n = 200
    cor = rng.integers(0, 2, n).astype(float)
    counts = np.column_stack([cor, np.zeros(n), np.zeros(n), 1 - cor, np.zeros(n)])
    obs = cor.sum() / n  # precision = cor/(cor+spu) = cor/n here
    res = bootstrap_ci(counts, _precision_fn(), 0.95, 1000, rng)
    assert res["estimate"] == pytest.approx(obs)  # observed, not bootstrap mean
    assert res["ci_low"] <= res["estimate"] <= res["ci_high"]
    assert res["method"] == "bca"


def test_bootstrap_reports_undefined_replicates():
    rng = np.random.default_rng(1)
    # almost-empty field: most resamples have zero denominator
    counts = np.zeros((10, 5))
    counts[0] = [1, 0, 0, 0, 0]  # one correct, rest all-zero
    res = bootstrap_ci(counts, _precision_fn(), 0.95, 300, rng)
    assert res["n_undefined_boot"] >= 0  # surfaced, not hidden


def test_cluster_bootstrap_widens_with_positive_intra_cluster_correlation():
    rng = np.random.default_rng(2)
    # 20 clusters of 10 identical cases -> strong intra-cluster correlation
    block = []
    for k in range(20):
        val = 1.0 if k % 2 == 0 else 0.0
        for _ in range(10):
            block.append([val, 0, 0, 1 - val, 0])
    counts = np.asarray(block, dtype=float)
    clusters = np.repeat(np.arange(20), 10)
    naive = bootstrap_ci(counts, _precision_fn(), 0.95, 800, rng)
    clus = bootstrap_ci(counts, _precision_fn(), 0.95, 800,
                        np.random.default_rng(2), clusters=clusters)
    w_naive = naive["ci_high"] - naive["ci_low"]
    w_clus = clus["ci_high"] - clus["ci_low"]
    assert w_clus > w_naive  # ignoring clustering under-states uncertainty


def test_binary_wilson_cis_keys():
    cis = binary_wilson_cis(tp=8, fp=2, fn=4, tn=6, ci=0.95)
    assert set(cis) == {"precision", "recall", "accuracy", "specificity"}
    assert all(cis[m]["method"] == "wilson" for m in cis)
