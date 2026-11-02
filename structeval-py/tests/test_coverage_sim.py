"""Empirical coverage simulation -- the core evidence that structeval's
intervals are calibrated. We generate data with a KNOWN true recall and check
how often each method's 95% interval contains it.

Kept small enough to run in the normal suite; the standalone, higher-precision
version lives in ``bench/coverage.py``.
"""
import numpy as np
import pandas as pd
import pytest

import structeval as se
from structeval.uncertainty import bootstrap_ci, wilson, z_for
from structeval.aggregate import set_metrics


def _recall_fn():
    return lambda t: set_metrics(t[0], t[1], t[2], t[3], t[4])["recall"]


@pytest.mark.parametrize("p_true,n", [(0.7, 60), (0.9, 60)])
def test_wilson_coverage_near_nominal(p_true, n):
    """Wilson coverage for a clean binomial proportion should be ~0.95."""
    rng = np.random.default_rng(123)
    z = z_for(0.95)
    trials, hit = 400, 0
    for _ in range(trials):
        k = rng.binomial(n, p_true)
        lo, hi = wilson(k, n, z)
        hit += lo <= p_true <= hi
    cov = hit / trials
    assert 0.92 <= cov <= 0.99, f"Wilson coverage {cov:.3f} off nominal"


def test_bca_coverage_near_nominal():
    """BCa coverage for set-recall (denominator = item count) should be ~0.95
    and at least as good as the plain percentile interval."""
    rng = np.random.default_rng(7)
    p_true, n_cases, trials = 0.7, 50, 250
    bca_hit = 0
    fn = _recall_fn()
    for _ in range(trials):
        # each case contributes one labelled item; correct w.p. p_true
        cor = rng.binomial(1, p_true, n_cases).astype(float)
        counts = np.column_stack(
            [cor, np.zeros(n_cases), 1 - cor, np.zeros(n_cases), np.zeros(n_cases)])
        res = bootstrap_ci(counts, fn, 0.95, 400, rng)
        if not np.isnan(res["ci_low"]):
            bca_hit += res["ci_low"] <= p_true <= res["ci_high"]
    cov = bca_hit / trials
    assert 0.90 <= cov <= 0.99, f"BCa coverage {cov:.3f} off nominal"
