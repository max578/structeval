import math

import numpy as np
import pytest

from structeval.aggregate import binary_metrics, fbeta, set_metrics


def test_fbeta_formulae():
    assert fbeta(0.8, 0.6, 1.0) == pytest.approx(2 * 0.8 * 0.6 / (0.8 + 0.6))
    assert fbeta(0.8, 0.6, 2.0) == pytest.approx(5 * 0.8 * 0.6 / (4 * 0.8 + 0.6))


def test_fbeta_zero_and_nan():
    assert fbeta(0.0, 0.0, 1.0) == 0.0
    assert math.isnan(fbeta(np.nan, 0.5, 1.0))
    assert math.isnan(fbeta(0.5, np.nan, 1.0))


def test_binary_metrics_values():
    m = binary_metrics(tp=8, fp=2, fn=4, tn=6, beta=1.0)
    assert m["precision"] == pytest.approx(0.8)
    assert m["recall"] == pytest.approx(8 / 12)
    assert m["accuracy"] == pytest.approx(14 / 20)
    assert m["specificity"] == pytest.approx(6 / 8)


def test_binary_undefined_recall_is_nan():
    m = binary_metrics(tp=0, fp=3, fn=0, tn=0)
    assert m["precision"] == 0.0
    assert math.isnan(m["recall"])
    # 0.1.1: F1 from counts, 2TP/(2TP+FP+FN), as in the R package
    assert m["f1"] == 0.0


def test_set_metrics_incorrect_counts_in_both_denominators():
    m = set_metrics(cor=3, inc=2, mis=1, spu=1, par=0)
    assert m["precision"] == pytest.approx(3 / 6)
    assert m["recall"] == pytest.approx(3 / 6)


def test_set_metrics_partial_weight():
    m = set_metrics(cor=0, inc=0, mis=0, spu=0, par=2, partial_weight=0.5)
    assert m["precision"] == pytest.approx(0.5)
    m2 = set_metrics(cor=0, inc=0, mis=0, spu=0, par=2, partial_weight=1.0)
    assert m2["precision"] == pytest.approx(1.0)


def test_set_metrics_zero_denominator_is_nan():
    m = set_metrics(cor=0, inc=0, mis=0, spu=0, par=0)
    assert math.isnan(m["precision"]) and math.isnan(m["recall"])
