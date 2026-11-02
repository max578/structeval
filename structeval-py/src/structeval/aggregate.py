"""Turn per-case count totals into metrics (micro and macro).

All functions take *summed* counts so the same code computes the point estimate
and every bootstrap replicate. Undefined metrics (zero denominator) return
``nan`` rather than a misleading zero.
"""
from __future__ import annotations

import numpy as np


def fbeta(precision: float, recall: float, beta: float = 1.0) -> float:
    """F-beta from precision and recall; ``nan`` if either is undefined."""
    if np.isnan(precision) or np.isnan(recall):
        return np.nan
    if precision == 0 and recall == 0:
        return 0.0
    b2 = beta * beta
    denom = b2 * precision + recall
    return (1 + b2) * precision * recall / denom if denom > 0 else np.nan


def _ratio(num: float, den: float) -> float:
    return num / den if den > 0 else np.nan


def binary_metrics(tp, fp, fn, tn, beta: float = 1.0) -> dict:
    """Precision/recall/F/accuracy/specificity from binary confusion totals."""
    precision = _ratio(tp, tp + fp)
    recall = _ratio(tp, tp + fn)
    return {
        "precision": precision,
        "recall": recall,
        f"f{beta:g}": fbeta(precision, recall, beta),
        "accuracy": _ratio(tp + tn, tp + fp + fn + tn),
        "specificity": _ratio(tn, tn + fp),
    }


def set_metrics(cor, inc, mis, spu, par, beta: float = 1.0,
                partial_weight: float = 0.5) -> dict:
    """Precision/recall/F from set-comparison totals (David-Batista convention).

    A wrong value (``inc``) counts against both precision and recall; a partial
    match earns ``partial_weight`` (default 0.5) of a correct match.
    """
    hit = cor + partial_weight * par
    precision = _ratio(hit, cor + inc + par + spu)
    recall = _ratio(hit, cor + inc + par + mis)
    return {
        "precision": precision,
        "recall": recall,
        f"f{beta:g}": fbeta(precision, recall, beta),
    }


# counts that define each field kind, in display order
BINARY_COUNTS = ("TP", "FP", "FN", "TN")
SET_COUNTS = ("cor", "inc", "mis", "spu", "par")


def metrics_from_totals(totals: dict, kind: str, beta: float,
                        partial_weight: float) -> dict:
    """Dispatch summed counts to the right metric function."""
    if kind == "binary":
        return binary_metrics(totals["TP"], totals["FP"], totals["FN"],
                              totals["TN"], beta)
    return set_metrics(totals["cor"], totals["inc"], totals["mis"],
                       totals["spu"], totals["par"], beta, partial_weight)
