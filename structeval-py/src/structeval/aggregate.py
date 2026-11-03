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


def fbeta_counts(hit: float, pos: float, act: float, beta: float = 1.0) -> float:
    """F-beta from counts, ``(1+b^2) hit / ((1+b^2) hit + b^2 (pos-hit) + (act-hit))``.

    ``pos`` is the number of reference items and ``act`` the number of predicted
    items. It equals :func:`fbeta` of precision and recall whenever both are
    defined, and is 0 (not ``nan``) when the model returns nothing for a field
    with reference items, as in the R package and scikit-learn's ``f1_score``.
    """
    b2 = beta * beta
    den = (1 + b2) * hit + b2 * (pos - hit) + (act - hit)
    return (1 + b2) * hit / den if den > 0 else np.nan


def _ratio(num: float, den: float) -> float:
    return num / den if den > 0 else np.nan


def binary_metrics(tp, fp, fn, tn, beta: float = 1.0) -> dict:
    """Precision/recall/F/accuracy/specificity from binary confusion totals."""
    precision = _ratio(tp, tp + fp)
    recall = _ratio(tp, tp + fn)
    return {
        "precision": precision,
        "recall": recall,
        f"f{beta:g}": fbeta_counts(tp, tp + fn, tp + fp, beta),
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
    pos, act = cor + inc + par + mis, cor + inc + par + spu
    return {
        "precision": _ratio(hit, act),
        "recall": _ratio(hit, pos),
        f"f{beta:g}": fbeta_counts(hit, pos, act, beta),
    }


def scalar_metrics(cor, inc, mis, spu, par, tn, beta: float = 1.0,
                   partial_weight: float = 0.5) -> dict:
    """:func:`set_metrics` plus specificity and accuracy for a scalar field.

    In a scalar field a spurious value occurs only on an empty reference, so
    specificity is ``tn / (tn + spu)`` and accuracy
    ``(cor + tn) / (cor + inc + mis + spu + par + tn)``. List fields have no
    specificity or accuracy (base paper Table 1).
    """
    out = set_metrics(cor, inc, mis, spu, par, beta, partial_weight)
    out["specificity"] = _ratio(tn, tn + spu)
    out["accuracy"] = _ratio(cor + tn, cor + inc + mis + spu + par + tn)
    return out


# counts that define each field kind, in display order
BINARY_COUNTS = ("TP", "FP", "FN", "TN")
SET_COUNTS = ("cor", "inc", "mis", "spu", "par", "tn")


def metrics_from_totals(totals: dict, kind: str, beta: float,
                        partial_weight: float) -> dict:
    """Dispatch summed counts to the right metric function."""
    if kind == "binary":
        return binary_metrics(totals["TP"], totals["FP"], totals["FN"],
                              totals["TN"], beta)
    if kind == "scalar":
        return scalar_metrics(totals["cor"], totals["inc"], totals["mis"],
                              totals["spu"], totals["par"], totals["tn"],
                              beta, partial_weight)
    return set_metrics(totals["cor"], totals["inc"], totals["mis"],
                       totals["spu"], totals["par"], beta, partial_weight)
