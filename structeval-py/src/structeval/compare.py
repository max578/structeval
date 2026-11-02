"""Per-case comparison of predicted vs. labelled structured fields.

Two field kinds:

* **binary** -- True/False detection, scored as a confusion matrix with an
  explicit policy for missing predictions.
* **set** -- a scalar is treated as a singleton set, a list as a set, so a
  single code path covers both. Counts are Correct / Incorrect / Missing /
  Spurious / Partial (the David-Batista named-entity convention).

Value conventions (matching the data contract of the original framework):

* unlabelled (``None`` / ``""`` / ``NaN``)  -> case is skipped for that field
* no-information (``"-"`` or ``[]``)         -> labelled-but-empty (evaluated)
* anything else                               -> a value / set of values
"""
from __future__ import annotations

import math
from ast import literal_eval

import pandas as pd

# --- value classification ------------------------------------------------

def is_unlabeled(value) -> bool:
    """True when a case carries no label for this field (skip it)."""
    if value is None or value == "":
        return True
    if isinstance(value, float):
        return math.isnan(value)
    return False


def _is_empty_scalar(value) -> bool:
    """True for an empty / no-information scalar (``None``/``""``/``"-"``/NaN)."""
    if value in (None, "", "-"):
        return True
    if isinstance(value, float):
        return math.isnan(value)
    return False


def normalize(x):
    """Case-fold strings and coerce numeric strings to float; pass others through.

    ``"4"`` -> ``4.0`` so a string label and a numeric prediction compare equal;
    ``"Female"`` -> ``"female"`` for case-insensitive matching.
    """
    if isinstance(x, str):
        x = x.strip()
        try:
            return float(x)
        except ValueError:
            return x.casefold()
    return x


def to_set(value) -> set:
    """Normalise any cell into a set of values (empty set for no-information)."""
    if isinstance(value, list):
        return {normalize(v) for v in value if not _is_empty_scalar(v)}
    if _is_empty_scalar(value):
        return set()
    return {normalize(value)}


def parse_listish(value):
    """Parse a stringified list ``"['a', 'b']"`` into a real list; else unchanged."""
    if isinstance(value, str) and value.strip().startswith("["):
        try:
            return literal_eval(value)
        except (ValueError, SyntaxError):
            return value
    return value


# --- binary --------------------------------------------------------------

_TRUE = {True, 1, "true", "yes", "y", "1", "t"}
_FALSE = {False, 0, "false", "no", "n", "0", "f"}


def coerce_bool(x):
    """Map common truthy/falsey encodings to a bool; return ``None`` if unknown."""
    if isinstance(x, bool):
        return x
    if isinstance(x, float) and math.isnan(x):
        return None
    key = x.strip().lower() if isinstance(x, str) else x
    if key in _TRUE:
        return True
    if key in _FALSE:
        return False
    return None


# how a *missing* prediction is scored against a present label
_MISSING_POLICIES = ("wrong", "abstain", "negative")


def compare_binary(expected, actual, missing: str = "wrong") -> dict | None:
    """Confusion-matrix counts for one binary case, or ``None`` if unlabelled.

    ``missing`` controls a prediction that could not be coerced to a bool:

    * ``"wrong"``   -- missing on a positive label is a FN, on a negative a FP
                       (reproduces the original framework's behaviour).
    * ``"abstain"`` -- missing is excluded from all four counts (no credit, no
                       penalty); it shows up only in coverage (``n`` vs ``n_present``).
    * ``"negative"``-- missing is treated as a predicted ``False``.
    """
    exp = coerce_bool(expected)
    if exp is None:
        return None  # unlabelled -> skip

    act = coerce_bool(actual)
    if act is None:
        if missing == "abstain":
            return {"TP": 0, "TN": 0, "FP": 0, "FN": 0, "present": 0}
        if missing == "negative":
            act = False
        # "wrong": leave act as None -> handled below as "not the label"
    present = 1
    tp = int(exp is True and act is True)
    tn = int(exp is False and act is False)
    fn = int(exp is True and act is not True)
    fp = int(exp is False and act is not False)
    return {"TP": tp, "TN": tn, "FP": fp, "FN": fn, "present": present}


# --- set / scalar --------------------------------------------------------

def compare_set(expected, actual, parents: dict | None = None,
                scalar: bool = False) -> dict | None:
    """Set-level counts for one case, or ``None`` if the label is unlabelled.

    ``parents`` maps a child value to its parent for partial-credit matching
    (predicting a parent when a child was expected -> Partial; predicting a
    child when the parent was expected -> Incorrect). When ``scalar`` is set and
    no parent map is in play, a singleton miss+spurious pair collapses to a
    single Incorrect, matching the scalar semantics of the original framework.
    """
    if is_unlabeled(expected):
        return None

    parents = {k.casefold(): v.casefold() for k, v in (parents or {}).items()}
    exp, act = to_set(expected), to_set(actual)

    correct = exp & act
    exp_parents = {parents[s] for s in exp if s in parents}
    act_parents = {parents[s] for s in act if s in parents}

    missing, partial = [], []
    for x in exp - act - act_parents:
        if x in parents and parents[x] in act:
            partial.append(x)
        else:
            missing.append(x)

    spurious, incorrect = [], []
    for x in act - exp - exp_parents:
        if x in parents and parents[x] in exp:
            incorrect.append(x)
        else:
            spurious.append(x)

    if scalar and not parents and len(missing) == 1 and len(spurious) == 1:
        incorrect, missing, spurious = spurious, [], []

    present = 1 if exp else 0  # label carried information
    return {
        "cor": len(correct), "inc": len(incorrect), "mis": len(missing),
        "spu": len(spurious), "par": len(partial), "present": present,
    }
