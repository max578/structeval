"""Per-case comparison of predicted vs. labelled structured fields.

Two field kinds:

* **binary** -- True/False detection, scored as a confusion matrix with an
  explicit policy for missing predictions.
* **set** -- a scalar is treated as a singleton set, a list as a set, so a
  single code path covers both. Counts are Correct / Incorrect / Missing /
  Spurious / Partial (the David-Batista named-entity convention).

Value conventions (the same as the R package ``structeval``):

* unlabelled (``None`` / NaN / a string that is blank after trimming white
  space)                                   -> case is skipped for that field
* no-information (``"-"``, an en or em dash, ``[]``, or a null token:
  ``null``, ``none``, ``na``, ``n/a`` in any case, so ``"NULL"``, ``"None"``,
  ``"NA"`` and ``"N/A"``)                 -> labelled-but-empty (evaluated)
  ``"unknown"`` and ``"NaN"`` stay values.
* anything else                            -> a value / set of values

Every item is normalised the same way on both sides: Unicode white space
trimmed and internal runs collapsed to one space; the simple one-to-one
lowercase mapping per character; a plain decimal (``-?(0|[1-9][0-9]*)
(\\.[0-9]+)?``) rewritten in one canonical form, so ``"4"`` == ``"4.0"``; every
other item, including ``"007"``, ``"1e3"``, ``"inf"`` and ``"NaN"``, is text.
"""
from __future__ import annotations

import math
import re

# --- value classification ------------------------------------------------

# The white-space characters of the R package (``.space_class``).
_SPACE = ("\t\n\x0b\x0c\r \x85\xa0\u1680\u2000-\u200a"
          "\u2028\u2029\u202f\u205f\u3000")
_SPACE_RUN = re.compile(f"[{_SPACE}]+")
_DECIMAL = re.compile(r"-?(0|[1-9][0-9]*)(\.[0-9]+)?")
DASHES = ("-", "\u2013", "\u2014")
NULL_TOKENS = ("null", "none", "na", "n/a")


_SPACE_ENDS = re.compile(f"^[{_SPACE}]+|[{_SPACE}]+$")


def _trim(x: str) -> str:
    return _SPACE_ENDS.sub("", x)


def _squish(x: str) -> str:
    return _SPACE_RUN.sub(" ", x).strip(" ")


def _lower(x: str) -> str:
    """Simple one-to-one lowercase mapping (``"\u0130"`` -> ``"i"``)."""
    return "".join((c.lower() or c)[0] for c in x)


def _canon_number(x: str) -> str:
    if not _DECIMAL.fullmatch(x):
        return x
    if "." in x:
        x = x.rstrip("0").rstrip(".")
    return "0" if x == "-0" else x


def _as_text(x) -> str | None:
    """Render a cell value as text; ``None`` for a missing value."""
    if x is None:
        return None
    if isinstance(x, bool):
        return "true" if x else "false"
    if isinstance(x, float):
        return None if math.isnan(x) else repr(x)
    return str(x)


def is_unlabeled(value) -> bool:
    """True when a case carries no label for this field (skip it)."""
    if isinstance(value, list):
        return False
    text = _as_text(value)
    return text is None or _squish(text) == ""


def normalize(x):
    """Normalised text of one item (see the module docstring)."""
    text = _as_text(x)
    if text is None:
        return None
    return _canon_number(_lower(_squish(text)))


def to_set(value, null_tokens=NULL_TOKENS) -> set:
    """Normalise any cell into a set of items (empty set for no-information)."""
    items = value if isinstance(value, list) else [value]
    empty = set(DASHES) | {normalize(t) for t in null_tokens}
    out = {normalize(v) for v in items}
    return {v for v in out if v not in empty and v not in (None, "")}


# Characters that, in a bare bracket list, show it is malformed.
_BARE_BAD = re.compile(r"['\"\[\]\\]")
_TOKEN = re.compile(r"'(?:[^'\\]|\\.)*'|\"(?:[^\"\\]|\\.)*\"")
_ESCAPE = re.compile(r"\\(.)", re.S)


def parse_listish(value):
    """Parse a bracketed list ``"['a', 'b']"`` into a list; else unchanged.

    Either every item is quoted (``'`` or ``"``, escapes ``\\'``, ``\\"`` and
    ``\\\\`` only) or none is (``"[a, b]"``). Blank items are dropped. Anything
    else raises ``ValueError``, as in the R reader.
    """
    if not isinstance(value, str):
        return value
    raw = _trim(value)
    if not raw.startswith("["):
        return value

    def fail(why):
        raise ValueError(f"Malformed bracket list '{raw}': {why}.")

    if not raw.endswith("]"):
        fail("text after the closing bracket, or no closing bracket")
    inner = _trim(raw[1:-1])
    if inner == "":
        return []
    tokens = _TOKEN.findall(inner)
    if not tokens:
        if _BARE_BAD.search(inner):
            fail("unbalanced quotes, brackets or backslashes in a bare list")
        return [t for t in (_trim(i) for i in inner.split(",")) if t]
    shape = re.sub(f"[{_SPACE}]+", "", _TOKEN.sub("Q", inner))
    if not re.fullmatch(r"Q(,Q)*", shape):
        fail("every item must be quoted and separated by one comma")
    items = []
    for t in tokens:
        body = t[1:-1]
        if any(e not in ("'", '"', "\\") for e in _ESCAPE.findall(body)):
            fail("only the escapes \\', \\\" and \\\\ are supported")
        item = _trim(_ESCAPE.sub(lambda m: m.group(1), body))
        if item:
            items.append(item)
    return items


# --- binary --------------------------------------------------------------

_TRUE = {"true", "t", "yes", "y", "1"}
_FALSE = {"false", "f", "no", "n", "0"}


def binary_state(x, null_tokens=NULL_TOKENS) -> str:
    """``"true"``, ``"false"``, ``"empty"`` (no value) or ``"other"``."""
    if is_unlabeled(x):
        return "empty"
    items = to_set(x, null_tokens)
    if not items:
        return "empty"
    if len(items) == 1:
        (item,) = items
        if item in _TRUE:
            return "true"
        if item in _FALSE:
            return "false"
    return "other"


def coerce_bool(x, null_tokens=NULL_TOKENS):
    """Map common truthy/falsey encodings to a bool; ``None`` otherwise."""
    state = binary_state(x, null_tokens)
    return {"true": True, "false": False}.get(state)


# how a *missing* prediction is scored against a present label
_MISSING_POLICIES = ("wrong", "abstain", "negative")


def compare_binary(expected, actual, missing: str = "wrong",
                   null_tokens=NULL_TOKENS) -> dict | None:
    """Confusion-matrix counts for one binary case, or ``None`` if skipped.

    A reference that is unlabelled, or holds no information (``"-"``, a dash
    or a null token), is skipped. A reference or prediction that is neither
    empty nor a recognised flag raises ``ValueError``, as in the R package.

    ``missing`` controls an empty prediction:

    * ``"wrong"``   -- missing on a positive label is a FN, on a negative a FP
                       (reproduces the original framework's behaviour).
    * ``"abstain"`` -- missing is excluded from all four counts (no credit, no
                       penalty); it shows up only in coverage (``n`` vs ``n_present``).
    * ``"negative"``-- missing is treated as a predicted ``False``.
    """
    t_state = binary_state(expected, null_tokens)
    if t_state == "empty":
        return None  # unlabelled or no information -> skip
    if t_state == "other":
        raise ValueError("Binary reference values must be TRUE or FALSE; "
                         f"got '{expected}'.")
    exp = t_state == "true"
    p_state = binary_state(actual, null_tokens)
    if p_state == "other":
        raise ValueError("Binary prediction values must be TRUE or FALSE; "
                         f"got '{actual}'. Recode them before scoring.")
    act = None if p_state == "empty" else p_state == "true"
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
                scalar: bool = False, null_tokens=NULL_TOKENS) -> dict | None:
    """Set-level counts for one case, or ``None`` if the label is unlabelled.

    ``parents`` maps a child value to its parent for partial-credit matching
    (predicting a parent when a child was expected -> Partial; predicting a
    child when the parent was expected -> Incorrect). When ``scalar`` is set and
    no parent map is in play, a singleton miss+spurious pair collapses to a
    single Incorrect, matching the scalar semantics of the original framework,
    and a cell with more than one item raises ``ValueError``. ``tn`` is 1 when
    both sides are empty.
    """
    if is_unlabeled(expected):
        return None
    if is_unlabeled(actual):
        actual = []

    parents = {normalize(k): normalize(v) for k, v in (parents or {}).items()}
    exp, act = to_set(expected, null_tokens), to_set(actual, null_tokens)
    if scalar and (len(exp) > 1 or len(act) > 1):
        raise ValueError("A scalar field holds more than one item: "
                         f"{sorted(exp | act)}.")

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
        "spu": len(spurious), "par": len(partial),
        "tn": int(not exp and not act), "present": present,
    }
