"""structeval: a structure-aware evaluator for structured fields extracted
from clinical notes. Python reference implementation.

One entry point :func:`evaluate` consumes the same data contract as the
Oncoshot extraction framework (a label column ``F`` and a prediction column
``Res: F`` per field, with ``"-"`` meaning "no information" and blank meaning
"unlabelled") and returns a tidy :class:`Report`: one row per
field x confidence x metric, each carrying a point estimate and a confidence
interval computed by the right method (Wilson for binomial proportions, BCa
bootstrap for everything else).
"""
from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .aggregate import (BINARY_COUNTS, SET_COUNTS, binary_metrics,
                        scalar_metrics, set_metrics)
from .compare import (NULL_TOKENS, binary_state, compare_binary, compare_set,
                      is_unlabeled, parse_listish)
from .read import read_oncoshot
from .uncertainty import binary_wilson_cis, bootstrap_ci, paired_difference_ci

__version__ = "0.1.1"

__all__ = ["evaluate", "compare_models", "read_oncoshot", "Report",
           "__version__"]


# --------------------------------------------------------------------------- #
# Report container
# --------------------------------------------------------------------------- #
@dataclass
class Report:
    """Result of :func:`evaluate`.

    Attributes
    ----------
    metrics : tidy DataFrame -- field, confidence, metric, estimate, ci_low,
        ci_high, method, n, n_present, n_undefined_boot.
    cases : per-case comparison detail (labels, predictions, counts).
    params : the resolved evaluation parameters.
    """
    metrics: pd.DataFrame
    cases: pd.DataFrame
    params: dict

    def to_csv(self, folder: str, prefix: str = "") -> None:
        """Write ``metrics.csv`` and ``cases.csv`` to ``folder``."""
        os.makedirs(folder, exist_ok=True)
        p = f"{prefix} " if prefix else ""
        self.metrics.to_csv(os.path.join(folder, f"{p}metrics.csv"), index=False)
        self.cases.to_csv(os.path.join(folder, f"{p}cases.csv"))

    def pivot(self) -> pd.DataFrame:
        """Wide view: estimate per metric, indexed by field x confidence."""
        rows = self.metrics[self.metrics.method != "count"]
        return rows.pivot_table(index=["field", "confidence"], columns="metric",
                               values="estimate", sort=False)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _infer_fields(df: pd.DataFrame, prefix: str) -> list[str]:
    return [c[len(prefix):] for c in df.columns
            if c.startswith(prefix) and c[len(prefix):] in df.columns]


_KINDS = ("binary", "scalar", "list")


def _infer_kind(field: str, labels: pd.Series, preds: pd.Series) -> str:
    """Field type from the data, refusing every ambiguous case.

    A label column of booleans is binary. A label column whose values all read
    as flags (``0``/``1``, ``"yes"``/``"no"``, ``"TRUE"``/``"-"``) but are not all
    booleans, or that mixes booleans with other values, raises: declare the
    field in ``types`` instead. Otherwise the field is a list when any cell on
    either side is a list, and a scalar if not.
    """
    vals = [v for v in labels if not is_unlabeled(v)]
    is_bool = [isinstance(v, (bool, np.bool_)) for v in vals]
    if vals and all(is_bool):
        return "binary"
    flaglike = [binary_state(v) in ("true", "false", "empty") for v in vals]
    if any(is_bool) or (vals and all(flaglike)
                        and any(binary_state(v) != "empty" for v in vals)):
        raise ValueError(
            f"cannot tell whether field '{field}' is binary from its values; "
            "declare it with types={'" + field + "': 'binary'} (or 'scalar').")
    both = pd.concat([labels, preds])
    return "list" if both.map(lambda x: isinstance(x, list)).any() else "scalar"


def _field_kind(field, labels, preds, types):
    if types is None:
        return _infer_kind(field, labels, preds)
    if field not in types:
        raise ValueError(f"types has no entry for field '{field}'")
    kind = types[field]
    if kind not in _KINDS:
        raise ValueError(f"field types must be one of {_KINDS}; got '{kind}'")
    return kind


def _metric_fns(kind: str, beta: float, partial_weight: float):
    """Return {metric_name: fn(totals_vector)} and the count-column order."""
    fkey = f"f{beta:g}"
    if kind == "binary":
        order = BINARY_COUNTS  # TP, FP, FN, TN
        def make(name):
            return lambda t: binary_metrics(t[0], t[1], t[2], t[3], beta)[name]
        names = ["precision", "recall", fkey, "accuracy", "specificity"]
    elif kind == "scalar":
        order = SET_COUNTS  # cor, inc, mis, spu, par, tn
        def make(name):
            return lambda t: scalar_metrics(t[0], t[1], t[2], t[3], t[4], t[5],
                                           beta, partial_weight)[name]
        names = ["precision", "recall", fkey, "specificity", "accuracy"]
    else:  # list: no specificity or accuracy (base paper Table 1)
        order = SET_COUNTS
        def make(name):
            return lambda t: set_metrics(t[0], t[1], t[2], t[3], t[4],
                                        beta, partial_weight)[name]
        names = ["precision", "recall", fkey]
    return {n: make(n) for n in names}, order, fkey


# --------------------------------------------------------------------------- #
# public API
# --------------------------------------------------------------------------- #
def evaluate(df: pd.DataFrame, fields: list[str] | None = None,
             types: dict | None = None, pred_prefix: str = "Res: ", beta: float = 1.0, ci: float = 0.95,
             n_boot: int = 2000, cluster: str | None = None,
             missing: str = "wrong", partial_weight: float = 0.5,
             parents: dict | None = None,
             confidence_suffix: str = " confidence",
             null_tokens=NULL_TOKENS,
             random_state: int | None = 42) -> Report:
    """Evaluate predicted vs. labelled fields with point estimates and CIs.

    Parameters
    ----------
    df : DataFrame with a unique index; a label column ``F`` and a prediction
        column ``{pred_prefix}F`` for every field evaluated.
    fields : fields to evaluate; ``None`` auto-infers from ``{pred_prefix}*``.
    types : mapping field -> ``"binary"`` | ``"scalar"`` | ``"list"``, as in the
        R package. ``None`` infers the type from the data and raises when the
        inference is ambiguous (for example 0/1 or TRUE/FALSE/"-" labels).
    beta : F-beta weight (1.0 = F1, 2.0 = recall-weighted F2).
    ci : confidence level (e.g. 0.95).
    n_boot : bootstrap replicates for non-binomial metrics.
    cluster : column of cluster ids for a cluster bootstrap (repeated measures).
    missing : missing-prediction policy for binary fields
        (``"wrong"`` | ``"abstain"`` | ``"negative"``); see
        :func:`structeval.compare.compare_binary`.
    partial_weight : credit for a partial (parent) match.
    parents : child -> parent map enabling partial-credit matching.
    null_tokens : items that mean an empty value, matched after case folding
        (default ``"null"``, ``"none"``, ``"na"``, ``"n/a"``); ``()`` keeps
        them as values.
    """
    if not (0.0 < ci < 1.0):
        raise ValueError("ci must be in (0, 1)")
    if df.index.has_duplicates:
        raise ValueError("df index must be unique")

    df = df.copy()
    for col in df.columns:
        df[col] = df[col].map(parse_listish)

    if fields is None:
        fields = _infer_fields(df, pred_prefix)
        if not fields:
            raise ValueError(
                f"cannot infer fields: no '{pred_prefix}F' columns with a "
                "matching 'F' column")
    missing_cols = [f for f in fields if f not in df.columns
                    or f"{pred_prefix}{f}" not in df.columns]
    if missing_cols:
        raise ValueError(f"missing label/prediction columns for: {missing_cols}")
    if cluster is not None and cluster not in df.columns:
        raise ValueError(f"cluster column '{cluster}' not in df")

    rng = np.random.default_rng(random_state)
    cases = df.copy()
    metric_rows: list[dict] = []

    for field in fields:
        labels, preds = df[field], df[f"{pred_prefix}{field}"]
        kind = _field_kind(field, labels, preds, types)
        scalar = kind == "scalar"
        count_cols = BINARY_COUNTS if kind == "binary" else SET_COUNTS

        # --- per-case comparison ---
        rec_idx, recs, present, conf, clus = [], [], [], [], []
        conf_col = f"{pred_prefix}{field}{confidence_suffix}"
        has_conf = conf_col in df.columns
        for idx in df.index:
            exp, act = labels.loc[idx], preds.loc[idx]
            r = (compare_binary(exp, act, missing, null_tokens)
                 if kind == "binary"
                 else compare_set(exp, act, parents, scalar, null_tokens))
            for c in count_cols:
                cases.loc[idx, f"{field} [{c}]"] = r[c] if r else np.nan
            if r is None:
                continue
            rec_idx.append(idx)
            recs.append([r[c] for c in count_cols])
            present.append(r["present"])
            conf.append(df.loc[idx, conf_col] if has_conf else "Overall")
            if cluster is not None:
                clus.append(df.loc[idx, cluster])

        if not recs:
            continue
        counts = np.asarray(recs, dtype=float)
        present = np.asarray(present)
        conf = np.asarray(conf, dtype=object)
        clus = np.asarray(clus, dtype=object) if cluster is not None else None

        fns, order, fkey = _metric_fns(kind, beta, partial_weight)

        # --- aggregate over Overall + each confidence level ---
        levels = ["Overall"]
        if has_conf:
            levels += sorted({c for c in conf
                              if not (c is None or (isinstance(c, float) and np.isnan(c)))},
                             key=str)
        for level in levels:
            mask = (np.ones(len(counts), bool) if level == "Overall"
                    else conf == level)
            if not mask.any():
                continue
            sub = counts[mask]
            sub_clus = clus[mask] if clus is not None else None
            totals = sub.sum(0)
            n, n_present = int(mask.sum()), int(present[mask].sum())

            # count rows (no CI -- a count is not a rate)
            for j, c in enumerate(order):
                metric_rows.append(_row(field, level, c, totals[j], np.nan,
                                        np.nan, "count", n, n_present, 0))

            # Wilson for binomial binary proportions
            wil = (binary_wilson_cis(*totals, ci) if kind == "binary" else {})
            for name, fn in fns.items():
                est = fn(totals)
                if name in wil:
                    w = wil[name]
                    metric_rows.append(_row(field, level, name, est,
                                            w["ci_low"], w["ci_high"],
                                            "wilson", n, n_present, 0))
                else:  # F-scores and all set metrics -> BCa bootstrap
                    b = bootstrap_ci(sub, fn, ci, n_boot, rng, sub_clus)
                    metric_rows.append(_row(field, level, name, est,
                                            b["ci_low"], b["ci_high"],
                                            b["method"], n, n_present,
                                            b["n_undefined_boot"]))

    metrics = pd.DataFrame(metric_rows, columns=[
        "field", "confidence", "metric", "estimate", "ci_low", "ci_high",
        "method", "n", "n_present", "n_undefined_boot"])
    if not metrics["confidence"].astype(str).ne("Overall").any():
        pass  # keep the column even when only Overall, for a stable schema
    params = {"ci": ci, "beta": beta, "n_boot": n_boot, "missing": missing,
              "partial_weight": partial_weight, "cluster": cluster,
              "random_state": random_state}
    return Report(metrics=metrics, cases=cases, params=params)


def compare_models(df: pd.DataFrame, fields: list[str], a_prefix: str,
                   b_prefix: str, metric: str = "f1", beta: float = 1.0,
                   types: dict | None = None,
                   ci: float = 0.95, n_boot: int = 2000, cluster: str | None = None,
                   missing: str = "wrong", partial_weight: float = 0.5,
                   parents: dict | None = None,
                   random_state: int | None = 42) -> pd.DataFrame:
    """Paired comparison of two models scored on the same cases.

    Both models predict into their own prefixed columns (``{a_prefix}F`` and
    ``{b_prefix}F``) against the shared label column ``F``. Returns a tidy frame
    with, per field, the metric difference (A minus B), its BCa interval and a
    two-sided bootstrap p-value. A difference is significant when the interval
    excludes zero.
    """
    if df.index.has_duplicates:
        raise ValueError("df index must be unique")
    df = df.copy()
    for col in df.columns:
        df[col] = df[col].map(parse_listish)
    rng = np.random.default_rng(random_state)
    rows = []
    for field in fields:
        labels = df[field]
        kind = _field_kind(field, labels,
                           df.get(f"{a_prefix}{field}", labels), types)
        scalar = kind == "scalar"
        count_cols = BINARY_COUNTS if kind == "binary" else SET_COUNTS
        fns, _, fkey = _metric_fns(kind, beta, partial_weight)
        mkey = fkey if metric in ("f", "fbeta", fkey) else metric
        fn = fns[mkey]

        a_counts, b_counts, clus = [], [], []
        for idx in df.index:
            exp = labels.loc[idx]
            pa, pb = df.loc[idx, f"{a_prefix}{field}"], df.loc[idx, f"{b_prefix}{field}"]
            ra = (compare_binary(exp, pa, missing) if kind == "binary"
                  else compare_set(exp, pa, parents, scalar))
            rb = (compare_binary(exp, pb, missing) if kind == "binary"
                  else compare_set(exp, pb, parents, scalar))
            if ra is None or rb is None:
                continue  # paired: case must be labelled for both
            a_counts.append([ra[c] for c in count_cols])
            b_counts.append([rb[c] for c in count_cols])
            if cluster is not None:
                clus.append(df.loc[idx, cluster])
        if not a_counts:
            continue
        res = paired_difference_ci(
            np.asarray(a_counts, float), np.asarray(b_counts, float), fn, ci,
            n_boot, rng, np.asarray(clus, object) if cluster is not None else None)
        rows.append({"field": field, "metric": mkey, "n": len(a_counts),
                     "diff_A_minus_B": res["estimate"], "ci_low": res["ci_low"],
                     "ci_high": res["ci_high"], "p_value": res["p_value"],
                     "significant": bool(res["ci_low"] > 0 or res["ci_high"] < 0)})
    return pd.DataFrame(rows)


def _row(field, conf, metric, est, lo, hi, method, n, n_present, n_undef):
    return {"field": field, "confidence": conf, "metric": metric,
            "estimate": float(est) if est is not None else np.nan,
            "ci_low": lo, "ci_high": hi, "method": method,
            "n": n, "n_present": n_present, "n_undefined_boot": n_undef}
