import math

import numpy as np
import pandas as pd
import pytest

import structeval as se


@pytest.fixture
def samples():
    return pd.read_csv("samples.csv", index_col="Patient ID")


def _est(rep, field, metric, conf="Overall"):
    m = rep.metrics
    return m[(m.field == field) & (m.confidence == conf)
            & (m.metric == metric)].estimate.iloc[0]


# --- parity with the original framework (hard-coded expected values) ------
PARITY = {
    ("Has metastasis", "precision"): 0.5,
    ("Has metastasis", "recall"): 0.5,
    ("Has metastasis", "f1"): 0.5,
    ("Has metastasis", "accuracy"): 0.5,
    ("First Primary Diagnosis", "precision"): 1 / 3,
    ("First Primary Diagnosis", "recall"): 1 / 3,
    ("Treatment Drugs", "precision"): 0.5,
    ("Treatment Drugs", "recall"): 0.5,
    ("Test Results", "recall"): 0.6,
}


def test_point_estimate_parity(samples):
    rep = se.evaluate(samples, n_boot=200, random_state=0)
    for (field, metric), expected in PARITY.items():
        assert _est(rep, field, metric) == pytest.approx(expected), (field, metric)


def test_count_parity(samples):
    rep = se.evaluate(samples, n_boot=10, random_state=0)
    assert _est(rep, "Treatment Drugs", "cor") == 2
    assert _est(rep, "Treatment Drugs", "mis") == 2
    assert _est(rep, "Treatment Drugs", "spu") == 2
    assert _est(rep, "First Primary Diagnosis", "inc") == 1
    for c, val in [("TP", 1), ("FP", 1), ("FN", 1), ("TN", 1)]:
        assert _est(rep, "Has metastasis", c) == val


# --- structural / API behaviour ------------------------------------------
def test_field_inference(samples):
    rep = se.evaluate(samples, n_boot=10)
    assert set(rep.metrics.field.unique()) == {
        "Has metastasis", "First Primary Diagnosis", "First Primary Histology",
        "Treatment Drugs", "Test Results"}


def test_counts_have_no_ci(samples):
    rep = se.evaluate(samples, n_boot=10)
    counts = rep.metrics[rep.metrics.method == "count"]
    assert counts["ci_low"].isna().all() and counts["ci_high"].isna().all()


def test_rate_rows_have_ci_and_ordering(samples):
    rep = se.evaluate(samples, n_boot=200, random_state=1)
    rates = rep.metrics[rep.metrics.method != "count"]
    ok = rates.dropna(subset=["ci_low", "ci_high", "estimate"])
    assert ((ok.ci_low <= ok.estimate + 1e-9) & (ok.estimate <= ok.ci_high + 1e-9)).all()


def test_binary_uses_wilson_set_uses_bca(samples):
    rep = se.evaluate(samples, n_boot=50)
    m = rep.metrics
    binr = m[(m.field == "Has metastasis") & (m.metric == "precision")].method.iloc[0]
    setr = m[(m.field == "Treatment Drugs") & (m.metric == "precision")].method.iloc[0]
    assert binr == "wilson" and setr == "bca"


def test_unique_index_enforced():
    df = pd.DataFrame({"f": [True, False], "Res: f": [True, True]}, index=[1, 1])
    with pytest.raises(ValueError, match="unique"):
        se.evaluate(df)


def test_ci_range_validation(samples):
    with pytest.raises(ValueError, match="ci must be"):
        se.evaluate(samples, ci=1.5)


def test_partial_labels_per_field():
    # field b is labelled only for some rows; should not raise and uses labelled rows
    df = pd.DataFrame({
        "a": [True, False, True, False],
        "Res: a": [True, False, False, False],
        "b": ["x", np.nan, "y", np.nan],
        "Res: b": ["x", "z", "y", "w"],
    })
    rep = se.evaluate(df, n_boot=50)
    assert _est(rep, "b", "cor") == 2  # only the two labelled rows count
    assert int(rep.metrics[(rep.metrics.field == "b")
               & (rep.metrics.metric == "cor")].n.iloc[0]) == 2


def test_confidence_stratification():
    df = pd.DataFrame({
        "x": ["a", "b", "c", "d"],
        "Res: x": ["a", "x", "c", "y"],
        "Res: x confidence": ["High", "Low", "High", "Low"],
    })
    rep = se.evaluate(df, n_boot=50)
    confs = set(rep.metrics.confidence.unique())
    assert {"Overall", "High", "Low"} <= confs
    # High subset: both correct -> precision 1.0
    assert _est(rep, "x", "precision", conf="High") == pytest.approx(1.0)


def test_missing_policy_changes_binary_score():
    df = pd.DataFrame({"f": [True, False], "Res: f": [np.nan, np.nan]})
    wrong = se.evaluate(df, missing="wrong", n_boot=10)
    abst = se.evaluate(df, missing="abstain", n_boot=10)
    # wrong: TP=0,FN=1 (pos), FP=1 (neg) -> precision 0; abstain: all zero counts
    assert _est(wrong, "f", "FN") == 1 and _est(wrong, "f", "FP") == 1
    assert _est(abst, "f", "FN") == 0 and _est(abst, "f", "FP") == 0


def test_beta_two(samples):
    rep = se.evaluate(samples, beta=2.0, n_boot=10)
    assert "f2" in set(rep.metrics.metric.unique())


def test_pivot_and_to_csv(tmp_path, samples):
    rep = se.evaluate(samples, n_boot=10)
    piv = rep.pivot()
    assert "precision" in piv.columns
    rep.to_csv(str(tmp_path), prefix="run1")
    assert (tmp_path / "run1 metrics.csv").exists()
    assert (tmp_path / "run1 cases.csv").exists()


def test_compare_models_detects_paired_difference():
    # 100 cases; model A correct on 85, model B on 70, with A dominating B
    rng = np.random.default_rng(0)
    n = 100
    truth = ["x"] * n
    a = ["x" if i < 85 else "z" for i in range(n)]
    b = ["x" if i < 70 else "z" for i in range(n)]
    df = pd.DataFrame({"d": truth, "A: d": a, "B: d": b})
    out = se.compare_models(df, ["d"], "A: ", "B: ", metric="f1", n_boot=1000)
    row = out.iloc[0]
    assert row["diff_A_minus_B"] > 0
    assert "p_value" in out.columns and 0.0 <= row["p_value"] <= 1.0


def test_compare_models_no_false_difference():
    # identical predictions -> diff 0, not significant
    df = pd.DataFrame({"d": ["x", "y", "z", "x"], "A: d": ["x", "y", "q", "x"],
                       "B: d": ["x", "y", "q", "x"]})
    out = se.compare_models(df, ["d"], "A: ", "B: ", n_boot=500)
    assert out.iloc[0]["diff_A_minus_B"] == 0
    assert not out.iloc[0]["significant"]
