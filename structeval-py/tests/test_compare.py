import math

import numpy as np
import pytest

from structeval.compare import (coerce_bool, compare_binary, compare_set,
                                is_unlabeled, normalize, parse_listish, to_set)


# --- classification -------------------------------------------------------
@pytest.mark.parametrize("v,exp", [
    (None, True), ("", True), (np.nan, True), (float("nan"), True),
    ("-", False), ([], False), ("x", False), (0, False), (False, False),
])
def test_is_unlabeled(v, exp):
    assert is_unlabeled(v) is exp


def test_normalize_numeric_and_case():
    assert normalize("4") == 4.0
    assert normalize("  Female ") == "female"
    assert normalize(7) == 7


def test_to_set_unifies_scalar_and_list():
    assert to_set("Apple") == {"apple"}
    assert to_set(["Apple", "banana", "apple"]) == {"apple", "banana"}
    assert to_set("-") == set()
    assert to_set([]) == set()
    assert to_set("4") == {4.0}


def test_parse_listish():
    assert parse_listish("['a', 'b']") == ["a", "b"]
    assert parse_listish("not a list") == "not a list"
    assert parse_listish(5) == 5


# --- binary ---------------------------------------------------------------
def test_binary_confusion_basic():
    assert compare_binary(True, True) == {"TP": 1, "TN": 0, "FP": 0, "FN": 0, "present": 1}
    assert compare_binary(False, False)["TN"] == 1
    assert compare_binary(False, True)["FP"] == 1
    assert compare_binary(True, False)["FN"] == 1


def test_binary_unlabeled_skipped():
    assert compare_binary(None, True) is None
    assert compare_binary(np.nan, False) is None


def test_binary_robust_coercion():
    # the original mis-scored these; structeval coerces common encodings
    assert compare_binary(True, "true")["TP"] == 1
    assert compare_binary(True, 1)["TP"] == 1
    assert compare_binary(False, "no")["TN"] == 1
    assert coerce_bool("YES") is True
    assert coerce_bool("garbage") is None


def test_binary_missing_policies():
    # positive label, missing prediction
    assert compare_binary(True, None, missing="wrong")["FN"] == 1
    assert compare_binary(True, None, missing="abstain") == {
        "TP": 0, "TN": 0, "FP": 0, "FN": 0, "present": 0}
    assert compare_binary(True, None, missing="negative")["FN"] == 1
    # negative label, missing prediction
    assert compare_binary(False, None, missing="wrong")["FP"] == 1
    assert compare_binary(False, None, missing="negative")["TN"] == 1
    assert compare_binary(False, None, missing="abstain")["present"] == 0


# --- set / scalar ---------------------------------------------------------
def test_set_unlabeled_skipped():
    assert compare_set(None, ["a"]) is None
    assert compare_set("", "a") is None


def test_set_list_counts():
    r = compare_set(["apple", "banana"], ["apple", "cherry"])
    assert (r["cor"], r["mis"], r["spu"], r["inc"], r["par"]) == (1, 1, 1, 0, 0)


def test_scalar_mismatch_collapses_to_incorrect():
    # both present and different -> one Incorrect (not mis+spu)
    r = compare_set("Cancer", "Diabetes", scalar=True)
    assert (r["cor"], r["inc"], r["mis"], r["spu"]) == (0, 1, 0, 0)


def test_scalar_empty_label_with_prediction_is_spurious():
    r = compare_set("-", "Cancer", scalar=True)
    assert r["spu"] == 1 and r["mis"] == 0


def test_scalar_present_label_no_prediction_is_missing():
    r = compare_set("Cancer", "-", scalar=True)
    assert r["mis"] == 1 and r["spu"] == 0


def test_both_empty_is_all_zero():
    r = compare_set("-", "-", scalar=True)
    assert (r["cor"], r["inc"], r["mis"], r["spu"], r["par"]) == (0, 0, 0, 0, 0)
    assert r["present"] == 0


def test_parents_partial_and_incorrect():
    # predict parent when child expected -> Partial
    r = compare_set(["apple"], ["tree"], parents={"apple": "tree"})
    assert r["par"] == 1 and r["mis"] == 0 and r["spu"] == 0
    # predict child when parent expected -> Incorrect
    r2 = compare_set(["tree"], ["apple"], parents={"apple": "tree"})
    assert r2["inc"] == 1 and r2["spu"] == 0


def test_numeric_string_matches_number():
    r = compare_set("4", 4, scalar=True)
    assert r["cor"] == 1
