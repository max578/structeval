"""Scoring rules shared with the R package (adversarial review SC-01..SC-18,
K1, K2). Each test names the finding it closes."""
import io
import math

import numpy as np
import pandas as pd
import pytest

import structeval as se
from structeval.aggregate import binary_metrics, set_metrics
from structeval.compare import (compare_binary, compare_set, is_unlabeled,
                                normalize, parse_listish, to_set)


def counts(t, p, scalar=True):
    r = compare_set(t, p, scalar=scalar)
    return (r["cor"], r["inc"], r["mis"], r["spu"], r["tn"])


def write(tmp_path, text, name="x.csv", raw=None):
    f = tmp_path / name
    if raw is not None:
        f.write_bytes(raw)
    else:
        f.write_text(text, encoding="utf-8")
    return str(f)


# --- K2 / SC-03: numbers written as text ---------------------------------
@pytest.mark.parametrize("t,p", [("4", "4.0"), ("12", "12.00"),
                                 ("4.50", "4.5"), ("-0", "0"),
                                 ("0.5", "0.50"), ("NaN", "nan")])
def test_k2_plain_decimals_and_nan_text_match(t, p):
    assert counts(t, p)[0] == 1


@pytest.mark.parametrize("t,p", [("007", "7"), ("02", "2"), ("1e3", "1000"),
                                 ("inf", "Infinity"), ("1_000", "1000"),
                                 (".5", "0.5"), ("+4", "4"), ("100", "100.")])
def test_k2_codes_exponents_and_specials_stay_text(t, p):
    assert counts(t, p)[1] == 1


def test_k2_normalize_returns_text():
    assert normalize("4") == "4"
    assert normalize("4.00") == "4"
    assert normalize(4.0) == "4"
    assert normalize("NaN") == "nan"


# --- SC-04: null-like tokens ---------------------------------------------
def test_sc04_null_and_none_are_empty():
    assert counts("Adeno", "null") == (0, 0, 1, 0, 0)
    assert counts("Adeno", "NONE") == (0, 0, 1, 0, 0)
    assert counts("None", "-") == (0, 0, 0, 0, 1)
    assert counts("-", "None") == (0, 0, 0, 0, 1)


def test_sc04_na_tokens_are_empty_unknown_stays():
    assert counts("Adeno", "N/A") == (0, 0, 1, 0, 0)
    assert counts("Adeno", "NA") == (0, 0, 1, 0, 0)
    assert counts("Adeno", "nA") == (0, 0, 1, 0, 0)
    assert counts("N/A", "n/a") == (0, 0, 0, 0, 1)
    assert counts("na", "B") == (0, 0, 0, 1, 0)
    assert counts("Adeno", "unknown") == (0, 1, 0, 0, 0)
    assert compare_binary("NA", "TRUE") is None
    assert compare_binary("TRUE", "N/A")["FN"] == 1


def test_sc04_na_through_reader(tmp_path):
    df = se.read_oncoshot(write(tmp_path,
        "Patient ID,H,Res: H\n1,Adeno,NA\n2,NA,n/a\n3,N/A,B\n"))
    m = se.evaluate(df, fields=["H"], types={"H": "scalar"}, n_boot=10,
                    random_state=1).metrics
    est = dict(zip(m.metric, m.estimate))
    assert (est["mis"], est["tn"], est["spu"]) == (1, 1, 1)


def test_sc04_dashes_mean_no_information():
    assert counts("\u2013", "-") == (0, 0, 0, 0, 1)
    assert counts("\u2014", "A") == (0, 0, 0, 1, 0)
    assert counts("A", "\u2013") == (0, 0, 1, 0, 0)


def test_sc04_null_tokens_can_be_switched_off():
    r = compare_set("Adeno", "None", scalar=True, null_tokens=())
    assert r["inc"] == 1


# --- SC-05: bracket lists -------------------------------------------------
@pytest.mark.parametrize("bad", ["['A', B]", "['A', 'B'].", "['A'; 'B']"])
def test_sc05_malformed_bracket_lists_raise(bad):
    with pytest.raises(ValueError, match="bracket"):
        parse_listish(bad)


def test_sc05_escapes_bare_and_blank_items():
    assert parse_listish("['O\\'Brien']") == ["O'Brien"]
    assert parse_listish("[\"O'Brien\", 'a\\\\b']") == ["O'Brien", "a\\b"]
    assert parse_listish("[A, b]") == ["A", "b"]
    assert parse_listish("['', 'A']") == ["A"]
    assert parse_listish("[]") == []


# --- SC-02: binary vocabulary ---------------------------------------------
def test_sc02_unrecognised_binary_prediction_raises():
    with pytest.raises(ValueError, match="Present"):
        compare_binary("TRUE", "Present")


def test_sc02_binary_spellings():
    assert compare_binary("TRUE", "1.0")["TP"] == 1
    assert compare_binary("FALSE", "0.0")["TN"] == 1
    assert compare_binary("TRUE", " Yes\u00a0")["TP"] == 1
    assert compare_binary("FALSE", "None", missing="wrong")["FP"] == 1
    assert compare_binary("TRUE", "\u2013", missing="wrong")["FN"] == 1


def test_sc02_unrecognised_binary_reference_raises():
    with pytest.raises(ValueError, match="maybe"):
        compare_binary("maybe", "TRUE")


# --- SC-07 / SC-16: blank and no-information references -------------------
def test_sc07_whitespace_only_is_unlabelled():
    assert is_unlabeled("   ")
    assert is_unlabeled("\u00a0")
    assert compare_set("  ", "X", scalar=True) is None


def test_sc16_dash_binary_reference_is_skipped():
    assert compare_binary("-", "TRUE") is None


# --- SC-08 / SC-11: case and whitespace -----------------------------------
def test_sc08_simple_lowercase_mapping():
    assert counts("\u00c9tape II", "\u00e9tape ii")[0] == 1
    assert counts("Stra\u00dfe", "STRASSE")[1] == 1   # no full case folding
    assert counts("\u039f\u0394\u039f\u03a3", "\u03bf\u03b4\u03bf\u03c3")[0] == 1
    assert normalize("\u0130") == "i"


def test_sc11_whitespace_trim_and_collapse():
    for t, p in [("Stage  IIA", "stage iia"), ("\u00a0Adeno", "Adeno"),
                 ("\vAdeno", "Adeno"), ("Adeno\u2003carcinoma", "Adeno carcinoma")]:
        assert counts(t, p)[0] == 1


# --- SC-06: declared field types ------------------------------------------
def _df(csv):
    return se.read_oncoshot(io.StringIO(csv))


def test_sc06_declared_binary_with_dash_label():
    df = _df("Patient ID,Met,Res: Met\n1,TRUE,TRUE\n2,FALSE,FALSE\n"
             "3,FALSE,FALSE\n4,TRUE,FALSE\n5,-,TRUE\n")
    m = se.evaluate(df, fields=["Met"], types={"Met": "binary"}, n_boot=20,
                    random_state=1).metrics
    est = dict(zip(m.metric, m.estimate))
    assert est["precision"] == 1.0 and est["recall"] == 0.5


def test_sc06_declared_binary_01():
    df = _df("Patient ID,Met,Res: Met\n1,1,1\n2,0,0\n3,0,0\n4,1,0\n")
    m = se.evaluate(df, fields=["Met"], types={"Met": "binary"}, n_boot=20,
                    random_state=1).metrics
    est = dict(zip(m.metric, m.estimate))
    assert est["recall"] == 0.5


def test_sc06_ambiguous_inference_raises():
    df = pd.DataFrame({"Met": ["TRUE", "FALSE", "-"],
                       "Res: Met": ["TRUE", "FALSE", "TRUE"]})
    with pytest.raises(ValueError, match="types"):
        se.evaluate(df, fields=["Met"], n_boot=10)
    df = pd.DataFrame({"Met": [1, 0, 1], "Res: Met": [1, 0, 0]})
    with pytest.raises(ValueError, match="types"):
        se.evaluate(df, fields=["Met"], n_boot=10)


def test_sc06_declared_scalar_with_bracket_cell():
    df = _df("Patient ID,H,Res: H\n1,A,['B']\n2,A,A\n")
    m = se.evaluate(df, fields=["H"], types={"H": "scalar"}, n_boot=20,
                    random_state=1).metrics
    est = dict(zip(m.metric, m.estimate))
    assert est["inc"] == 1 and est["mis"] == 0 and est["spu"] == 0


def test_sc06_scalar_with_two_items_raises():
    df = _df("Patient ID,H,Res: H\n1,A,\"['B', 'C']\"\n")
    with pytest.raises(ValueError, match="more than one item"):
        se.evaluate(df, fields=["H"], types={"H": "scalar"}, n_boot=10)


# --- SC-10: scalar TN, specificity and accuracy ---------------------------
def test_sc10_scalar_tn_specificity_accuracy():
    df = _df("Patient ID,H,Res: H\n1,A,A\n2,-,-\n3,-,B\n4,A,C\n5,-,-\n")
    m = se.evaluate(df, fields=["H"], types={"H": "scalar"}, n_boot=20,
                    random_state=1).metrics
    est = dict(zip(m.metric, m.estimate))
    assert est["tn"] == 2
    assert est["specificity"] == pytest.approx(2 / 3)
    assert est["accuracy"] == pytest.approx(3 / 5)


def test_k1_list_reports_no_specificity():
    df = _df("Patient ID,D,Res: D\n1,[],[]\n2,\"['A', 'B']\",\"['A', 'B']\"\n"
             "3,[],['C']\n")
    m = se.evaluate(df, fields=["D"], types={"D": "list"}, n_boot=20,
                    random_state=1).metrics
    assert "specificity" not in set(m.metric)
    assert "accuracy" not in set(m.metric)


# --- SC-14: F1 from counts ------------------------------------------------
def test_sc14_f1_zero_when_nothing_returned():
    m = set_metrics(cor=0, inc=0, mis=4, spu=0, par=0)
    assert math.isnan(m["precision"])
    assert m["recall"] == 0.0 and m["f1"] == 0.0
    b = binary_metrics(tp=0, fp=0, fn=3, tn=1)
    assert b["f1"] == 0.0


# --- SC-01 / SC-15: the reader ---------------------------------------------
def test_sc01_extra_field_raises(tmp_path):
    f = write(tmp_path, "Patient ID,D,Res: D\n1,['A'],['B']\n"
              "2,['B'],['B'],extra\n3,['C'],['D']\n")
    with pytest.raises(ValueError, match="line 3"):
        se.read_oncoshot(f)


def test_sc15_short_row_and_bad_quote_raise(tmp_path):
    with pytest.raises(ValueError, match="fields"):
        se.read_oncoshot(write(tmp_path, "Patient ID,D,Res: D\n1,A,A\n2,B\n"))
    with pytest.raises(ValueError, match="fields"):
        se.read_oncoshot(write(tmp_path, "Patient ID,D,Res: D\n1,\"A,A\n2,B,B\n"))


def test_sc15_duplicate_columns_utf8_empty_ids(tmp_path):
    with pytest.raises(ValueError, match="more than once"):
        se.read_oncoshot(write(tmp_path, "Patient ID,D,Res: D,Res: D\n1,A,A,B\n"))
    raw = b"Patient ID,H,Res: H\n1,Sj\xf6gren,x\n"
    with pytest.raises(ValueError, match="UTF-8"):
        se.read_oncoshot(write(tmp_path, None, raw=raw))
    with pytest.raises(ValueError, match="empty"):
        se.read_oncoshot(write(tmp_path, "Patient ID,H,Res: H\n1,A,A\n ,B,B\n"))


def test_reader_keeps_null_like_text(tmp_path):
    df = se.read_oncoshot(write(tmp_path,
        "\ufeffPatient ID,H,Res: H\n1,NaN,NA\n2,,None\n3,N/A,-\n"))
    assert list(df.index) == ["1", "2", "3"]
    assert df.loc["1", "H"] == "NaN" and df.loc["1", "Res: H"] == "NA"
    assert is_unlabeled(df.loc["2", "H"])
