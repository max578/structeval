"""Run the Python reference implementation on a cross-check fixture.

Usage (from the project root), for each fixture:
    python3 -I structeval-r/tests/fixtures/make_crosscheck_python.py \
        src structeval-r/tests/fixtures/crosscheck.csv \
        structeval-r/tests/fixtures/crosscheck_python.json
    python3 -I structeval-r/tests/fixtures/make_crosscheck_python.py \
        src structeval-r/tests/fixtures/crosscheck_adversarial.csv \
        structeval-r/tests/fixtures/crosscheck_adversarial_python.json

The file is read with structeval.read_oncoshot() and the field types are
declared, as in R. Writes per-note counts (including TN for scalar and list
fields), point estimates and Wilson bounds for every field. Wilson bounds
come from evaluate() for binary fields and from structeval.uncertainty.wilson()
on the pooled (k, n) for scalar precision, recall, specificity and accuracy
(note-level proportions) and for list precision, recall and F1 (the naive
item-level interval).
"""
import json
import math
import os
import sys

src, csv_path, json_path = sys.argv[1:4]
sys.path.insert(0, src)

import pandas as pd  # noqa: E402

from structeval import evaluate, read_oncoshot  # noqa: E402
from structeval.uncertainty import wilson, z_for  # noqa: E402

FIXTURES = {
    "crosscheck.csv": {
        "Has metastasis": "binary",
        "First Primary Diagnosis": "scalar",
        "First Primary Histology": "scalar",
        "Treatment Drugs": "list",
        "Radiotherapy Sites": "list",
    },
    "crosscheck_adversarial.csv": {
        "Met": "binary", "Hist": "scalar", "Size": "scalar", "Drugs": "list",
    },
}
FIELDS = FIXTURES[os.path.basename(csv_path)]


def num(x):
    x = float(x)
    return None if math.isnan(x) else x


df = read_oncoshot(csv_path)
rep = evaluate(df, fields=list(FIELDS), types=FIELDS, n_boot=50,
               random_state=1)
z = z_for(0.95)
out = {"python": sys.version.split()[0], "pandas": pd.__version__,
       "fields": {}}
m = rep.metrics[rep.metrics.confidence == "Overall"]

for field, kind in FIELDS.items():
    rows = m[m.field == field]
    names = (["TP", "FP", "FN", "TN"] if kind == "binary"
             else ["cor", "inc", "mis", "spu", "tn"])
    per_note = {}
    for note in rep.cases.index:
        vals = [rep.cases.loc[note, f"{field} [{c}]"] for c in names]
        if any(isinstance(v, float) and math.isnan(v) for v in vals):
            continue
        per_note[str(note)] = dict(zip(names, (int(v) for v in vals)))
    tot = {c: sum(r[c] for r in per_note.values()) for c in names}
    est = {r.metric: num(r.estimate) for r in rows.itertuples()
           if r.method != "count"}
    wil = {}
    if kind == "binary":
        for r in rows.itertuples():
            if r.method == "wilson":
                wil[r.metric] = [num(r.ci_low), num(r.ci_high)]
    else:
        c, i, mi, s = tot["cor"], tot["inc"], tot["mis"], tot["spu"]
        tn = tot["tn"]
        kn = {"precision": (c, c + s + i), "recall": (c, c + mi + i)}
        if kind == "scalar":
            kn["specificity"] = (tn, tn + s)
            kn["accuracy"] = (c + tn, c + tn + i + mi + s)
        if kind == "list":
            kn["f1"] = (2 * c, 2 * c + mi + s + 2 * i)
        for name, (k, n) in kn.items():
            lo, hi = wilson(k, n, z)
            wil[name] = [num(lo), num(hi)]
    out["fields"][field] = {"type": kind, "per_note": per_note,
                            "totals": tot, "estimates": est, "wilson": wil}

with open(json_path, "w") as fh:
    json.dump(out, fh, indent=1, sort_keys=True)
print("wrote", json_path, "fields:", len(out["fields"]))
