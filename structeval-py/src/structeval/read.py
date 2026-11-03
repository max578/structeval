"""Strict reader for the Oncoshot wide CSV, mirroring ``read_oncoshot()`` in R.

Every cell is read as text. Null-like strings (``"NA"``, ``"N/A"``,
``"NaN"``, ``"None"``) stay text here, as in the R reader, and the scoring
rules then decide: ``"NA"``, ``"N/A"`` and ``"None"`` are empty values and
``"NaN"`` is a value. pandas' default ``na_values`` would instead turn them
into unlabelled cells. A blank cell is returned as ``""`` and is unlabelled
for :func:`structeval.evaluate`.
"""
from __future__ import annotations

import csv
import io
import os

import pandas as pd

from .compare import _trim


def read_oncoshot(file, id: str = "Patient ID") -> pd.DataFrame:
    """Read a wide CSV (path or text stream) into a DataFrame indexed by note.

    Raises ``ValueError`` when the file is not valid UTF-8, when a record has
    a different number of fields from the header (an unquoted comma would
    otherwise shift the columns), on an unterminated quote, when a column name
    appears more than once, or when a note identifier is empty or repeated.
    """
    if isinstance(file, (str, os.PathLike)):
        with open(file, "rb") as fh:
            raw = fh.read()
        try:
            text = raw.decode("utf-8-sig")
        except UnicodeDecodeError as err:
            line = raw[:err.start].count(b"\n") + 1
            raise ValueError(f"file is not valid UTF-8 on line {line}.") from err
    else:
        text = file.read()
        if text.startswith("\ufeff"):
            text = text[1:]

    reader = csv.reader(io.StringIO(text, newline=""), strict=True)
    rows = []
    try:
        header = next(reader)
        for rec in reader:
            if not rec:
                continue  # blank line, skipped as by read.csv()
            if len(rec) != len(header):
                raise ValueError(
                    f"Every record must have {len(header)} fields, as the "
                    f"header has; line {reader.line_num} has {len(rec)}. "
                    "Quote cells that contain a comma.")
            rows.append(rec)
    except csv.Error as err:
        raise ValueError(f"file could not be read cleanly (line "
                         f"{reader.line_num}): {err}; check the number of "
                         "fields and the quotes.") from err

    repeated = sorted({c for c in header if header.count(c) > 1})
    if repeated:
        raise ValueError(f"Columns appear more than once: {repeated}.")
    if id not in header:
        raise ValueError(f"Column '{id}' not found.")
    df = pd.DataFrame(rows, columns=header, dtype=object)
    notes = df[id].map(_trim)
    if (notes == "").any():
        raise ValueError(f"Note identifiers in column '{id}' are empty in "
                         f"data row(s) {list(notes.index[notes == ''] + 1)}.")
    if notes.duplicated().any():
        raise ValueError(f"Note identifiers in column '{id}' are not unique.")
    df.index = notes.astype(str)
    return df.drop(columns=[id])
