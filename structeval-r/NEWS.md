# structeval 0.1.1

## Scoring and reading

* List fields no longer report specificity, and `metrics()`,
  `test_threshold()` and the `power_*()` functions stop when it is
  requested, as for accuracy. The base paper reports neither for list
  fields. The `tn` count of a list field is kept.
* `score()` compares plain decimals as numbers, so `"4"`, `"4.0"` and
  `"4.00"` match. Codes such as `"007"`, exponent forms such as `"1e3"`,
  and `"NaN"` stay text.
* `score()` treats `"null"`, `"none"`, `"NA"` and `"N/A"` (any case) as
  empty values, the base paper's "empty or null", and en and em dashes as
  `"-"`. `"unknown"` stays a value. New argument `null_tokens`;
  `character(0)` keeps them as values.
* `score()` trims and collapses every Unicode space, including no-break
  and vertical spaces, and stops when non-ASCII text would be case-folded
  outside a UTF-8 locale.
* A binary prediction outside the TRUE/FALSE vocabulary now stops with the
  values listed, as a reference already did. `"1.0"` and `"0.0"` are read
  as flags.
* A blank or white-space reference is unlabelled for every field type.
  Duplicated note and field rows stop `score()`.
* `score()` returns a `coverage` attribute with the rows skipped as
  unlabelled, skipped as binary references with no information, with a
  missing prediction and dropped under `missing = "skip"`, and reports
  skipped rows by a message. `metrics()` carries the attribute over.
* `read_oncoshot()` stops on a record with the wrong number of fields
  (an unquoted comma shifted every column), on invalid UTF-8, on a
  parser warning, on repeated columns and on empty note identifiers.
  Bracketed lists must be all quoted or all bare; escaped quotes are
  unescaped.

## Intervals and planning

* `metrics()` and `test_threshold()`: the `cluster_bca` and
  `cluster_delta` lower bound is now capped at the Wilson lower bound of a
  perfect record on the same notes (and the upper bound floored at that of
  an empty record), so a record with errors never certifies more than a
  perfect one. In 0.1.0, 200 of 200 items in 50 notes failed to certify
  0.95 while 199 of 200 passed. Now neither does.
* `n_threshold()` documents that its design effect assumes equal note
  sizes and can under-plan when sizes vary.
* `n_paired()` now returns the first n, from Connor's formula upward, at
  which the exact McNemar test of `test_paired()` reaches the requested
  power. Connor's formula alone gave about 76% to 78% power at a target of
  80%. New columns `n_connor` and `exact_power`.
* `test_asymmetry(method = "cluster_t")` on a field with one note returns
  `NA` with a warning instead of an error.

# structeval 0.1.0

* First release. `read_oncoshot()` and `score()` count correct, incorrect,
  missed, spurious and true-negative outcomes per note and field.
  `metrics()` reports precision, recall, F1, specificity and accuracy with
  Wilson intervals for note-level proportions and a note-level cluster
  bootstrap (BCa) for list fields, with a stated rule at the boundary.
* Tests: `test_threshold()`, `test_paired()` and `test_asymmetry()`.
* Planning: `n_threshold()`, `n_paired()`, `min_certifiable()`,
  `simulate_counts()` and the `power_*()` functions.
* `plot_errors()` draws the error decomposition of each field.
