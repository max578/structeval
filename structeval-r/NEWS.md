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
