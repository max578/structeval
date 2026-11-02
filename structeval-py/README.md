# structeval, Python reference implementation

This is the Python counterpart of the R package in `../structeval-r/`. It exists to cross-check the R results and to serve readers who work in Python. The R package is the reference for the methods paper.

`evaluate()` reads the extraction table convention of the Oncoshot framework: one label column `F` and one prediction column `Res: F` per field, with `-` meaning no information and a blank meaning unlabelled. It returns a tidy report with one row per field, confidence level and metric, each carrying the observed point estimate and an interval computed by the method that fits the field: Wilson for binomial proportions, a bias-corrected and accelerated bootstrap for ratios of sums. Counts are reported as counts, without an interval.

```python
import pandas as pd
import structeval as se

df = pd.read_csv("data.csv", index_col="Patient ID")
report = se.evaluate(df)        # fields inferred from the "Res: " columns
print(report.metrics)           # field x confidence x metric, with intervals
print(report.pivot())           # wide view of the point estimates
report.to_csv("results/")       # metrics.csv and cases.csv
```

## Install and test

From the repository root, so that the tests find `samples.csv`:

```sh
pip install ./structeval-py
python -m pytest structeval-py/tests
```

Python 3.11 or later, pandas 2 and numpy are required.

## Design choices

| Concern | Choice |
|---|---|
| Point estimate | The observed statistic, never a bootstrap mean |
| Binary proportions | Wilson score interval |
| F-scores and set metrics | Bias-corrected and accelerated bootstrap |
| Counts | Reported as counts, no interval |
| Undefined bootstrap replicates | Counted and reported |
| Missing binary prediction | Explicit `missing=` policy |
| Repeated measures | Optional cluster bootstrap through `cluster=` |
| Partial labels | Handled per field |

Every interval carries the `method` that produced it (`wilson`, `bca` or `count`).

## Cross-check with R

`../structeval-r/tests/fixtures/make_crosscheck_python.py` runs this implementation on the shared fixture and writes the values that the R test suite compares against, to a tolerance of 1e-12.

## Licence

MIT. See `../LICENSE`.
