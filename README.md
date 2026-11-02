# structeval

Structure-aware evaluation of fields that a model extracted from clinical notes.

Every evaluation metric here is a ratio of note-level sums, and the interval is chosen by that structure. Binary and scalar fields are note-level proportions and receive Wilson score intervals. List fields pool several items per note and receive a note-level cluster bootstrap with bias-corrected and accelerated percentiles, or a cluster-robust delta interval, with a stated rule at the boundary. On top of the intervals sit three tests and a set of sample-size planners: a certification test against a threshold, a paired comparison of two systems on the same notes, and a test of whether a system omits more than it fabricates.

A methods paper describing these procedures is in preparation. The worked example in the package vignette uses the aggregate counts published in Moldovan et al. (2026), medRxiv 2026.05.18.26353541, <https://doi.org/10.64898/2026.05.18.26353541>.

## What is in this repository

| Folder | Contents |
|---|---|
| `structeval-r/` | The R package, version 0.1.0. Source, tests, vignette. |
| `structeval-py/` | The Python reference implementation used to cross-check the R results, with its tests. |
| `sim/` | The size-and-power simulation study: code, seeds, logs and result tables. |
| `release/` | The built source package `structeval_0.1.0.tar.gz` with its checksum. |
| `samples.csv` | A 14-row synthetic table in the extraction layout, one row per outcome type. The Python tests and the R README example read it. |

## Install the R package

From this repository:

```r
install.packages("remotes")
remotes::install_github("max578/structeval", subdir = "structeval-r")
```

From the built source package:

```r
install.packages("release/structeval_0.1.0.tar.gz", repos = NULL, type = "source")
```

The package depends on base R only. The vignette and the tests use knitr, rmarkdown, testthat, ggplot2 and jsonlite.

## Install the Python reference implementation

```sh
pip install ./structeval-py
python -m pytest structeval-py/tests
```

It needs Python 3.11 or later, pandas and numpy.

## Cross-check between the two implementations

`structeval-r/tests/fixtures/make_crosscheck_python.py` runs the Python implementation on the shared fixture `crosscheck.csv` and writes `crosscheck_python.json`. The R test `test-crosscheck.R` recomputes the same quantities and compares them to that file. The fixture and the Python output ship with the package, so the comparison runs as part of the R tests without Python installed.

## Reproduce the simulation study

From the repository root:

```sh
Rscript sim/run_sim.R 2000
Rscript sim/figures/fig2_sizepower.R
```

The full study takes about 14 minutes on 5 cores. Finished cells are cached in `sim/results/cells_reps2000/` and skipped on a re-run, so the shipped cache makes the second command run from the tables alone. Seeds, grid design and the result files are described in `sim/README.md`.

## Citation

See `CITATION.cff`. Until the methods paper appears, cite the software.

## Licence

MIT. See `LICENSE`.
