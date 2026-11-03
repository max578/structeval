# structeval

structeval scores fields that a model extracted from clinical notes against
reference labels. It reports precision, recall, F1, specificity and accuracy
with intervals that match each field's structure. Binary and scalar fields
are note-level proportions and get Wilson intervals. List fields pool
several items per note and get a note-level cluster bootstrap. The package
also tests a threshold, compares two systems, tests whether a system misses
more than it adds, and plans sample sizes.

## Install

From GitHub:

```r
install.packages("remotes")
remotes::install_github("max578/structeval", subdir = "structeval-r")
```

From the built source package:

```r
install.packages("structeval_0.1.1.tar.gz", repos = NULL, type = "source")
```

## Example

```r
library(structeval)
long <- read_oncoshot("samples.csv", types = c(
  "Has metastasis" = "binary", "First Primary Diagnosis" = "scalar",
  "Treatment Drugs" = "list"
))
counts <- score(long)
metrics(counts, seed = 1)
test_threshold(counts, "Treatment Drugs", "f1", theta0 = 0.90, seed = 1)
test_asymmetry(counts, "Treatment Drugs", seed = 1)
n_threshold(theta0 = 0.90, theta1 = 0.95, m = 4, rho = 0.2)
```

## Tests and planners

| Function | Question | Methods |
|---|---|---|
| `test_threshold()` | Does a metric exceed a threshold? | Wilson, cluster BCa, cluster delta |
| `test_paired()` | Do two systems differ on the same notes? | cluster delta, cluster BCa, sign flip, McNemar |
| `test_asymmetry()` | Are missed and spurious items balanced? | sign flip, cluster t, binomial |
| `n_threshold()` | How many notes certify a threshold? | closed form with design effect |
| `n_paired()` | How many notes detect a paired difference? | exact McNemar, searched up from Connor (1987) |

`power_threshold()`, `power_paired()` and `power_asymmetry()` give Monte
Carlo power from `simulate_counts()`. `min_certifiable()` gives the smallest
count that certifies a threshold at a given sample size.

## Paper

The deployment example in the vignette uses the counts published in
Moldovan et al. (2026), medRxiv 2026.05.18.26353541,
<https://doi.org/10.64898/2026.05.18.26353541>.

## Licence

MIT.
