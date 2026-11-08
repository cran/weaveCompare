# weaveCompare

Periodic binary drawdown comparison for R. Rows are picks, columns are ends;
1 means warp above weft. Normalize software viewing conventions before use.

```r
install.packages("weaveCompare_0.1.0.tar.gz", repos = NULL, type = "source")
library(weaveCompare)
a <- matrix(c(1,0,0,1), 2)
b <- shift_drawdown(a, 1)
compare_drawdowns(a, b)
duplicate_drawdowns(list(original=a, shifted=b))
```

With weaveDraft installed, pass a weave_draft object directly or its drawdown
matrix. WIF import remains the responsibility of weaveDraft. This package
has no dependency on an unpublished package, MATLAB or WeavePoint.

Cyclic phase shifts are always allowed. Reflection, transpose and inversion
are opt-in. Transpose alone exchanges axes without exchanging bits: users
requiring physical warp/weft relabelling must consider inversion as well.
The combined transpose/invert option searches both independently.
Only complete rectangular repeats are reduced. No image registration,
thread colour comparison, yarn geometry or physical-property prediction is
performed. This is not a cultural attribution or authentication tool.

The exhaustive algorithms target modest repeat matrices, not full fabric scans.
Work limits guard search; increasing them may cause large memory/runtime use.
Collection budgets are per item/pair, not an overall runtime limit.

All examples are synthetic. No institutional or community data are distributed.
The methods implement established periodic matrix operations; novelty and
CRAN acceptance are not asserted. Author/licensing follow the prior package.
