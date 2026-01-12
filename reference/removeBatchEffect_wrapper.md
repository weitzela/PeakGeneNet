# Wrapper: Remove Batch Effect This supplies the experimental model to remove technical variables, while protecting for the experimental variables.

Wrapper: Remove Batch Effect This supplies the experimental model to
remove technical variables, while protecting for the experimental
variables.

## Usage

``` r
removeBatchEffect_wrapper(
  counts,
  df,
  var_adj,
  var_protect = c("line", "sex", "grp", "train.sed")
)
```
