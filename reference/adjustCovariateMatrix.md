# Use of matrix functions to adjust data via linear model

This is a wrapper function that adjusts for multiple additional
experimental variables while protecting the contrast(s) of interest.

## Usage

``` r
adjustCovariateMatrix(
  counts,
  covariate_df,
  vars_to_protect,
  return_coefficients = FALSE
)
```

## Arguments

- counts:

  matrix of transformed counts (vst, log, inverse rank normalized,
  etc.). Should have sample IDs as rownames, feature ID as colnames

- covariate_df:

  dataframe with all covariates, sample names should be rownames. should
  be in the same order as the input matrix

- vars_to_protect:

  character vector of column names that should be protected. the
  coefficients calculated for these variables will not be subtracted
  from the returned values

- return_coefficients:

  logcal. default=FALSE

## Value

matrix of adjusted counts
