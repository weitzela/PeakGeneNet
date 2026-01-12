# Correlation of two matrices for more efficient calculations

Correlation of two matrices for more efficient calculations

## Usage

``` r
matrixCorrelation(
  mat1,
  mat2,
  type = c("pearson", "spearman"),
  rank_data = FALSE,
  return_pvalue = TRUE
)
```

## Arguments

- mat1:

  first matrix of data to correlate

- mat2:

  must have same dimensions as mat1

- type:

  method of correlation

- rank_data:

  logical

- return_pvalue:

  default=TRUE

## Value

A dataframe with five columns containing the pairs of data that were
correlated, the correlation coefficient, p-value, and number of samples
included in the correlation.

## Examples

``` r
# example code
```
