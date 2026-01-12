# Wrapper function for corWithDistalPeak to iterate through the correlations

Wrapper function for corWithDistalPeak to iterate through the
correlations

## Usage

``` r
corWithDistalPeak_wrapper(
  full_cor_df,
  count_mat,
  v = TRUE,
  .iter_batch_cutoff = 5
)
```

## Arguments

- full_cor_df:

  dataframe of full correlation results

- count_mat:

  matrix of counts supplied to correlation analysis

- v:

  logical. print informative messages

## Value

dataframe with iterated correlation results until a significant
correlation was identified with other distal peaks
