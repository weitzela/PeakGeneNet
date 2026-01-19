# A wrapper function to carry out correlations in smaller chunks to help save on memory

A wrapper function to carry out correlations in smaller chunks to help
save on memory

## Usage

``` r
correlateByChromosome(
  count_mat,
  correlation_pairs,
  grp_contrast = NULL,
  rds_fn = NULL
)
```

## Arguments

- count_mat:

  numeric matrix with features as columns and samples as rows

- correlation_pairs:

  dataframe supplying all pairs of features

- grp_contrast:

  a way to filter to include samples from only certain groups

- rds_fn:

  filename to save the output to for future access

## Value

a dataframe of full results
