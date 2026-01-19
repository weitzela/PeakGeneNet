# Format input

This function pastes the modality at the end of each feature ID and then
combines the matrices into one matrix object. The matrices must have
common sample IDs as the rownames to join by.

## Usage

``` r
formatMatrixForCorrelation(gene_counts, peak_counts)
```

## Arguments

- gene_counts:

  count matrix where rownames are sample IDs and colnames are feature ID
