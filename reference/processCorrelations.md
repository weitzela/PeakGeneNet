# Process correlation results

Genes with a significant correlation with either a promoter or distal
peak are retained and focused on further.

## Usage

``` r
processCorrelations(
  count_mat,
  cor_res,
  correlation_pairs,
  p2g_info,
  gene_inclusion_thresh = 0.05,
  BH_thresh = 0.05
)
```
