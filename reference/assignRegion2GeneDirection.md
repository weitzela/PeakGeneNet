# The results from all correlations are used to help assign direction using correlations that have a significant nominal pvalue

Directional assignments are categorized into 3 levels: 1) Most
confident: assign direction using FDR significant links and they are all
in agreement 2) Directional disagreement existed, but decision was made
by taking the most significant FDR 3) Directional disagreement existed,
but decision was made by taking most significant nominal pvalue that was
included in comparison 4) These are looser directions associated with
the gene - neither distal or promoter region is nominally correlated
with the gene, so the direction is based on how that peak is associated
with distal peaks that are significantly linked to the gene with high
confidence. The iteration of identifying a significant correlation
between the peak and another distal peak is divided by 1000 and then
added to 4. This shows that the peak:gene relationship exists at
confidence level \#4, and notes how many distal genes it had to try to
correlate with before finding a significant (pval \< 0.05) connection.

## Usage

``` r
assignRegion2GeneDirection(sig_cor_res, all_cor_res, count_mat)
```

## Arguments

- sig_cor_res:

  dataframe with the significant results to annotate link direction

- all_cor_res:

  all resulting correlations

- count_mat:

  matrix of counts supplied to correlation analysis

## Value

dataframe with peak-gene link annotation directions with the degree of
confidence
