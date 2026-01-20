# workflow

``` r
# install.packages("pak")
# pak::pak("weitzela/PeakGeneNet")

library(tidyverse)
library(PeakGeneNet)
```

## 1. Peak-Gene Link Preparation

Use these objects to create peak-gene pairs for inclusion in the
correlation analysis. The `p2g_info` provides information for individual
peaks per gene they will be correlated with. The `correlation_pairs`
object is used for reference for all of the correlations that will be
run.

``` r
p2g_ls = createPeak2GeneObjects(
  rownames(gene_counts), peak_counts |> map(rownames), 
  biomaRt::useEnsembl(biomart = "genes", 
                      dataset = "rnorvegicus_gene_ensembl", version = 109), 
  "rn7"
)
p2g_ls$p2g_info[1:5,] |> kableExtra::kable()
```

| unique_id                      | region_id              | modality | ensembl_gene_id    | dist_to_tss | promoter_peak |
|:-------------------------------|:-----------------------|:---------|:-------------------|------------:|:--------------|
| chr3:78634005:78634640_ATACSeq | chr3:78634005:78634640 | ATACSeq  | ENSRNOG00000000008 |     -977078 | FALSE         |
| chr3:78652401:78653178_H3K4me1 | chr3:78652401:78653178 | H3K4me1  | ENSRNOG00000000008 |     -958540 | FALSE         |
| chr3:78652614:78653063_ATACSeq | chr3:78652614:78653063 | ATACSeq  | ENSRNOG00000000008 |     -958655 | FALSE         |
| chr3:78683169:78683832_H3K27ac | chr3:78683169:78683832 | H3K27ac  | ENSRNOG00000000008 |     -927886 | FALSE         |
| chr3:78683236:78683697_H3K4me1 | chr3:78683236:78683697 | H3K4me1  | ENSRNOG00000000008 |     -928021 | FALSE         |

``` r
p2g_ls$correlation_pairs[1:5,] |> kableExtra::kable()
```

| ensembl_gene_id    | regulatory_element             | target_id                      | link_label                   | chr  | modality_pair   |
|:-------------------|:-------------------------------|:-------------------------------|:-----------------------------|:-----|:----------------|
| ENSRNOG00000000008 | chr3:78634005:78634640_ATACSeq | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 | ATACSeq-H3K4me3 |
| ENSRNOG00000000008 | chr3:78652401:78653178_H3K4me1 | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 | H3K4me1-H3K4me3 |
| ENSRNOG00000000008 | chr3:78652614:78653063_ATACSeq | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 | ATACSeq-H3K4me3 |
| ENSRNOG00000000008 | chr3:78683169:78683832_H3K27ac | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 | H3K27ac-H3K4me3 |
| ENSRNOG00000000008 | chr3:78683236:78683697_H3K4me1 | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 | H3K4me1-H3K4me3 |

See how many relationships will be tested within each category.

``` r
p2g_ls$correlation_pairs |> 
  select(-ensembl_gene_id) |> 
  distinct() |> 
  group_by(link_label, modality_pair) |> 
  tally(name = "No. of Correlations") |> 
  janitor::adorn_totals(where = c("row")) |> 
  mutate(across(where(is.numeric), scales::comma)) |> 
  kableExtra::kable()
```

[TABLE]

## 2. Data Adjustment

If your experimental design is simple, you have a lot of samples, and
your data does not show signs of variation, then you can skip this step
and move onto the correlation analysis.

Correlation-based analyses are sensitive to unwanted variation. In
multi-factor experimental designs, failing to account for technical
effects or biological variables can reduce power, inflate variance, or
introduce type I errors. To prepare the data for use in the PeakGeneNet
pipeline, we recommend an adjustment strategy targeting both technical
and experimental correction while preserving the design and contrast of
interest. This approach allows users to ask multiple questions from the
same dataset while ensuring that the signal relevant to a contrast of
interest is preserved. Different contrasts can be explored by re-running
this adjustment step with different protected variables.

Data provided at this stage should **not** be raw counts. Instead, input
matrices should be transofrmed (e.g., VST, log, or inverse-rank
normalization) to improve comparability across samples. Adjustment is
performed in two steps, following standard best practices used in
RNA-Seq and epigenomic analyses that preserve biological signal while
removing unwanted variation:

1.  Removal of technical and latent sources of variation, while
    explicitly preserving the experimental design.
    - Latent sources of variation in this example are variables prefixed
      by “W\_” and included because they capture technical effects.  
2.  Regression of additional experimental covariates, while protecting
    one or more contrasts of interest.
    - Protected variables should correspond directly to the hypothesis
      being tested.

The example data files in this vignette contain a `sample_info`
attribute where libraries are ordered in the same order as what is
contained in the count matrix. The objects are set up in this manner to
easily apply the same adjustment strategy to all data. This example
study has four experimental variables of rat lines, sex, acute exercise
groups, and a training group.

``` r
# this is a way to apply the same adjustment strategy across all modalities
experimental_vars = c("line", "sex", "grp", "train.sed") 
adj_counts = c(list(RNASeq = gene_counts), peak_counts) |> 
  map(function(.count_mat) {
    pheno_df = attr(.count_mat, "samp")
    # step 1: remove technical variation
    adj_tmp = limma::removeBatchEffect(
      .count_mat, 
      batch = if ("prep_date" %in% colnames(pheno_df)) pheno_df[["prep_date"]] else pheno_df[[""]],
      batch2 = if ("flowcell" %in% colnames(pheno_df)) pheno_df[["flowcell"]] else pheno_df[[""]],
      covariates = pheno_df |> select(starts_with("W_")),
      design = model.matrix(reformulate(experimental_vars), data = pheno_df)
    )
    # step 2: regress out additional experimental covariates, while retaining effects attributed to the experimental group of interest
    adj_counts = adjustCovariateMatrix(
      counts = t(adj_tmp),
      covariate_df = pheno_df |> select(all_of(experimental_vars)), 
      vars_to_protect = "grp" 
    )
    # replace library IDs with common sample IDs that can be matched across matrices 
    rownames(adj_counts) = attr(adj_counts, "samp")[rownames(adj_counts), "samp_id"]
    return(adj_counts)
  })
```

## 3. Correlations

Then, carry out and process the correlations using the following
functions. The count matrices used as the input for the correlation
analysis should have samples in rows and features in columns.

``` r
count_mat = formatMatrixForCorrelation(adj_counts$RNASeq, adj_counts |> discard_at("RNASeq"))
# count_mat[1:5,1:5] |> kableExtra::kable()

peak_counts_t = map(peak_counts, function(.x) {
  colnames(.x) = attr(.x, "samp")[colnames(.x), "samp_id"]
  return(t(.x))
}); count_mat = formatMatrixForCorrelation(t(gene_counts), peak_counts_t)

full_cor_res = correlateByChromosome(count_mat, p2g_ls$correlation_pairs)
#> Starting chr1
#> 0.02 sec elapsed
#> Starting chr2
#> 0.012 sec elapsed
#> Starting chr3
#> 0.007 sec elapsed
#> Joining with `by = join_by(regulatory_element, target_id, link_label, chr,
#> modality_pair, r, n, P)`
full_cor_res |> head(5) |> kableExtra::kable()
```

| ensembl_gene_id    | regulatory_element             | target_id                      | link_label                   | chr  | modality_pair   |          r |   n |         P |        BH |
|:-------------------|:-------------------------------|:-------------------------------|:-----------------------------|:-----|:----------------|-----------:|----:|----------:|----------:|
| ENSRNOG00000000008 | chr3:78634005:78634640_ATACSeq | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 | ATACSeq-H3K4me3 |  0.1656186 |  56 | 0.2225149 | 0.8032780 |
| ENSRNOG00000000008 | chr3:78652401:78653178_H3K4me1 | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 | H3K4me1-H3K4me3 |  0.0949883 |  59 | 0.4742207 | 0.9788157 |
| ENSRNOG00000000008 | chr3:78652614:78653063_ATACSeq | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 | ATACSeq-H3K4me3 | -0.0216678 |  56 | 0.8740562 | 0.9817481 |
| ENSRNOG00000000008 | chr3:78683169:78683832_H3K27ac | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 | H3K27ac-H3K4me3 |  0.0444185 |  59 | 0.7383421 | 0.9916836 |
| ENSRNOG00000000008 | chr3:78683236:78683697_H3K4me1 | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 | H3K4me1-H3K4me3 |  0.1593074 |  59 | 0.2281274 | 0.9237553 |

``` r

# this will not run because of the limited genes and peaks included in the example data
# but this is a way to process and interpret the results 
# processed_cor = processCorrelations(count_mat, full_cor_res, p2g_ls$correlation_pairs, p2g_ls$p2g_info)
# processed_cor |> head(5) |> kableExtra::kable()
```

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] PeakGeneNet_0.1.0 lubridate_1.9.4   forcats_1.0.1     stringr_1.6.0    
#>  [5] dplyr_1.1.4       purrr_1.2.1       readr_2.1.6       tidyr_1.3.2      
#>  [9] tibble_3.3.1      ggplot2_4.0.1     tidyverse_2.0.0  
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1     viridisLite_0.4.2    farver_2.1.2        
#>  [4] blob_1.3.0           filelock_1.0.3       Biostrings_2.78.0   
#>  [7] S7_0.2.1             fastmap_1.2.0        BiocFileCache_3.0.0 
#> [10] janitor_2.2.1        digest_0.6.39        timechange_0.3.0    
#> [13] lifecycle_1.0.5      statmod_1.5.1        KEGGREST_1.50.0     
#> [16] RSQLite_2.4.5        magrittr_2.0.4       compiler_4.5.2      
#> [19] rlang_1.1.7          sass_0.4.10          progress_1.2.3      
#> [22] tools_4.5.2          yaml_2.3.12          knitr_1.51          
#> [25] prettyunits_1.2.0    bit_4.6.0            curl_7.0.0          
#> [28] xml2_1.5.2           RColorBrewer_1.1-3   withr_3.0.2         
#> [31] BiocGenerics_0.56.0  desc_1.4.3           grid_4.5.2          
#> [34] stats4_4.5.2         scales_1.4.0         gtools_3.9.5        
#> [37] biomaRt_2.66.0       cli_3.6.5            rmarkdown_2.30      
#> [40] crayon_1.5.3         ragg_1.5.0           generics_0.1.4      
#> [43] rstudioapi_0.18.0    httr_1.4.7           tzdb_0.5.0          
#> [46] DBI_1.2.3            cachem_1.1.0         AnnotationDbi_1.72.0
#> [49] XVector_0.50.0       matrixStats_1.5.0    vctrs_0.7.0         
#> [52] jsonlite_2.0.0       IRanges_2.44.0       hms_1.1.4           
#> [55] S4Vectors_0.48.0     bit64_4.6.0-1        systemfonts_1.3.1   
#> [58] limma_3.66.0         jquerylib_0.1.4      glue_1.8.0          
#> [61] pkgdown_2.2.0        stringi_1.8.7        gtable_0.3.6        
#> [64] GenomeInfoDb_1.46.2  GenomicRanges_1.62.1 UCSC.utils_1.6.1    
#> [67] pillar_1.11.1        rappdirs_0.3.4       htmltools_0.5.9     
#> [70] Seqinfo_1.0.0        R6_2.6.1             dbplyr_2.5.1        
#> [73] httr2_1.2.2          textshaping_1.0.4    evaluate_1.0.5      
#> [76] kableExtra_1.4.0     Biobase_2.70.0       png_0.1-8           
#> [79] tictoc_1.2.1         snakecase_0.11.1     memoise_2.0.1       
#> [82] bslib_0.9.0          svglite_2.2.2        xfun_0.56           
#> [85] fs_1.6.6             pkgconfig_2.0.3
```
