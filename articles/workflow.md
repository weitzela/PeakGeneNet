# workflow

``` r
# install.packages("pak")
# pak::pak("weitzela/PeakGeneNet")

library(tidyverse)
library(PeakGeneNet)
```

## 1. Data Preparation

Required objects:

- Gene counts: this is a matrix where the row names should be ensembl
  IDs

- Peak counts: a list object with one or more matrices. Matrices should
  be in same format as gene counts, with row names as genomic region
  where the chromosome names, start, and end locations are separaated by
  a punctuation mark (e.g., chr1:123:456 or chr1:123-456).

Use these objects to create peak-gene pairs for inclusion in the
correlation analysis. The `p2g_info` provides information for individual
peaks per gene they will be correlated with. The `correlation_pairs`
object is used for reference for all of the correlations that will be
run.

``` r
p2g_ls = createPeak2GeneObjects(
  gene_counts, peak_counts, 
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

| ensembl_gene_id    | regulatory_element             | target_id                      | link_label                   | chr  |
|:-------------------|:-------------------------------|:-------------------------------|:-----------------------------|:-----|
| ENSRNOG00000000008 | chr3:78634005:78634640_ATACSeq | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 |
| ENSRNOG00000000008 | chr3:78652401:78653178_H3K4me1 | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 |
| ENSRNOG00000000008 | chr3:78652614:78653063_ATACSeq | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 |
| ENSRNOG00000000008 | chr3:78683169:78683832_H3K27ac | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 |
| ENSRNOG00000000008 | chr3:78683236:78683697_H3K4me1 | chr3:79610922:79613046_H3K4me3 | distal_peak_to_promoter_peak | chr3 |

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
#>  [1] KEGGREST_1.50.0      gtable_0.3.6         httr2_1.2.2         
#>  [4] xfun_0.56            bslib_0.9.0          Biobase_2.70.0      
#>  [7] tzdb_0.5.0           vctrs_0.7.0          tools_4.5.2         
#> [10] generics_0.1.4       curl_7.0.0           stats4_4.5.2        
#> [13] AnnotationDbi_1.72.0 RSQLite_2.4.5        blob_1.3.0          
#> [16] pkgconfig_2.0.3      dbplyr_2.5.1         RColorBrewer_1.1-3  
#> [19] S7_0.2.1             desc_1.4.3           S4Vectors_0.48.0    
#> [22] lifecycle_1.0.5      compiler_4.5.2       farver_2.1.2        
#> [25] Biostrings_2.78.0    progress_1.2.3       tictoc_1.2.1        
#> [28] textshaping_1.0.4    Seqinfo_1.0.0        GenomeInfoDb_1.46.2 
#> [31] htmltools_0.5.9      sass_0.4.10          yaml_2.3.12         
#> [34] crayon_1.5.3         pillar_1.11.1        pkgdown_2.2.0       
#> [37] jquerylib_0.1.4      cachem_1.1.0         gtools_3.9.5        
#> [40] tidyselect_1.2.1     digest_0.6.39        stringi_1.8.7       
#> [43] biomaRt_2.66.0       fastmap_1.2.0        grid_4.5.2          
#> [46] cli_3.6.5            magrittr_2.0.4       withr_3.0.2         
#> [49] prettyunits_1.2.0    filelock_1.0.3       rappdirs_0.3.4      
#> [52] scales_1.4.0         UCSC.utils_1.6.1     bit64_4.6.0-1       
#> [55] timechange_0.3.0     XVector_0.50.0       rmarkdown_2.30      
#> [58] httr_1.4.7           bit_4.6.0            png_0.1-8           
#> [61] ragg_1.5.0           hms_1.1.4            memoise_2.0.1       
#> [64] kableExtra_1.4.0     evaluate_1.0.5       knitr_1.51          
#> [67] GenomicRanges_1.62.1 IRanges_2.44.0       BiocFileCache_3.0.0 
#> [70] viridisLite_0.4.2    rlang_1.1.7          glue_1.8.0          
#> [73] DBI_1.2.3            xml2_1.5.2           BiocGenerics_0.56.0 
#> [76] svglite_2.2.2        rstudioapi_0.18.0    jsonlite_2.0.0      
#> [79] R6_2.6.1             systemfonts_1.3.1    fs_1.6.6
```
