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
  |                    |   R18804 |   R18805 |   R18809 |
  |:-------------------|---------:|---------:|---------:|
  | ENSRNOG00000000008 | 6.226779 | 6.150745 | 6.388451 |
  | ENSRNOG00000000082 | 7.267295 | 7.324605 | 7.220238 |
  | ENSRNOG00000001489 | 7.768652 | 8.032395 | 7.814327 |
- Peak counts: a list object with one or more matrices. Matrices should
  be in same format as gene counts, with row names as genomic region
  where the chromosome names, start, and end locations are separaated by
  a punctuation mark (e.g., chr1:123:456 or chr1:123-456).
  |                          | PL119617 | PL119618 | PL119619 |
  |:-------------------------|---------:|---------:|---------:|
  | chr2:103325658:103326019 | 6.616630 | 6.816892 | 6.210107 |
  | chr2:103337000:103337431 | 7.643895 | 7.417112 | 7.204352 |
  | chr2:103483778:103484423 | 7.449124 | 7.336309 | 7.031643 |
  | chr3:78634005:78634640   | 6.754136 | 6.966854 | 7.031643 |
  | chr3:78652614:78653063   | 6.707064 | 6.966854 | 6.674483 |
  | chr3:78702369:78702924   | 6.890375 | 7.033393 | 7.272677 |

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

Then, carry out and process the correlations using the following
functions. The count matrices to for used as the input for the
correlation analysis should have samples in rows and features in
columns.

``` r
tmp_mat = t(gene_counts)[,c(1,1)]
!duplicated(colnames(tmp_mat))
#> [1]  TRUE FALSE
tmp_mat[,-duplicated(colnames(tmp_mat)),drop=FALSE]
#>        ENSRNOG00000000008
#> R18804           6.226779
#> R18805           6.150745
#> R18809           6.388451
#> R18810           6.201142
#> R18819           6.328846
#> R18820           6.287093
#> R18821           5.827359
#> R18827           6.009111
#> R18828           6.206376
#> R18829           6.216744
#> R18839           6.284814
#> R18840           6.560351
#> R18841           6.306194
#> R18842           6.200470
#> R18843           6.302628
#> R18850           6.300516
#> R18851           6.154064
#> R18852           6.218144
#> R18856           6.116110
#> R18857           6.163643
#> R18858           6.473600
#> R18859           6.092026
#> R18860           6.256806
#> R18861           6.135792
#> R18879           6.300636
#> R18880           6.510249
#> R18885           6.139304
#> R18886           6.323827
#> R18887           6.240320
#> R18888           5.841920
#> R18889           6.464813
#> R18890           5.969445
#> R18900           6.208546
#> R18901           6.240876
#> R18909           6.061753
#> R18910           6.390905
#> R18917           6.200783
#> R18918           6.373255
#> R18919           6.354534
#> R18933           6.389630
#> R18934           6.294382
#> R18935           6.357744
#> R18936           6.019776
#> R18949           5.526674
#> R18950           6.372019
#> R18951           6.678617
#> R18954           6.494188
#> R18955           6.221581
#> R18956           6.262311
#> R18962           6.570335
#> R18963           6.221500
#> R18972           6.246718
#> R18973           6.244364
#> R18974           6.323325
#> R18975           6.462049
#> R18976           6.361749
#> R18979           6.204799
#> R18980           6.346200
#> R18981           6.470114
#> R18997           6.323925
#> R18998           6.120808
#> R18999           6.019391
#> R19014           6.140557
#> R19015           6.569781
#> R28311           6.017178
#> R28312           6.152096
#> R28313           6.167043
#> R28321           6.262532
#> R28322           6.291001
#> R28327           6.331873
#> R28328           6.373952
#> R28329           6.014881
#> R28336           6.164190
#> R28337           6.246948
#> R28338           5.911523
#> R28339           6.107326
#> R28340           5.953431
#> R28348           6.476898
#> R28349           6.215214
#> R28353           6.198157
#> R28354           6.139912
#> R28358           6.213172
#> R28359           6.239303
#> R28375           6.463560
#> R28376           6.542750
#> R28377           6.251136
#> R28381           6.312130
#> R28382           5.994931
#> R28383           5.735034
#> R28387           6.129441
#> R28388           6.119432
#> R28393           6.417523
#> R28394           6.306894
#> R28395           6.238548
#> R28400           6.182375
#> R28401           6.135212
#> R28402           6.032087
#> R28408           6.298229
#> R28409           6.053847
#> R28416           5.998940
#> R28417           6.138375
#> R28423           6.222035
#> R28424           6.326649
#> R28425           6.135861
#> R28427           6.282609
#> R28434           6.427184
#> R28436           6.434721
#> R28437           6.395132
#> R28453           6.285133
#> R28454           5.911657
#> R28456           6.371622
#> R28457           6.239782
#> R28463           6.127067
#> R28464           6.287524
#> R28470           6.235805
#> R28471           6.389448
#> R28481           6.335875
#> R28483           6.081721
#> R28489           6.195184
#> R28490           6.361663
#> R28503           6.260248
#> R28504           6.331142
#> R28511           6.439461
#> R28512           6.186302
#> R28514           6.175235
#> R28515           6.028437
#> R28518           6.257430
#> R28519           5.991977
```

``` r
peak_counts_t = map(peak_counts, function(.x) {
  colnames(.x) = attr(.x, "samp")[colnames(.x), "samp_id"]
  return(t(.x))
})
count_mat = formatMatrixForCorrelation(t(gene_counts), peak_counts_t)
count_mat[1:5,1:5] |> kableExtra::kable()
```

|        | ENSRNOG00000000008_RNASeq | ENSRNOG00000000082_RNASeq | ENSRNOG00000001489_RNASeq | chr3:78683169:78683832_H3K27ac | chr3:78708054:78708512_H3K27ac |
|:-------|--------------------------:|--------------------------:|--------------------------:|-------------------------------:|-------------------------------:|
| R18804 |                  6.226779 |                  7.267295 |                  7.768652 |                       5.974274 |                       5.830707 |
| R18805 |                  6.150745 |                  7.324605 |                  8.032395 |                       6.015659 |                       5.922595 |
| R18809 |                  6.388451 |                  7.220238 |                  7.814327 |                       6.018896 |                       5.880515 |
| R18810 |                  6.201142 |                  7.341687 |                  7.850113 |                       6.069272 |                       6.599047 |
| R18819 |                  6.328846 |                  7.358075 |                  7.973985 |                             NA |                             NA |

``` r

full_cor_res = correlateByChromosome(count_mat, p2g_ls$correlation_pairs)
#> Starting chr1
#> 0.021 sec elapsed
#> Starting chr2
#> 0.007 sec elapsed
#> Starting chr3
#> 0.331 sec elapsed
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

processed_cor = processCorrelations(count_mat, full_cor_res, p2g_ls$correlation_pairs, p2g_ls$p2g_info)
#> Joining with `by = join_by(ensembl_gene_id, regulatory_element)`
#> Joining with `by = join_by(ensembl_gene_id, regulatory_element, target_id)`
#> Joining with `by = join_by(ensembl_gene_id, regulatory_element, target_id)`
#> ### Correlating ambiguous peaks with other distal peaks assocaited with gene
#> ### pre fxn: disclude: 74 neg: 22 pos: 22
#> round 1: disclude: 2 drop: 27 neg: 44 pos: 45
#> round 2: disclude: 1 drop: 27 neg: 44 pos: 46
#> round 3: disclude: 1 drop: 27 neg: 44 pos: 46
#> round 4: disclude: 1 drop: 27 neg: 44 pos: 46
#> round 5: drop: 28 neg: 44 pos: 46
processed_cor |> head(5) |> kableExtra::kable()
```

| ensembl_gene_id    | regulatory_element             | link_label                   | modality_pair   |          r |        BH | dist_to_tss | region2gene_dir | region2gene_conf | target_id                      | chr  |   n |         P | sig_cor | prom_peak_sig_cor_w_gene | distal_no_promoter |
|:-------------------|:-------------------------------|:-----------------------------|:----------------|-----------:|----------:|------------:|:----------------|-----------------:|:-------------------------------|:-----|----:|----------:|:--------|:-------------------------|:-------------------|
| ENSRNOG00000000008 | chr3:78634005:78634640_ATACSeq | distal_peak_to_promoter_peak | ATACSeq-ATACSeq |  0.3013569 | 0.0043280 |     -977078 | pos             |            4.001 | chr3:79611428:79612928_ATACSeq | chr3 | 123 | 0.0007061 | TRUE    | FALSE                    | FALSE              |
| ENSRNOG00000000008 | chr3:78702369:78702924_ATACSeq | distal_peak_to_promoter_peak | ATACSeq-ATACSeq |  0.3539229 | 0.0006360 |     -908794 | pos             |            4.001 | chr3:79611428:79612928_ATACSeq | chr3 | 123 | 0.0000592 | TRUE    | FALSE                    | FALSE              |
| ENSRNOG00000000008 | chr3:78713352:78714085_ATACSeq | distal_peak_to_promoter_peak | ATACSeq-ATACSeq |  0.2480879 | 0.0226489 |     -897633 | pos             |            4.001 | chr3:79611428:79612928_ATACSeq | chr3 | 123 | 0.0056622 | TRUE    | FALSE                    | FALSE              |
| ENSRNOG00000000008 | chr3:78734426:78734839_ATACSeq | distal_peak_to_promoter_peak | ATACSeq-ATACSeq | -0.3866776 | 0.0001516 |     -876879 | neg             |            4.001 | chr3:79611428:79612928_ATACSeq | chr3 | 123 | 0.0000100 | TRUE    | FALSE                    | FALSE              |
| ENSRNOG00000000008 | chr3:78783524:78783997_ATACSeq | distal_peak_to_promoter_peak | ATACSeq-ATACSeq |  0.2642814 | 0.0135542 |     -827721 | pos             |            4.001 | chr3:79611428:79612928_ATACSeq | chr3 | 123 | 0.0031389 | TRUE    | FALSE                    | FALSE              |

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
#> [13] lifecycle_1.0.5      KEGGREST_1.50.0      RSQLite_2.4.5       
#> [16] magrittr_2.0.4       compiler_4.5.2       rlang_1.1.7         
#> [19] sass_0.4.10          progress_1.2.3       tools_4.5.2         
#> [22] yaml_2.3.12          knitr_1.51           prettyunits_1.2.0   
#> [25] bit_4.6.0            curl_7.0.0           xml2_1.5.2          
#> [28] RColorBrewer_1.1-3   withr_3.0.2          BiocGenerics_0.56.0 
#> [31] desc_1.4.3           grid_4.5.2           stats4_4.5.2        
#> [34] scales_1.4.0         gtools_3.9.5         biomaRt_2.66.0      
#> [37] cli_3.6.5            rmarkdown_2.30       crayon_1.5.3        
#> [40] ragg_1.5.0           generics_0.1.4       rstudioapi_0.18.0   
#> [43] httr_1.4.7           tzdb_0.5.0           DBI_1.2.3           
#> [46] cachem_1.1.0         AnnotationDbi_1.72.0 XVector_0.50.0      
#> [49] matrixStats_1.5.0    vctrs_0.7.0          jsonlite_2.0.0      
#> [52] IRanges_2.44.0       hms_1.1.4            S4Vectors_0.48.0    
#> [55] bit64_4.6.0-1        systemfonts_1.3.1    jquerylib_0.1.4     
#> [58] glue_1.8.0           pkgdown_2.2.0        stringi_1.8.7       
#> [61] gtable_0.3.6         GenomeInfoDb_1.46.2  GenomicRanges_1.62.1
#> [64] UCSC.utils_1.6.1     pillar_1.11.1        rappdirs_0.3.4      
#> [67] htmltools_0.5.9      Seqinfo_1.0.0        R6_2.6.1            
#> [70] dbplyr_2.5.1         httr2_1.2.2          textshaping_1.0.4   
#> [73] evaluate_1.0.5       kableExtra_1.4.0     Biobase_2.70.0      
#> [76] png_0.1-8            tictoc_1.2.1         memoise_2.0.1       
#> [79] snakecase_0.11.1     bslib_0.9.0          svglite_2.2.2       
#> [82] xfun_0.56            fs_1.6.6             pkgconfig_2.0.3
```
