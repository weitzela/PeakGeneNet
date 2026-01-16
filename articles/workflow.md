# workflow

``` r
library(tidyverse)
library(PeakGeneNet)
```

``` r
p2g_ls = createPeak2GeneObjects(
  gene_counts, peak_counts, 
  biomaRt::useEnsembl(biomart = "genes", 
                      dataset = "rnorvegicus_gene_ensembl", version = 109), 
  "rn7"
)
p2g_ls$paired_df[1:5,]
#>                        unique_id              region_id modality
#> 1 chr3:78634005:78634640_ATACSeq chr3:78634005:78634640  ATACSeq
#> 2 chr3:78652401:78653178_H3K4me1 chr3:78652401:78653178  H3K4me1
#> 3 chr3:78652614:78653063_ATACSeq chr3:78652614:78653063  ATACSeq
#> 4 chr3:78683169:78683832_H3K27ac chr3:78683169:78683832  H3K27ac
#> 5 chr3:78683236:78683697_H3K4me1 chr3:78683236:78683697  H3K4me1
#>      ensembl_gene_id dist_to_tss promoter_peak
#> 1 ENSRNOG00000000008     -977078         FALSE
#> 2 ENSRNOG00000000008     -958540         FALSE
#> 3 ENSRNOG00000000008     -958655         FALSE
#> 4 ENSRNOG00000000008     -927886         FALSE
#> 5 ENSRNOG00000000008     -928021         FALSE
p2g_ls$all_pairs[1:5,]
#>      ensembl_gene_id             regulatory_element
#> 1 ENSRNOG00000000008 chr3:78634005:78634640_ATACSeq
#> 2 ENSRNOG00000000008 chr3:78652401:78653178_H3K4me1
#> 3 ENSRNOG00000000008 chr3:78652614:78653063_ATACSeq
#> 4 ENSRNOG00000000008 chr3:78683169:78683832_H3K27ac
#> 5 ENSRNOG00000000008 chr3:78683236:78683697_H3K4me1
#>                        target_id                   link_label  chr
#> 1 chr3:79610922:79613046_H3K4me3 distal_peak_to_promoter_peak chr3
#> 2 chr3:79610922:79613046_H3K4me3 distal_peak_to_promoter_peak chr3
#> 3 chr3:79610922:79613046_H3K4me3 distal_peak_to_promoter_peak chr3
#> 4 chr3:79610922:79613046_H3K4me3 distal_peak_to_promoter_peak chr3
#> 5 chr3:79610922:79613046_H3K4me3 distal_peak_to_promoter_peak chr3
```
