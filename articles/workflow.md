# workflow

``` r
library(tidyverse)
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.6
#> ✔ forcats   1.0.1     ✔ stringr   1.6.0
#> ✔ ggplot2   4.0.1     ✔ tibble    3.3.1
#> ✔ lubridate 1.9.4     ✔ tidyr     1.3.2
#> ✔ purrr     1.2.1     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(PeakGeneNet)
#> Warning: replacing previous import 'dplyr::between' by 'plyranges::between'
#> when loading 'PeakGeneNet'
#> Warning: replacing previous import 'dplyr::n_distinct' by
#> 'plyranges::n_distinct' when loading 'PeakGeneNet'
#> Warning: replacing previous import 'dplyr::n' by 'plyranges::n' when loading
#> 'PeakGeneNet'
#> Warning: replacing previous import 'plyranges::between' by 'dplyr::between'
#> when loading 'PeakGeneNet'
#> Warning: replacing previous import 'plyranges::n' by 'dplyr::n' when loading
#> 'PeakGeneNet'
#> Warning: replacing previous import 'plyranges::n_distinct' by
#> 'dplyr::n_distinct' when loading 'PeakGeneNet'
```

``` r
p2g_ls = createPeak2GeneObjects(gene_counts, peak_counts, biomaRt::useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl", version = 109), "rn7")
```
