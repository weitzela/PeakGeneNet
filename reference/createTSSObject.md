# Setup gene TSS object

This creates a dataframe with canonical TSS for each gene that is
included in the analysis

## Usage

``` r
createTSSObject(ensembl_ids, biomart_ensembl)
```

## Examples

``` r
# createTSSObject(c("ENSR00129", "ENSR00139"), biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 105))
```
