# Setup gene TSS object

This creates a dataframe with canonical TSS for each gene that is
included in the analysis

## Usage

``` r
createTSSGr(ensembl_ids, biomart_ensembl, ucsc_genome)
```

## Arguments

- ucsc_genome:

  genome release formatting in ucsc style. For example: hg38, rn7, etc.

## Examples

``` r
# createTSSObject(c("ENSR00129", "ENSR00139"), biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 105), "hg38")
```
