# Create Peak-Gene Links Set up genomicranges objects for the data included in the correlation analysis. Data matrices must have *unique* feature IDs.

Create Peak-Gene Links Set up genomicranges objects for the data
included in the correlation analysis. Data matrices must have *unique*
feature IDs.

## Usage

``` r
createPeak2GeneObjects(genes, peaks, biomart_ensembl, ucsc_genome)
```

## Arguments

- genes:

  either a character vector of ensembl IDs or a matrix with the row or
  column names as the ensembl IDs

- peaks:

  a named list of character vectors containing renomic region. Genomic
  regions should begin with "chr". The list can also contain count
  matrices, where row names or column names are genomic regions.

## Examples
