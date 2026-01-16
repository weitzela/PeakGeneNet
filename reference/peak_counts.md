# Example Peak Counts

A list of matrices containing VST transformed counts for multiple
modalities to be integrated. The matrices include features (peaks) as
rows and samples as columns. Each matrix contains a "sample_info"
attribute that is a dataframe which provides metadata information to be
used for data normalization. The rownames should be identical to the
column names for the respective matrix. The dataframe should only
include numeric and factor data.

## Usage

``` r
peak_counts
```

## Format

A list of matrices

- H3K27ac:

  A matrix of H3K27ac count data with 124 rows and 59 columns

- H3K4me1:

  A matrix of H3K4me1 count data with 302 rows and 59 columns

- H3K4me3:

  A matrix of H3K4me3 count data with 5 rows and 59 columns

- ATACSeq:

  A matrix of ATACSeq count data with 296 rows and 123 columns

## Source

University of Michigan
