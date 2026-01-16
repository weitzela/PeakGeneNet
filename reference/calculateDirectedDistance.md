# Genomic distance with the direction

This function assigns a negative value if the subject is upstream from
the target and positive value if the subject is downstream from the
target, while taking into account if the gene is located on the negative
or positive strand. This function is created because the
GenomicRanges::distance function itself will not account for direction
from TSS (upstream (-) or downstream (+)).

## Usage

``` r
calculateDirectedDistance(subject, target)
```

## Value

A numeric vector of signed distances

## Examples

``` r
# example code
```
