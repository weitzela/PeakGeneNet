# Package index

## All functions

- [`adjustCounts()`](https://weitzela.github.io/PeakGeneNet/reference/adjustCounts.md)
  : Adjust Counts
- [`adjustCountsMultipleContrasts()`](https://weitzela.github.io/PeakGeneNet/reference/adjustCountsMultipleContrasts.md)
  : Adjust counts and protect for multiple biological contrasts Wrapper
  function to iterate through multiple biological contrasts to protect.
  This first uses the limma::removeBatchEffect function to remove
  technical variables such as preparation batches and flowcell. Then, it
  applies the adjustCovariateMatrix to further adjust for other known
  biological contrasts
- [`adjustCovariateMatrix()`](https://weitzela.github.io/PeakGeneNet/reference/adjustCovariateMatrix.md)
  : Use of matrix functions to adjust data via linear model
- [`assignRegion2GeneDirection()`](https://weitzela.github.io/PeakGeneNet/reference/assignRegion2GeneDirection.md)
  : The results from all correlations are used to help assign direction
  using correlations that have a significant nominal pvalue
- [`categorizeAndAnnotateCorrelationResults()`](https://weitzela.github.io/PeakGeneNet/reference/categorizeAndAnnotateCorrelationResults.md)
  : Annotation of genes that are significant and provide some additional
  downstream information for interpreting linked networks
- [`corWithDistalPeak()`](https://weitzela.github.io/PeakGeneNet/reference/corWithDistalPeak.md)
  : Iteration of correlations between distal peaks and other peaks
  residing near it with significant relationships with said gene
- [`corWithDistalPeak_wrapper()`](https://weitzela.github.io/PeakGeneNet/reference/corWithDistalPeak_wrapper.md)
  : Wrapper function for corWithDistalPeak to iterate through the
  correlations
- [`correlateByChromosome()`](https://weitzela.github.io/PeakGeneNet/reference/correlateByChromosome.md)
  : A wrapper function to carry out correlations in smaller chunks to
  help save on memory
- [`createPeak2GeneObjects()`](https://weitzela.github.io/PeakGeneNet/reference/createPeak2GeneObjects.md)
  : Set up genomicranges objects for the data included in the
  correlation analysis
- [`editCorResFields()`](https://weitzela.github.io/PeakGeneNet/reference/editCorResFields.md)
  : Annotation of significant promoter peak information
- [`gene_counts`](https://weitzela.github.io/PeakGeneNet/reference/gene_counts.md)
  : Example Gene Counts
- [`matrixCorrelation()`](https://weitzela.github.io/PeakGeneNet/reference/matrixCorrelation.md)
  : Correlation of two matrices for more efficient calculations
- [`peak_counts`](https://weitzela.github.io/PeakGeneNet/reference/peak_counts.md)
  : Example Peak Counts
- [`removeBatchEffect_wrapper()`](https://weitzela.github.io/PeakGeneNet/reference/removeBatchEffect_wrapper.md)
  : Wrapper: Remove Batch Effect This supplies the experimental model to
  remove technical variables, while protecting for the experimental
  variables.
