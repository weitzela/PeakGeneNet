# Adjust counts and protect for multiple biological contrasts Wrapper function to iterate through multiple biological contrasts to protect. This first uses the limma::removeBatchEffect function to remove technical variables such as preparation batches and flowcell. Then, it applies the adjustCovariateMatrix to further adjust for other known biological contrasts

Adjust counts and protect for multiple biological contrasts Wrapper
function to iterate through multiple biological contrasts to protect.
This first uses the limma::removeBatchEffect function to remove
technical variables such as preparation batches and flowcell. Then, it
applies the adjustCovariateMatrix to further adjust for other known
biological contrasts

## Usage

``` r
adjustCountsMultipleContrasts(
  count_mat,
  data_label,
  contrast_list = c("line", "sex", "grp"),
  concat_data_label = TRUE,
  remove_batch_effect = TRUE
)
```
