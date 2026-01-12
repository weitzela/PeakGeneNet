#' Use of matrix functions to adjust data via linear model
#' @param counts matrix of transformed counts (vst, log, inverse rank normalized, etc.). Should have sample IDs as rownames, feature ID as colnames
#' @param covariate_df dataframe with all covariates, sample names should be rownames. should be in the same order as the input matrix
#' @param vars_to_protect character vector of column names that should be protected. the coefficients calculated for these variables will not be subtracted from the returned values
#' @param return_coefficients logcal. default=FALSE
#' @return matrix of adjusted counts 
#' @export
adjustCovariateMatrix = function(counts, covariate_df, vars_to_protect, return_coefficients = FALSE) {
  # counts: matrix, should have sample IDs as rownames, feature ID as colnames
  # covariate_df: dataframe with all covariates, sample names should be rownames. should be in the same order as the input matrix
  vars_to_disclude_from_adj = paste0(unique(c("intercept", vars_to_protect)), collapse = "|")
  model_matrix = covariate_df %>% 
    droplevels() %>%
    model.matrix(reformulate(colnames(.)), data = .)
  coefficients = lm.fit(model_matrix, counts)$coefficients
  coefficients_sums = (model_matrix %>% .[,-grep(vars_to_disclude_from_adj, colnames(.), ignore.case = TRUE),drop = FALSE]) %*% (coefficients %>% .[-grep(vars_to_disclude_from_adj, rownames(.), ignore.case = TRUE),,drop = FALSE])
  adj_counts = counts - coefficients_sums
  if (return_coefficients) {
    adj_counts = adj_counts %>% 
      `attr<-`("coefficients", coefficients)
  }
  return(adj_counts)
}

#' Wrapper: Remove Batch Effect
#' This supplies the experimental model to remove technical variables, while protecting for the experimental variables.
removeBatchEffect_wrapper = function(counts, df, var_adj, var_protect = c("line", "sex", "grp", "train.sed")) {
  covariate_df = df %>% select(any_of(var_adj)) #%>% mutate(across(where(is.factor), as.numeric))
  adj_formula =  paste0(colnames(covariate_df), collapse = " + ") %>% paste0("~ ", .)
  if (!all(var_adj %in% colnames(covariate_df))) {
    warning("COVARIATES: ", knitr::combine_words(setdiff(var_adj, colnames(covariate_df))), " are not columns in the data.frame provided and will not be adjusted for.\n")
    if (ncol(covariate_df) < 1) {
      warning("COVARIATES: There are no variables to adjust for.\n")
      covariate_df = NULL
      adj_formula = NULL
    }
  } 
  model_df = df %>% select(any_of(var_protect))
  if (!all(var_protect %in% colnames(model_df))) {
    warning("MODEL MATRIX: ", knitr::combine_words(setdiff(var_protect, colnames(model_df))), " are not columns in the data.frame provided and will not be protected for.\n")
  }
  formula = paste0(colnames(model_df), collapse = " + ")
  # covariate_df = df %>% select()
  counts_adj = limma::removeBatchEffect(x = counts, 
                                        batch = if ("prep_date" %in% colnames(covariate_df)) covariate_df[["prep_date"]] else covariate_df[[""]], 
                                        batch2 = if ("flowcell" %in% colnames(covariate_df)) covariate_df[["flowcell"]] else covariate_df[[""]],
                                        covariates = covariate_df %>% {if (inherits(., "NULL")) . else select(., -any_of(c("prep_date", "flowcell")))}, design = model.matrix(reformulate(formula), data = model_df))
  counts_adj = counts_adj %>% 
    `attr<-`("covariate_df", covariate_df) %>% 
    `attr<-`("adj_formula", adj_formula) %>% 
    `attr<-`("model_matrix_formula", formula)
  return(counts_adj)
}

#' Adjust counts and protect for multiple biological contrasts
#' Wrapper function to iterate through multiple biological contrasts to protect. This first uses the limma::removeBatchEffect function to remove technical variables such as preparation batches and flowcell. Then, it applies the adjustCovariateMatrix to further adjust for other known biological contrasts
adjustCountsMultipleContrasts = function(count_mat, data_label, contrast_list = c("line", "sex", "grp"), concat_data_label = TRUE, remove_batch_effect = TRUE) {
  # count_mat: formatted with features as the rownames and sample IDs as the column name. the 'sample_info' attribute must have rownames that match up with the colnames and be in the same order
  if (remove_batch_effect & !is.null(attr(count_mat, "remove_batch_effect"))) {
    if (attr(count_mat, "remove_batch_effect") == FALSE) remove_batch_effect = FALSE
  }
  if (remove_batch_effect) {
    adj_counts = removeBatchEffect_wrapper(count_mat, attr(count_mat, "sample_info"), 
                                           c("prep_date", "flowcell", colnames(select(attr(count_mat, "sample_info"), starts_with("W_")))))
  } else {
    # prep date and flowcell are not adjusted for when putting transformed posterior probabilities through this function
    message(data_label, ": skipping fxn removeBatchEffect_wrapper")
    adj_counts = count_mat
  }
  adj_counts = purrr::set_names(contrast_list) |> 
    map(function(.protect_var) {
      vars_to_adjust = setdiff(c("line", "sex", "grp", "train.sed"), .protect_var)
      message("protecting ", .protect_var, " and adjusting for ", knitr::combine_words(vars_to_adjust))
      adj_counts = adjustCovariateMatrix(t(adj_counts), attr(adj_counts, "sample_info") %>% select(any_of(vars_to_adjust)), vars_to_protect = .protect_var)
      if (concat_data_label) {
        colnames(adj_counts) = paste0(colnames(adj_counts), "_", data_label)
      }
      if (!((data_label == "RNASeq") | (get_modality(data_label) == "ChromHMM"))) rownames(adj_counts) = pull(metadata, samp_id, name = "library")[rownames(adj_counts)]
      adj_counts = rownames_to_column(as.data.frame(adj_counts), var = "samp_id")
      return(adj_counts)
    })
  if (length(contrast_list) == 1) {
    adj_counts = adj_counts[[1]]
  }
  return(adj_counts)
}

#' Adjust Counts
#' @export
adjustCounts = function(count_ls, contrast_list = c("line", "sex", "grp"), .remove_batch_effect = TRUE) {
  adj_counts = imap(count_ls, ~ adjustCountsMultipleContrasts(.x, .y, contrast_list, remove_batch_effect = .remove_batch_effect)) |> 
    purrr::transpose() |> 
    map(function(.x) {
      .x |> 
        purrr::reduce(dplyr::full_join, by = "samp_id") |>
        column_to_rownames(var = "samp_id") |>
        as.matrix()
    })
  return(adj_counts)
}