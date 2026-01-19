#' Correlation of two matrices for more efficient calculations
#' 
#' @param mat1 first matrix of data to correlate
#' @param mat2 must have same dimensions as mat1
#' @param type method of correlation
#' @param rank_data logical 
#' @param return_pvalue default=TRUE
#' @return A dataframe with five columns containing the pairs of data that were correlated, the correlation coefficient, p-value, and number of samples included in the correlation.
#' @export
#' @examples
#' # example code
matrixCorrelation = function(mat1, mat2, type = c("pearson", "spearman"), rank_data = FALSE, return_pvalue = TRUE) {
  # mat_ls = list of matrices to conduct pairwise correlation on
  # room for improvement: 
  # - check for same dimensions of matrices
  
  # example of test to compare calculations are accurate:
  # test = matrixCorrelation(counts[,vars1], counts[,vars2])
  # test2 = rcorr(counts)
  # test3 = test |> rowwise() |> mutate(r_compare = test2$r[var1, var2], p_compare = test2$P[var1, var2])
  
  type = match.arg(type)
  
  mat1[is.na(mat2)] = NA
  mat2[is.na(mat1)] = NA
  
  mat_ls = list(mat1, mat2)
  
  if (isTRUE(rank_data) | (type == "spearman")) {
    # in the stats::cor source code, it looks like they just rank the data and then treat it as a pearsons correlation from there. i had trouble getting the rho value output to match up in that sense, so i adjusted this function to be able to carry out the calculation in a qay that matched up with the function I would have looped in the first place. i think this problem was occuring when my data had many 0s. https://github.com/SurajGupta/r-source/blob/a28e609e72ed7c47f6ddfbb86c85279a0750f0b7/src/library/stats/R/cor.R#L70
    mat_ls = lapply(mat_ls, function(x) matrixStats::colRanks(x, ties.method = "average", preserveShape = TRUE))
  }
  
  if (type == "pearson") {
    diff_mat = lapply(mat_ls, function(x) sweep(x, 2, colMeans(x, na.rm = TRUE)))
    sum_sq_diff_mat = lapply(diff_mat, function(x) colSums(x ^ 2, na.rm = TRUE))
    r = colSums(reduce(diff_mat, `*`), na.rm = TRUE) / sqrt(reduce(sum_sq_diff_mat, `*`))
    n = colSums(!is.na(mat1))
  } else if (type == "spearman") {
    # spearman does the manual calculation, whereas ranked pearson ranks the values then puts it into pearson and the coefficient is the same as what returns from cor and rcorr
    # ranked_mats = lapply(list(mat1, mat2), function(x) matrixStats::colRanks(x, ties.method = "average", preserveShape = TRUE))
    diff_mats = reduce(mat_ls, `-`)^2
    n = colSums(!is.na(diff_mats))
    r = 1 - ((6 * (colSums(diff_mats, na.rm = TRUE))) / (n * (n^2 - 1)))
  }
  
  results_df = data.frame(r = r, n = n)
  if (return_pvalue) {
    # calculate pvalue
    # note: p-values for spearman match rcorr, but not cor.test because of different way of calculating test statistic shown in https://github.com/SurajGupta/r-source/blob/a28e609e72ed7c47f6ddfbb86c85279a0750f0b7/src/library/stats/R/cor.test.R#L155
    t_stat = r * sqrt((n - 2) / (1 - r^2))
    p_values = 2 * pt(-abs(t_stat), df = n - 2)
    results_df$P = p_values
  }
  
  if (all(!is.null(c(colnames(mat1), colnames(mat2))))) {
    results_df = cbind(data.frame(var1 = colnames(mat1), var2 = colnames(mat2)), results_df)
  }
  return(results_df)
}

#' Format input
#' 
#' This function pastes the modality at the end of each feature ID and then combines the matrices into one matrix object
#' @param gene_counts count matrix where rownames are sample IDs and colnames are feature ID
#' @export
formatMatrixForCorrelation = function(gene_counts, peak_counts) {
  stopifnot(inherits(gene_counts, c("matrix", "numeric")))
  stopifnot(inherits(peak_counts, "list"))
  stopifnot(any(purrr::map_lgl(peak_counts, ~ inherits(.x, c("matrix", "numeric")))))
  count_mat = c(list(RNASeq = gene_counts), peak_counts) |> 
    purrr::imap(function(.mat, .modality) {
      colnames(.mat) = paste0(colnames(.mat), "_", .modality)
      .mat = as.data.frame(.mat) |> 
        tibble::rownames_to_column(var = "samp_id")
      return(.mat)
    }) |> 
    purrr::reduce(dplyr::full_join, by = "samp_id") |> 
    tibble::column_to_rownames(var = "samp_id") |> 
    as.matrix()
  return(count_mat)
}

#' A wrapper function to carry out correlations in smaller chunks to help save on memory
#' 
#' @param count_mat numeric matrix with features as columns and samples as rows
#' @param correlation_pairs dataframe supplying all pairs of features
#' @param grp_contrast a way to filter to include samples from only certain groups
#' @param rds_fn filename to save the output to for future access
#' @return a dataframe of full results
#' @export
correlateByChromosome = function(count_mat, correlation_pairs, grp_contrast = NULL, rds_fn = NULL) {
  # correlation_pairs is a dataframe with the following columns: ensembl_gene_id, reegulatory element, target_id, link_label, chr, modality_pair
  # split pairs by chromosome to and create distinct list of region pairs to correlate
  chr_pair_ls = correlation_pairs |> select(-any_of(c("ensembl_gene_id"))) |> 
    distinct() |> 
    (\(x) split(x, x$chr))()
  if (!is.null(grp_contrast)) {
    # for acute exercise groups, only include the samples present in the comparison groups
    contrast_var = str_remove_all(grp_contrast, "_.*$")
    groups_compared = formatContrastNames(grp_contrast, "")
    samples_to_retain = .GlobalEnv$sample_info |> 
      filter(!!rlang::sym(contrast_var) %in% groups_compared) |> 
      pull(samp_id)
    if (!all(rownames(count_mat) %in% samples_to_retain)) {
      count_mat = count_mat[samples_to_retain,]
      message("Filtered samples retain ", nrow(count_mat), " samples present in contrast ", grp_contrast)
    }
  }
  # carry out correlation
  cor_results = imap(chr_pair_ls, function(.df, .chr) {
    message("Starting ", .chr)
    tictoc::tic()
    chunk_cutoff = 2e6
    if (nrow(.df) > chunk_cutoff) {
      start_ranges = seq(1, nrow(.df), chunk_cutoff) 
      end_ranges = c((start_ranges - 1)[-1], nrow(.df))
      message("\t...too many pairs, carrying out correlation in ", length(start_ranges), " chunks")
      tmp_cor_res = map2(start_ranges, end_ranges, function(.x, .y) {
        tmp_mat1 = count_mat[,.df$regulatory_element[c(.x:.y)], drop = FALSE]
        tmp_mat2 = count_mat[,.df$target_id[c(.x:.y)], drop = FALSE]
        tmp_cor_res = matrixCorrelation(tmp_mat1, tmp_mat2, "spearman")
      }) |> 
        bind_rows()
    } else {
      tmp_mat1 = count_mat[,.df$regulatory_element, drop = FALSE]
      tmp_mat2 = count_mat[,.df$target_id, drop = FALSE]
      tmp_cor_res = matrixCorrelation(tmp_mat1, tmp_mat2, "spearman")
    }
    tictoc::toc()
    return(tmp_cor_res)
  })
  cor_results = bind_rows(cor_results)
  if (!is.null(rds_fn)) {
    saveRDS(cor_results, rds_fn, compress = FALSE)
  }
  # join the original pair information to the paired correlation results
  full_cor_res = correlation_pairs |> 
    left_join(cor_results |> dplyr::rename("regulatory_element" = "var1", "target_id" = "var2"), by = c("regulatory_element", "target_id"))
  # calculate adjusted pvalue within modality-pair:link-label groups
  p_adj = full_cor_res |>
    select(-ensembl_gene_id) |>
    distinct() |>
    group_by(link_label, modality_pair) |>
    mutate(BH = p.adjust(P, method = "BH"))
  # add those adjusted pvalues back to the full results dataframe
  full_cor_res = full_cor_res |>
    left_join(p_adj)
  return(full_cor_res)
}

#' Annotation of significant promoter peak information
#' @param df the correlation results
#' @return an annotated dataframe 
editCorResFields = function(df) {
  sig_promoter_peaks = df |> 
    filter(BH < 0.05, link_label == "promoter_peak_to_gene") |> 
    select(ensembl_gene_id, regulatory_element) |> 
    distinct() |> 
    mutate(prom_peak_sig_cor_w_gene = TRUE)
  df = df |> 
    left_join(sig_promoter_peaks, by = c("ensembl_gene_id", "target_id" = "regulatory_element")) |> 
    mutate(prom_peak_sig_cor_w_gene = case_when(link_label == "promoter_peak_to_gene" ~ TRUE,
                                                is.na(prom_peak_sig_cor_w_gene) ~ FALSE,
                                                TRUE ~ prom_peak_sig_cor_w_gene))
  return(df)
}

#' Annotation of genes that are significant and provide some additional downstream information for interpreting linked networks
#' @param .sig_cor_genes dataframe that includes genes with significant correlations
#' @param p2g_info dataframe
#' @return an annotated dataframe
categorizeAndAnnotateCorrelationResults = function(.sig_cor_genes, p2g_info) {
  sig_cor_res = .sig_cor_genes |> 
    filter(BH < 0.05) |> 
    left_join(select(p2g_info, ensembl_gene_id, regulatory_element = unique_id, dist_to_tss)) |> 
    editCorResFields() 
  # distal peaks that are correlated with genes AND also correlated with a promoter peak that is correlated with a gene will also have a TRUE value in the prom_peak_sig_cor_w_gene column
  distal_gene_prom_connection = sig_cor_res |> 
    filter(link_label == "distal_peak_to_gene" | prom_peak_sig_cor_w_gene) |> 
    group_by(ensembl_gene_id, regulatory_element) |> 
    mutate(distal_connection = "distal_peak_to_gene" %in% unique(link_label) && sum(prom_peak_sig_cor_w_gene) > 0) |> 
    ungroup() |> 
    filter(link_label == "distal_peak_to_gene", distal_connection) |> 
    select(ensembl_gene_id, regulatory_element, target_id, distal_connection)
  
  # genes where there is a significant correlation between the distal peak and gene, and the distal peak and promoter peak(s), but there is NOT a significant correlation between the promoter peak and gene. in this case, the "distal_no_promoter" column will be TRUE
  distal_no_promoter = .sig_cor_genes |> 
    editCorResFields() |> 
    filter(!prom_peak_sig_cor_w_gene) |> 
    filter(((link_label %in% c("distal_peak_to_gene", "distal_peak_to_promoter_peak")) & BH < 0.05)) |> 
    group_by(ensembl_gene_id, regulatory_element) |> 
    filter(length(unique(link_label)) > 1) |> 
    select(ensembl_gene_id, regulatory_element, target_id) |> 
    mutate(distal_no_promoter = TRUE) |> 
    ungroup()
  
  sig_cor_res = sig_cor_res |> 
    left_join(distal_gene_prom_connection) |> 
    replace_na(list(distal_connection = FALSE)) |> 
    mutate(prom_peak_sig_cor_w_gene = ifelse(link_label == "distal_peak_to_gene" & distal_connection, TRUE, prom_peak_sig_cor_w_gene)) |> 
    select(-distal_connection) |> 
    left_join(distal_no_promoter) |> 
    replace_na(list(distal_no_promoter = FALSE))
  df_info = list(column_descriptions = c("distal_no_promoter" = "This column is TRUE for genes where there is a significant correlation between the distal peak and gene, and the distal peak and promoter peak(s), but there is NOT a significant correlation between the promoter peak and gene.", "prom_peak_sig_cor_w_gene" = "As described by the column name, the value will be TRUE if the relationship between the promoter peak and gene is significant. Additionally, distal peaks that are correlated with genes AND also correlated with a promoter peak that is correlated with a gene will also have a TRUE value in this column"))
  sig_cor_res = sig_cor_res |> 
    `attr<-`("info", df_info)
  return(sig_cor_res)
}

#' Iteration of correlations between distal peaks and other peaks residing near it with significant relationships with said gene
#' @param df dataframe of significant correlation results
#' @param full_cor_df dataframe of full correlation results
#' @param count_mat matrix of counts supplied to correlation analysis
#' @return dataframe with correlation results for that particular iteration
corWithDistalPeak = function(df, full_cor_df, count_mat, iter = 1, iter_batch_cutoff = 5) {
  # peakUD = peak with an unknown direction of association with the gene itself.
  # for a peak with an unknown direction of association with a gene (peakUD), assign a direction based on the correlation between peakUD (that is likely linked to a promoter peak) and other distal peaks that are most significantly linked to the gene expression. the direction of the peakUD to gene is assigned based on the direction of correlation with the distal peak with a more confident directional assignment
  # if iter = 1, peakUD will be correlated with the distal peak that returns the most significant relationship with the gene in question. however, peakUD may not be significantly (pval < 0.05) correlated with the top distal region associated with the gene. in that case, the function should be iterated to the next 'n' peaks until it finds a region that peak UD significantly correlates with. The function should be iterated to until no more correlation options exist.
  # full_cor_df = df BEFORE the first iteration of this function. this should only include significant (BH < 0.05) correlations
  # increase the batch cutoff if the correlation matrix at that point consumes too much memory
  df_disclude = df |> 
    filter(region2gene_dir == "disclude")
  most_sig_peak_per_gene = full_cor_df |> 
    filter(ensembl_gene_id %in% unique(df_disclude$ensembl_gene_id), region2gene_dir != "disclude") |> 
    arrange(ensembl_gene_id, BH) |>
    group_by(ensembl_gene_id) |> 
    (\(x) if (iter < iter_batch_cutoff) filter(x, row_number() == iter) else mutate(x, sig_rank = row_number()) |> filter(sig_rank >= iter))() |> 
    ungroup() |> 
    select(ensembl_gene_id, gene_sig_peak = regulatory_element, gene_sig_dir = region2gene_dir, any_of(c("sig_rank")))
  
  ambiguous_peak_paired_with_top_peak = df_disclude |> 
    select(ensembl_gene_id:region2gene_dir) |> 
    left_join(most_sig_peak_per_gene, by = "ensembl_gene_id", relationship = "many-to-many")
  peaks_to_correlate = ambiguous_peak_paired_with_top_peak |> 
    distinct(regulatory_element, gene_sig_peak) |> 
    drop_na()
  if (nrow(peaks_to_correlate) == 0) {
    return(df)
  }
  cor_res = matrixCorrelation(count_mat[,peaks_to_correlate$regulatory_element,drop=FALSE], 
                              count_mat[,peaks_to_correlate$gene_sig_peak,drop=FALSE]) |> 
    select(regulatory_element = var1, gene_sig_peak = var2, gene_sig_r = r, gene_sig_P = P)
  ambiguous_peak_paired_with_top_peak = ambiguous_peak_paired_with_top_peak |> 
    left_join(cor_res, by = c("regulatory_element", "gene_sig_peak")) |> 
    mutate(region2gene_dir = case_when(gene_sig_P >= 0.05 ~ "disclude",
                                       (gene_sig_r < 0) & (gene_sig_dir == "neg") ~ "pos",
                                       (gene_sig_r > 0) & (gene_sig_dir == "pos") ~ "pos",
                                       (gene_sig_r > 0) & (gene_sig_dir == "neg") ~ "neg",
                                       (gene_sig_r < 0) & (gene_sig_dir == "pos") ~ "neg",
                                       TRUE ~ region2gene_dir))
  region2gene_dir_new = ambiguous_peak_paired_with_top_peak |> 
    filter(region2gene_dir != "disclude") |> 
    (\(x) if (iter < iter_batch_cutoff) mutate(x, iter_n_new = iter / 1000) else x |> slice_min(sig_rank, n = 1, by = c("ensembl_gene_id", "regulatory_element")) |>  mutate(iter_n_new = sig_rank / 1000))() |> 
    select(ensembl_gene_id:target_id, region2gene_dir_new = region2gene_dir, iter_n_new)
  if (!("iter_n" %in% colnames(df))) {
    df$iter_n = NA
  }
  df = df |> 
    left_join(region2gene_dir_new, by = c("ensembl_gene_id", "regulatory_element", "target_id")) |> 
    mutate(region2gene_dir = ifelse(is.na(region2gene_dir_new), region2gene_dir, region2gene_dir_new),
           iter_n = ifelse(is.na(iter_n_new), iter_n, iter_n_new)) |>
    # annotate drop to the ones that didn't have a region to correlate with during this iteration
    left_join(ambiguous_peak_paired_with_top_peak |> filter(is.na(gene_sig_peak)) |> select(1:3) |> distinct() |> mutate(no_match = TRUE), by = c("ensembl_gene_id", "regulatory_element", "target_id")) |> 
    replace_na(list("no_match" = FALSE)) |> 
    mutate(region2gene_dir = ifelse((region2gene_dir == "disclude") & no_match, "drop", region2gene_dir)) |> 
    select(-c(region2gene_dir_new, iter_n_new, no_match))
  if (iter >= iter_batch_cutoff) {
    df = df |> 
      mutate(region2gene_dir = ifelse(region2gene_dir == "disclude", "drop", region2gene_dir))
  }
  return(df)
}

#' Wrapper function for corWithDistalPeak to iterate through the correlations
#' @param full_cor_df dataframe of full correlation results
#' @param count_mat matrix of counts supplied to correlation analysis
#' @param v logical. print informative messages
#' @return dataframe with iterated correlation results until a significant correlation was identified with other distal peaks
corWithDistalPeak_wrapper = function(full_cor_df, count_mat, v = TRUE, .iter_batch_cutoff = 5) {
  df = full_cor_df
  i = 1
  n_discluded = df |> filter(region2gene_dir == "disclude") |> nrow()
  if (v) message("### Correlating ambiguous peaks with other distal peaks assocaited with gene ###\npre fxn: ", df$region2gene_dir |> table() |> (\(x) paste0(names(x), ":\t", scales::comma(as.numeric(x)), collapse = "\t"))())
  while ((n_discluded > 0) & i < 20) {
    df = corWithDistalPeak(df, full_cor_df, count_mat, i, iter_batch_cutoff = .iter_batch_cutoff)
    if (v) message("round ", i, ":\t", df$region2gene_dir |> table() |> (\(x) paste0(names(x), ": ", scales::comma(as.numeric(x)), collapse = "\t"))())
    i = i + 1
    n_discluded = df |> filter(region2gene_dir == "disclude") |> nrow()
  }
  return(df)
}

#' The results from all correlations are used to help assign direction using correlations that have a significant nominal pvalue
#' 
#' Directional assignments are categorized into 3 levels: 
#' 1) Most confident: assign direction using FDR significant links and they are all in agreement
#' 2) Directional disagreement existed, but decision was made by taking the most significant FDR
#' 3) Directional disagreement existed, but decision was made by taking most significant nominal pvalue that was included in comparison
#' 4) These are looser directions associated with the gene - neither distal or promoter region is nominally correlated with the gene, so the direction is based on how that peak is associated with distal peaks that are significantly linked to the gene with high confidence. The iteration of identifying a significant correlation between the peak and another distal peak is divided by 1000 and then added to 4. This shows that the peak:gene relationship exists at confidence level #4, and notes how many distal genes it had to try to correlate with before finding a significant (pval < 0.05) connection.
#' @param sig_cor_res dataframe with the significant results to annotate link direction
#' @param all_cor_res all resulting correlations
#' @param count_mat matrix of counts supplied to correlation analysis
#' @return dataframe with peak-gene link annotation directions with the degree of confidence 
#' @export
assignRegion2GeneDirection = function(sig_cor_res, all_cor_res, count_mat) {
  # the results from all correlations are used to help assign direction using correlations that have a significant nominal pvalue
  
  # directional assignments are categorized into 3 levels: 
  ### 1) Most confident: assign direction using FDR significant links and they are all in agreement
  ### 2) Directional disagreement existed, but decision was made by taking the most significant FDR
  ### 3) Directional disagreement existed, but decision was made by taking most significant nominal pvalue that was included in comparison
  ### 4) These are looser directions associated with the gene - neither distal or promoter region is nominally correlated with the gene, so the direction is based on how that peak is associated with distal peaks that are significantly linked to the gene with high confidence. The iteration of identifying a significant correlation between the peak and another distal peak is divided by 1000 and then added to 4. This shows that the peak:gene relationship exists at confidence level #4, and notes how many distal genes it had to try to correlate with before finding a significant (pval < 0.05) connection.  
  
  #### Step 1: Assign initial relationships based on most accessible information available ####
  full_dir_df = sig_cor_res |> 
    # work with the important info to assign direction
    select(ensembl_gene_id, regulatory_element, target_id, link_label, r, P, BH) |> 
    # add downstream correlation results of promoter peak to gene
    left_join(all_cor_res |> filter(link_label == "promoter_peak_to_gene") |> select(ensembl_gene_id, promoter_peak = regulatory_element, pp_r = r, pp_pval = P, pp_BH = BH), by = c("ensembl_gene_id", "target_id" = "promoter_peak")) |> 
    # add correlation results of regulatory peak itself to the gene. this will be helpful for assigning direction for peaks linked to promoter peaks that are not significantly correlated with gene expression
    left_join(all_cor_res |> filter(link_label %in% c("distal_peak_to_gene", "promoter_peak_to_gene")) |> select(ensembl_gene_id, distal_peak = regulatory_element, dp_r = r, dp_pval = P, dp_BH = BH), by = c("ensembl_gene_id", "regulatory_element" = "distal_peak")) |> 
    mutate(cor_to_compare = case_when(pp_BH < 0.05 ~ "pp",
                                      dp_pval < pp_pval ~ "dp",
                                      TRUE ~ "pp"),
           r_comp = ifelse(cor_to_compare == "pp", pp_r, dp_r),
           pval_comp = ifelse(cor_to_compare == "pp", pp_pval, dp_pval),
           pval_comp = ifelse(is.na(pval_comp), P, pval_comp)) |> 
    # add column that states direction between peak and gene
    mutate(region2gene_dir = case_when(
      # first, if there is a significant (BH < 0.05) correlation between the peak and gene expression, direction should be assigned based on that
      (link_label %in% c("promoter_peak_to_gene", "distal_peak_to_gene")) & (r > 0) ~ "pos",
      (link_label %in% c("promoter_peak_to_gene", "distal_peak_to_gene")) & (r < 0) ~ "neg",
      # if neither the regulatory element of the promoter peak are significantly correlated with gene expression, then you can't get a sense of how the peak is related to the gene. using nominal pvalue because although the values didnt pass the significance threshold for use moving forward, the individual test significance can still provide directional information. we will try to tease direction out of these links later
      (pp_pval >= 0.05) & (dp_pval >= 0.05) ~ "disclude",
      # if distal peak is positively correlated with a promoter peak, but the promoter peak is negatively correlated with gene expression, then the distal peak is negatively associated with the gene, and vice versa. 
      (cor_to_compare == "dp") & (dp_r > 0) ~ "pos",
      (cor_to_compare == "dp") & (dp_r < 0) ~ "neg",
      ((r > 0) & (pp_r > 0)) | ((r < 0) & (pp_r < 0)) ~ "pos",
      ((r > 0) & (pp_r < 0)) | ((r < 0) & (pp_r > 0)) ~ "neg",
      TRUE ~ NA
    ), .before = "pp_r") |> 
    # if there is disagreement between peak:gene direction relationships, then note that here
    mutate(recheck = all(c("pos", "neg") %in% region2gene_dir), .by = c(ensembl_gene_id, regulatory_element),
           dir_high_conf = !recheck & (region2gene_dir != "disclude"))
  
  # region2gene_dir_df is a running dataframe that collects the region:gene relationships with a solidified directional link, noting the confidence score
  region2gene_dir_df = full_dir_df |> 
    filter(!recheck, region2gene_dir != "disclude") |> 
    select(ensembl_gene_id, regulatory_element, region2gene_dir) |> 
    distinct() |> 
    mutate(region2gene_conf = 1)
  
  #### Step 2: Fix directional disagreement ####
  # If there is directional disagreement within a gene:regulatory element group, then correct it by honoring the most significant (BH < 0.05) one
  most_sig_dir = full_dir_df |> 
    filter(recheck) |> 
    group_by(ensembl_gene_id, regulatory_element) |> 
    arrange(BH, .by_group = TRUE) |> 
    # if there is a distal-gene and distal-promoter links, then take the direction from the most FDR sig correlation of that group
    mutate(region2gene_dir = ifelse(all(c("distal_peak_to_gene", "distal_peak_to_promoter_peak") %in% link_label) & (row_number() != 1), NA, region2gene_dir)) |> 
    fill(region2gene_dir) |> 
    ungroup()
  full_dir_df = full_dir_df |> 
    left_join(most_sig_dir |> select(ensembl_gene_id:target_id, region2gene_dir_new = region2gene_dir), by = c("ensembl_gene_id", "regulatory_element", "target_id")) |> 
    mutate(region2gene_dir = ifelse(!is.na(region2gene_dir_new), region2gene_dir_new, region2gene_dir)) |> 
    mutate(recheck = all(c("pos", "neg") %in% region2gene_dir), .by = c(ensembl_gene_id, regulatory_element)) |> 
    select(-region2gene_dir_new)
  region2gene_dir_df = bind_rows(region2gene_dir_df |> select(1:3),
                                 full_dir_df |> filter(!recheck, region2gene_dir != "disclude") |> select(ensembl_gene_id, regulatory_element, region2gene_dir)) |> 
    distinct() |> 
    left_join(region2gene_dir_df, by = c("ensembl_gene_id", "regulatory_element", "region2gene_dir")) |> 
    replace_na(list("region2gene_conf" = 2))
  # for the remaining regulatory elements that have disagreement in direction of the gene relationship, take the direction that has the most significant nominal pvalue 
  most_sig_compare = full_dir_df |> 
    filter(recheck) |> 
    group_by(ensembl_gene_id, regulatory_element) |> 
    arrange(pval_comp, .by_group = TRUE) |> 
    mutate(region2gene_dir_new = ifelse((row_number() != 1), NA, region2gene_dir)) |> 
    fill(region2gene_dir_new) |> 
    ungroup()
  full_dir_df = full_dir_df |> 
    left_join(most_sig_compare |> select(ensembl_gene_id:target_id, region2gene_dir_new), by = c("ensembl_gene_id", "regulatory_element", "target_id")) |> 
    mutate(region2gene_dir = ifelse(!is.na(region2gene_dir_new), region2gene_dir_new, region2gene_dir)) |> 
    mutate(recheck = all(c("pos", "neg") %in% region2gene_dir), .by = c(ensembl_gene_id, regulatory_element)) |> 
    select(-region2gene_dir_new)
  # if there is a link labeled as discluded, but the distal region is labeled elsewhere, then use the other direction label
  discluded_with_dir = full_dir_df |> filter(any(region2gene_dir == "disclude"), .by = c(ensembl_gene_id, regulatory_element)) |> 
    filter(any(c("pos", "neg") %in% region2gene_dir), .by = c(ensembl_gene_id, regulatory_element)) |> 
    filter(region2gene_dir != "disclude") |> 
    distinct(ensembl_gene_id, regulatory_element, region2gene_dir)
  
  full_dir_df = full_dir_df |> 
    left_join(discluded_with_dir |> select(ensembl_gene_id, regulatory_element, region2gene_dir_new = region2gene_dir), by = c("ensembl_gene_id", "regulatory_element")) |> 
    mutate(region2gene_dir = ifelse(!is.na(region2gene_dir_new), region2gene_dir_new, region2gene_dir)) |> 
    select(-region2gene_dir_new)
  
  region2gene_dir_df = bind_rows(region2gene_dir_df |> select(1:3),
                                 full_dir_df |> filter(!recheck, region2gene_dir != "disclude") |> select(ensembl_gene_id, regulatory_element, region2gene_dir)) |> 
    distinct() |> 
    left_join(region2gene_dir_df, by = c("ensembl_gene_id", "regulatory_element", "region2gene_dir")) |> 
    replace_na(list("region2gene_conf" = 3))
  
  full_dir_df = full_dir_df |> 
    select(-c(pp_r:last_col())) |> 
    left_join(region2gene_dir_df, by = c("ensembl_gene_id", "regulatory_element", "region2gene_dir"))
  
  #### Step 3: Salvage discluded relationships ###
  # these are ones that originally end in a dead end, where distal peaks are significantly correlated with promoter peaks that are positioned near the TSS of said gene, but neither the distal peak or the promoter peak are significantly (FDR or nominally) associated with gene expression. however, the distal:promoter relationship still made it into the final results set.
  full_dir_df = corWithDistalPeak_wrapper(full_dir_df, count_mat)
  # update the master dataframe that is cataloging the peak:gene direction and the level of confidence that assignment comes along with
  region2gene_dir_df = bind_rows(region2gene_dir_df,
                                 full_dir_df |> filter(!is.na(iter_n)) |> select(ensembl_gene_id, regulatory_element, region2gene_dir, region2gene_conf = iter_n) |> drop_na() |> mutate(region2gene_conf = region2gene_conf + 4)) |> 
    distinct() 
  region2gene_dir_df = bind_rows(region2gene_dir_df |> select(1:3),
                                 full_dir_df |> select(1:2, region2gene_dir)) |> 
    distinct() |> 
    left_join(region2gene_dir_df, by = c("ensembl_gene_id", "regulatory_element", "region2gene_dir"))
  return(region2gene_dir_df)
}

#' Process correlation results
#' 
#' Genes with a significant correlation with either a promoter or distal peak are retained and focused on further.
#' @export
processCorrelations = function(count_mat, cor_res, correlation_pairs, p2g_info) {
  genes_with_sig_anno = cor_res |> 
    filter((link_label %in% grep("to_gene$", levels(link_label), value = TRUE)) & (P < 0.01)) |> 
    pull(ensembl_gene_id) |> 
    unique()
  sig_cor = cor_res |> 
    filter(ensembl_gene_id %in% genes_with_sig_anno) |> 
    mutate(sig_cor = BH < 0.05) 
  anno_cor_res = categorizeAndAnnotateCorrelationResults(sig_cor, p2g_info)
  region2gene_dir = assignRegion2GeneDirection(anno_cor_res, cor_res, count_mat)
  sig_cor_res = anno_cor_res |> 
    left_join(region2gene_dir, by = c("ensembl_gene_id", "regulatory_element")) |> 
    relocate(ensembl_gene_id, regulatory_element, link_label, modality_pair, r, BH, dist_to_tss, starts_with("region2gene")) #|> 
    # filter(region2gene_dir != "drop")
  return(sig_cor_res)
}
