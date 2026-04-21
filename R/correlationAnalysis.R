#' Compute pairwise correlations between two aligned matrices
#'
#' Correlates column `i` of `mat1` against column `i` of `mat2` for all `i`,
#' returning results in a single vectorised pass rather than looping over
#' individual pairs. This is substantially faster than calling `cor()` or
#' `cor.test()` in a loop when there are thousands of peak–gene pairs.
#'
#' @param mat1 Numeric matrix. Features as columns, samples as rows.
#' @param mat2 Numeric matrix. Must have the same dimensions as `mat1` and the
#'   same row order (samples aligned).
#' @param type Correlation method: `"pearson"` or `"spearman"`. Spearman is
#'   used throughout the PeakGeneNet pipeline. Default: `"pearson"`.
#' @param rank_data Logical. If `TRUE`, columns are rank-transformed before
#'   correlation regardless of `type`. When `type = "spearman"` this is done
#'   automatically; set `rank_data = TRUE` to force ranking with
#'   `type = "pearson"`. Default: `FALSE`.
#' @param return_pvalue Logical. If `TRUE` (default), p-values are computed
#'   from the t-statistic and appended as column `P`. P-values for Spearman
#'   match `Hmisc::rcorr()` but differ from `cor.test()` due to a different
#'   test-statistic formula.
#'
#' @return A data frame with columns `var1`, `var2` (column names of `mat1`
#'   and `mat2`), `r` (correlation coefficient), `n` (number of non-missing
#'   sample pairs), and `P` (p-value, if `return_pvalue = TRUE`).
#' @export
matrixCorrelation = function(mat1, mat2, type = c("pearson", "spearman"), rank_data = FALSE, return_pvalue = TRUE) {
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

#' Combine gene and peak count matrices into a single matrix for correlation
#'
#' Appends the modality name to each feature ID (e.g. `"ENSRNOG00000000008"`
#' becomes `"ENSRNOG00000000008_RNASeq"`, and peak IDs gain their modality
#' suffix) and then joins all matrices on shared sample IDs. The resulting
#' combined matrix is the `count_mat` input expected by
#' [correlateByChromosome()].
#'
#' @param gene_counts Numeric matrix of gene expression counts with sample IDs
#'   as row names and Ensembl gene IDs as column names. Input matrices should
#'   be transformed (e.g., VST, log, or inverse-rank normalisation) before
#'   calling this function.
#' @param peak_counts A named list of numeric peak count matrices, each with
#'   sample IDs as row names and peak region IDs as column names. Names must
#'   match the modality labels used in [createPeak2GeneObjects()]. Input
#'   matrices should be transformed consistently with `gene_counts`.
#'
#' @return A single numeric matrix with samples as rows and all features
#'   (genes and peaks, modality-suffixed) as columns. Samples present in only
#'   a subset of matrices will have `NA` values for missing modalities.
#' @export
formatMatrixForCorrelation = function(gene_counts, peak_counts) {
  stopifnot(inherits(gene_counts, c("matrix", "numeric")))
  stopifnot(inherits(peak_counts, "list"))
  stopifnot(any(purrr::map_lgl(peak_counts, ~ inherits(.x, c("matrix", "numeric")))))
  count_mat = c(list(RNASeq = gene_counts), peak_counts) |> 
    purrr::imap(function(.mat, .modality) {
      duplicated_features = duplicated(colnames(.mat))
      if (any(duplicated_features)) {
        message(sum(duplicated_features), " duplicated feature(s) from the ", .modality, " matrix removed, based on repeated column IDs.")
        .mat = .mat[,!duplicated_features,drop=FALSE]
      }
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

#' Run Spearman correlations for all peak–gene pairs, chromosome by chromosome
#'
#' Iterates over chromosomes, correlating each `regulatory_element`–`target_id`
#' pair in `correlation_pairs` using [matrixCorrelation()]. Large chromosomes
#' are processed in chunks of 2 million pairs to limit peak memory use. After
#' all correlations are computed, BH-adjusted p-values are calculated within
#' groups defined by `link_label` × `modality_pair` — these groups define the
#' exchangeable family for multiple testing correction.
#'
#' @param count_mat Numeric matrix returned by [formatMatrixForCorrelation()]:
#'   samples as rows, modality-suffixed feature IDs as columns.
#' @param correlation_pairs Data frame returned by [createPeak2GeneObjects()]:
#'   must have columns `ensembl_gene_id`, `regulatory_element`, `target_id`,
#'   `link_label`, `chr`, and `modality_pair`.
#' @param grp_contrast Optional character string specifying a contrast used to
#'   subset samples before correlating (e.g. `"treatment_A_vs_B"`). When
#'   provided, the function looks up matching samples from a `sample_info` data
#'   frame in the global environment via `formatContrastNames()`. Default:
#'   `NULL` (all samples in `count_mat` are used).
#' @param rds_fn Optional file path. If provided, the raw correlation results
#'   (before joining back to `correlation_pairs`) are saved with `saveRDS()`
#'   for later retrieval. Default: `NULL`.
#'
#' @return A data frame with one row per `regulatory_element`–`target_id` pair
#'   per gene. Columns:
#'   \describe{
#'     \item{`ensembl_gene_id`}{Ensembl gene ID for the gene this link belongs
#'       to.}
#'     \item{`regulatory_element`}{Unique peak ID (`unique_id` from
#'       [createPeak2GeneObjects()]) of the peak acting as the regulatory
#'       element.}
#'     \item{`target_id`}{ID of the correlation target: either a peak
#'       `unique_id` (for peak–peak links) or an Ensembl gene ID suffixed with
#'       `"_RNASeq"` (for peak–gene links).}
#'     \item{`link_label`}{Factor indicating the link type. One of
#'       `"promoter_peak_to_gene"`, `"distal_peak_to_gene"`,
#'       `"distal_peak_to_promoter_peak"`, `"promoter_peak_to_promoter_peak"`.}
#'     \item{`chr`}{Chromosome of the gene TSS.}
#'     \item{`modality_pair`}{Factor. Hyphen-separated modality labels for the
#'       regulatory element and target, alphabetically ordered
#'       (e.g. `"ATACSeq-RNASeq"`). Used to define the multiple-testing
#'       correction groups.}
#'     \item{`r`}{Spearman correlation coefficient.}
#'     \item{`n`}{Number of sample pairs used in the correlation (non-missing
#'       observations).}
#'     \item{`P`}{Nominal (unadjusted) p-value.}
#'     \item{`BH`}{FDR-adjusted p-value (Benjamini–Hochberg), corrected within
#'       each `link_label`–`modality_pair` group.}
#'   }
#' @export
correlateByChromosome = function(count_mat, correlation_pairs, grp_contrast = NULL, rds_fn = NULL) {
  chr_pair_ls = correlation_pairs |> select(-any_of(c("ensembl_gene_id"))) |>
    distinct() |>
    (\(x) split(x, x$chr))()
  if (!is.null(grp_contrast)) {
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

#' Annotate whether each link's associated promoter peak is significantly correlated with its gene
#'
#' Adds a `prom_peak_sig_cor_w_gene` column to `df`. The flag indicates
#' whether the promoter peak involved in a given link has a significant
#' `promoter_peak_to_gene` correlation (BH < `BH_thresh`) *within the same
#' input data frame*. This is used in two ways downstream:
#' (1) on the BH-filtered results to identify `distal_peak_to_promoter_peak`
#' links where the promoter anchor is also significantly connected to gene
#' expression, supporting the full distal → promoter → gene chain; and
#' (2) on the broader significant-gene set to identify cases where a
#' distal→promoter link exists but the promoter peak itself is NOT
#' significantly correlated with the gene (`distal_no_promoter` annotation).
#'
#' @param df A correlation results data frame containing columns
#'   `ensembl_gene_id`, `regulatory_element`, `target_id`, `link_label`, and
#'   `BH`.
#' @param BH_thresh FDR threshold used to determine whether a
#'   `promoter_peak_to_gene` link is significant. Default: `0.05`.
#'
#' @return `df` with an additional logical column `prom_peak_sig_cor_w_gene`.
editCorResFields = function(df, BH_thresh = 0.05) {
  sig_promoter_peaks = df |> 
    filter(BH < BH_thresh, link_label == "promoter_peak_to_gene") |> 
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

#' Categorise and annotate significant correlation results
#'
#' Filters to BH-significant links, joins `dist_to_tss` from `p2g_info`, and
#' adds two interpretive columns:
#' \describe{
#'   \item{`prom_peak_sig_cor_w_gene`}{`TRUE` if the link involves a promoter
#'     peak with a significant `promoter_peak_to_gene` correlation, or if the
#'     distal peak is correlated with both the gene and such a promoter peak —
#'     indicating the full distal → promoter → gene chain is supported.}
#'   \item{`distal_no_promoter`}{`TRUE` for distal peaks that are significantly
#'     correlated with both the gene and a proximal promoter peak, but whose
#'     associated promoter peak is itself *not* significantly correlated with
#'     gene expression. These represent distal regulatory candidates that bypass
#'     the expected promoter anchor.}
#' }
#' Column descriptions are also stored in an `"info"` attribute on the
#' returned data frame for reference.
#'
#' @param .sig_cor_genes Data frame of correlation results for genes that have
#'   at least one significant peak-to-gene link (output of the gene-inclusion
#'   filter in [processCorrelations()]). Must contain all link types for those
#'   genes, not just significant ones.
#' @param p2g_info The `p2g_info` data frame from [createPeak2GeneObjects()],
#'   used to join `dist_to_tss` onto the results.
#' @param BH_thresh FDR threshold for defining significance. Default: `0.05`.
#'
#' @return A BH-filtered data frame with `dist_to_tss`, `prom_peak_sig_cor_w_gene`,
#'   and `distal_no_promoter` columns added, and an `"info"` attribute
#'   containing column descriptions.
categorizeAndAnnotateCorrelationResults = function(.sig_cor_genes, p2g_info, BH_thresh = 0.05) {
  sig_cor_res = .sig_cor_genes |> 
    filter(BH < BH_thresh) |> 
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
    filter(((link_label %in% c("distal_peak_to_gene", "distal_peak_to_promoter_peak")) & BH < BH_thresh)) |> 
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

#' Assign direction to ambiguous peaks by correlating with well-anchored distal peaks
#'
#' For peaks whose direction of association with gene expression cannot be
#' determined directly ("disclude" in `region2gene_dir`), attempts to infer a
#' direction by correlating each ambiguous peak against the most significantly
#' gene-linked distal peak for the same gene. If the two peaks are
#' significantly correlated (p < 0.05) the direction is inferred transitively:
#' a positive correlation with a positively-associated distal peak implies a
#' positive association with the gene, and vice versa.
#'
#' One call handles a single iteration — the `n`th-most-significant distal
#' anchor per gene. [corWithDistalPeak_wrapper()] calls this iteratively until
#' no "disclude" peaks remain or the iteration limit is reached.
#'
#' @param df Data frame of significant correlation results with a
#'   `region2gene_dir` column (values: `"pos"`, `"neg"`, `"disclude"`,
#'   `"drop"`).
#' @param full_cor_df The same data frame passed to the first call of this
#'   function — used as a stable reference of significant correlations for
#'   selecting the `n`th-best anchor. Should not be updated between iterations.
#' @param count_mat The combined count matrix from [formatMatrixForCorrelation()].
#' @param iter Integer. Which rank of distal anchor peak to use for this
#'   iteration (1 = most significant). Default: `1`.
#' @param iter_batch_cutoff Integer. Below this iteration the function takes
#'   the single `iter`th peak; at or above it, all remaining candidates are
#'   evaluated together to limit memory use. Increase if the correlation matrix
#'   consumes too much memory. Default: `5`.
#'
#' @return `df` with updated `region2gene_dir` values and an `iter_n` column
#'   recording which iteration resolved each peak.
corWithDistalPeak = function(df, full_cor_df, count_mat, iter = 1, iter_batch_cutoff = 5) {
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

#' Iteratively resolve ambiguous peak–gene direction assignments
#'
#' Calls [corWithDistalPeak()] in a loop, advancing the anchor rank by one
#' each iteration, until no peaks remain labeled `"disclude"` or 20 iterations
#' are reached. Called as Step 3 of [assignRegion2GeneDirection()] to salvage
#' directional assignments for peaks that could not be resolved through the
#' promoter-peak and direct distal-gene correlation routes.
#'
#' @param full_cor_df Data frame with a `region2gene_dir` column, passed
#'   as-is to each [corWithDistalPeak()] call as the stable reference.
#' @param count_mat The combined count matrix from [formatMatrixForCorrelation()].
#' @param v Logical. If `TRUE` (default), prints per-iteration summaries of
#'   `region2gene_dir` value counts.
#' @param .iter_batch_cutoff Passed to [corWithDistalPeak()]. Default: `5`.
#'
#' @return `full_cor_df` with `region2gene_dir` updated and an `iter_n` column
#'   recording which iteration resolved each previously ambiguous peak.
corWithDistalPeak_wrapper = function(full_cor_df, count_mat, v = TRUE, .iter_batch_cutoff = 5) {
  df = full_cor_df
  i = 1
  n_discluded = df |> filter(region2gene_dir == "disclude") |> nrow()
  if (n_discluded == 0) return(df |> filter(region2gene_dir == "disclude") |> mutate(iter_n = NA))
  if (v) message("### Correlating ambiguous peaks with other distal peaks associated with gene ###\npre fxn: ", df$region2gene_dir |> table() |> (\(x) paste0(names(x), ":\t", scales::comma(as.numeric(x)), collapse = "\t"))())
  while ((n_discluded > 0) & i < 20) {
    df = corWithDistalPeak(df, full_cor_df, count_mat, i, iter_batch_cutoff = .iter_batch_cutoff)
    if (v) message("round ", i, ":\t", df$region2gene_dir |> table() |> (\(x) paste0(names(x), ": ", scales::comma(as.numeric(x)), collapse = "\t"))())
    i = i + 1
    n_discluded = df |> filter(region2gene_dir == "disclude") |> nrow()
  }
  return(df)
}

#' Assign directional sign and confidence score to each peak–gene link
#'
#' Consolidates all peak–gene relationships into a consensus direction
#' (`"pos"` or `"neg"`) using a prioritised four-tier approach. The full
#' correlation table (including non-significant pairs) is used so that nominal
#' p-values can inform directionality even when FDR thresholds are not met.
#'
#' **Tier 1 (highest confidence):** Direction is assigned from FDR-significant
#' links. Nominal p-values (< 0.05) are used to break ties when more than one
#' link type contributes directional evidence. All peak–gene directions are in
#' agreement.
#'
#' **Tier 2:** For links with conflicting directional evidence where all
#' contributing links are FDR-significant, the direction from the most
#' FDR-significant link is honored. For example, if the distal peak→gene and
#' distal peak→promoter peak→gene paths disagree, the FDR values of those two
#' link types are compared and the more significant one wins.
#'
#' **Tier 3:** For links with conflicting direction where not all links meet
#' FDR significance, the direction with the most nominally significant
#' p-value (< 0.05) is honored, following the same logic as Tier 2.
#'
#' **Tier 4:** For peaks where the distal peak→promoter peak link is
#' FDR-significant but neither the promoter peak→gene nor distal peak→gene
#' correlations reach nominal significance (p < 0.05), direction is inferred
#' by correlating the ambiguous peak against other peaks already confidently
#' linked to the gene (FDR < 0.05). The most FDR-significant anchors are
#' tried first. The confidence score is `4 + (iter / 1000)`, where `iter` is
#' the rank of the anchor peak that produced a significant (p < 0.05)
#' correlation. Tier 4 links are considered lower confidence.
#'
#' @param sig_cor_res Data frame of annotated significant correlations from
#'   [categorizeAndAnnotateCorrelationResults()].
#' @param all_cor_res The full correlation results data frame (all tested
#'   pairs, not filtered by significance) from [correlateByChromosome()].
#'   Nominal p-values from this table are used for directional inference when
#'   FDR-significant evidence is unavailable.
#' @param count_mat The combined count matrix from [formatMatrixForCorrelation()],
#'   used for the iterative correlation step (Tier 4).
#' @param BH_thresh FDR threshold defining significance. Default: `0.05`.
#'
#' @return A data frame with columns `ensembl_gene_id`, `regulatory_element`,
#'   `region2gene_dir` (`"pos"`, `"neg"`, or `"drop"`), and `region2gene_conf`
#'   (numeric confidence tier as described above).
#' @export
assignRegion2GeneDirection = function(sig_cor_res, all_cor_res, count_mat, BH_thresh = 0.05) {
  #### Step 1: Assign initial relationships based on most accessible information available ####
  full_dir_df = sig_cor_res |> 
    # work with the important info to assign direction
    select(ensembl_gene_id, regulatory_element, target_id, link_label, r, P, BH) |> 
    # add downstream correlation results of promoter peak to gene
    left_join(all_cor_res |> filter(link_label == "promoter_peak_to_gene") |> select(ensembl_gene_id, promoter_peak = regulatory_element, pp_r = r, pp_pval = P, pp_BH = BH), by = c("ensembl_gene_id", "target_id" = "promoter_peak")) |> 
    # add correlation results of regulatory peak itself to the gene. this will be helpful for assigning direction for peaks linked to promoter peaks that are not significantly correlated with gene expression
    left_join(all_cor_res |> filter(link_label %in% c("distal_peak_to_gene", "promoter_peak_to_gene")) |> select(ensembl_gene_id, distal_peak = regulatory_element, dp_r = r, dp_pval = P, dp_BH = BH), by = c("ensembl_gene_id", "regulatory_element" = "distal_peak")) |> 
    mutate(cor_to_compare = case_when(pp_BH < BH_thresh ~ "pp",
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

#' Filter, annotate, and assign directionality to significant peak–gene correlations
#'
#' Orchestrates the post-correlation processing pipeline: retains genes with
#' at least one significant peak-to-gene link, annotates the full regulatory
#' chain context for each link via [categorizeAndAnnotateCorrelationResults()],
#' and assigns a directional sign and confidence tier to every link via
#' [assignRegion2GeneDirection()].
#'
#' @param count_mat The combined count matrix from [formatMatrixForCorrelation()].
#' @param cor_res The full correlation results data frame from
#'   [correlateByChromosome()] (all tested pairs, not pre-filtered).
#' @param correlation_pairs The `correlation_pairs` data frame from
#'   [createPeak2GeneObjects()]. Currently unused in the function body but
#'   retained for potential downstream use.
#' @param p2g_info The `p2g_info` data frame from [createPeak2GeneObjects()],
#'   used to join `dist_to_tss` onto results.
#' @param gene_inclusion_thresh FDR threshold used to determine which genes
#'   are included in further processing: a gene must have at least one
#'   peak-to-gene link (any link ending in `_to_gene`) with BH below this
#'   threshold. Default: `0.05`.
#' @param BH_thresh FDR threshold used for significance within the annotation
#'   and directional assignment steps. Default: `0.05`.
#'
#' @return A data frame filtered to BH-significant links for genes with at
#'   least one significant peak-to-gene correlation. Contains all columns from
#'   [correlateByChromosome()] plus the following:
#'   \describe{
#'     \item{`dist_to_tss`}{Signed distance (bp) from the regulatory element
#'       to the gene TSS, joined from [createPeak2GeneObjects()]. Negative =
#'       upstream, positive = downstream (strand-aware; see
#'       [calculateDirectedDistance()]).}
#'     \item{`prom_peak_sig_cor_w_gene`}{Logical. `TRUE` if the promoter peak
#'       involved in this link has a significant (`BH < BH_thresh`)
#'       `promoter_peak_to_gene` correlation for this gene. Also `TRUE` for
#'       distal peaks significantly correlated with both the gene and such a
#'       promoter peak, indicating the full distal → promoter → gene chain is
#'       supported. Set by [categorizeAndAnnotateCorrelationResults()].}
#'     \item{`distal_no_promoter`}{Logical. `TRUE` for distal peaks that are
#'       significantly correlated with both the gene and a promoter peak, but
#'       whose associated promoter peak is itself not significantly correlated
#'       with gene expression. These links lack a supported promoter anchor.
#'       Set by [categorizeAndAnnotateCorrelationResults()].}
#'     \item{`region2gene_dir`}{Consensus direction of the peak's association
#'       with gene expression: `"pos"` (positive), `"neg"` (negative), or
#'       `"drop"` (direction could not be resolved). Assigned by
#'       [assignRegion2GeneDirection()].}
#'     \item{`region2gene_conf`}{Numeric confidence tier (1–4+) from
#'       [assignRegion2GeneDirection()]. Tiers 1–3 reflect decreasing
#'       evidentiary strength from FDR-significant agreement down to nominal
#'       p-value arbitration. Tier 4+ links required iterative correlation
#'       with other distal peaks to infer direction and are generally excluded
#'       from downstream motif analyses.}
#'   }
#' @export
processCorrelations = function(count_mat, cor_res, correlation_pairs, p2g_info, gene_inclusion_thresh = 0.05, BH_thresh = 0.05) {
  genes_with_sig_anno = cor_res |> 
    filter((link_label %in% grep("to_gene$", levels(link_label), value = TRUE)) & (BH < gene_inclusion_thresh)) |> 
    pull(ensembl_gene_id) |> 
    unique()
  if (length(genes_with_sig_anno) == 0) stop("No significant peak to gene correlations detected.")
  sig_cor = cor_res |> 
    filter(ensembl_gene_id %in% genes_with_sig_anno) |> 
    mutate(sig_cor = BH < BH_thresh) 
  anno_cor_res = categorizeAndAnnotateCorrelationResults(sig_cor, p2g_info, BH_thresh = BH_thresh)
  # return(anno_cor_res)
  region2gene_dir = assignRegion2GeneDirection(anno_cor_res, cor_res, count_mat, BH_thresh = BH_thresh)
  sig_cor_res = anno_cor_res |> 
    left_join(region2gene_dir, by = c("ensembl_gene_id", "regulatory_element")) |> 
    relocate(ensembl_gene_id, regulatory_element, link_label, modality_pair, r, BH, dist_to_tss, starts_with("region2gene")) |> 
    select(-sig_cor)
    # filter(region2gene_dir != "drop")
  return(sig_cor_res)
}
