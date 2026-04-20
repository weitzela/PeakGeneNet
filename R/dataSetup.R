#' Genomic distance with the direction 
#' 
#' This function assigns a negative value if the subject is upstream from the target and positive value if the subject is downstream from the target, while taking into account if the gene is located on the negative or positive strand. This function is created because the GenomicRanges::distance function itself will not account for direction from TSS (upstream (-) or downstream (+)).
#' 
#' @return A numeric vector of signed distances
#' @export
#' @examples
#' # example code
#' @keywords internal
calculateDirectedDistance = function(subject, target) {
  dist = GenomicRanges::distance(subject, target)
  subject_upstream = GenomicRanges::start(subject) < GenomicRanges::start(target)
  dist = ifelse(as.logical(GenomicRanges::strand(target) == "+") & subject_upstream | (as.logical(GenomicRanges::strand(target) == "-") & !subject_upstream), dist * -1, dist)
  return(dist)
}

#' Setup gene TSS object
#' 
#' This creates a dataframe with canonical TSS for each gene that is included in the analysis
#' @param ucsc_genome genome release formatting in ucsc style. For example: hg38, rn7, etc.
#' @param promoter_upstream number of bases upstream of the TSS to include in the promoter window. Default: 2000.
#' @param promoter_downstream number of bases downstream of the TSS to include in the promoter window. Default: 1000.
#' @examples
#' # createTSSObject(c("ENSR00129", "ENSR00139"), biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 105), "hg38")
#' @keywords internal
createTSSGr = function(ensembl_ids, biomart_ensembl, ucsc_genome,
                       promoter_upstream = 2000, promoter_downstream = 1000) {
  if (inherits(ensembl_ids, "matrix")) {
    # if provided a count matrix, detect whether the IDs of the count matrix are row names or column names
    if (all(grepl("^ENS", rownames(ensembl_ids)))) {
      ensembl_ids = rownames(ensembl_ids)
    } else if (all(grepl("^ENS", colnames(ensembl_ids)))) {
      ensembl_ids = colnames(ensembl_ids)
    }
  } 
  stopifnot(inherits(ensembl_ids, "character"))
  ensembl_ids = unique(ensembl_ids)
    
  gene_gr = biomaRt::getBM(attributes = c("transcription_start_site", "transcript_is_canonical", "ensembl_gene_id", "chromosome_name", "strand"), filters = "ensembl_gene_id", values = ensembl_ids, mart = biomart_ensembl) |> 
    dplyr::rename("start" = "transcription_start_site", "chr" = "chromosome_name") |> 
    dplyr::mutate(end = start) |> 
    tidyr::drop_na(transcript_is_canonical) |> #remove TSS associated with any non-canonical transcripts
    dplyr::select(-transcript_is_canonical) |> 
    dplyr::filter(chr %in% as.character(1:22)) |> 
    dplyr::mutate(strand = ifelse(strand == 1, "+", "-"),
           annotation = "TSS") |> 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |> 
    GenomicRanges::sort(ignore.strand = TRUE) |> 
    GenomeInfoDb::`seqlevelsStyle<-`("UCSC")
  names(gene_gr) = gene_gr$ensembl_gene_id
  GenomeInfoDb::genome(gene_gr) = ucsc_genome
  GenomeInfoDb::seqinfo(gene_gr) = GenomeInfoDb::Seqinfo(genome = ucsc_genome)[as.character(GenomicRanges::seqnames(gene_gr)) |> unique(),]
  
  promoter_gr = gene_gr |> 
    IRanges::promoters(upstream = promoter_upstream, downstream = promoter_downstream)
  promoter_gr$annotation = "Promoter"
  
  gene_gr = c(gene_gr, promoter_gr) |> 
    GenomicRanges::sort(ignore.strand = TRUE)
  return(gene_gr)
}

createPeakGr = function(peaks, ucsc_genome) {
  if (!inherits(peaks, "list")) {
    stop("Peak regions must be input within a list. E.g., list(ATACSeq = c('chr1:1-2', 'chr1:4-7'), H3K4me1 = c('chr1:2-5', 'chr1:6-8'))")
  }
  peak_gr = purrr::imap(peaks, function(.x, .y) {
    if (inherits(.x, "matrix")) {
      if (all(grepl("^chr", rownames(.x), ignore.case = TRUE))) {
        df = data.frame(region_id = rownames(.x))
      }
      if (all(grepl("^chr", colnames(.x), ignore.case = TRUE))) {
        df = data.frame(region_id = colnames(.x))
      }
    } else if (inherits(.x, "character")) {
      df = data.frame(region_id = .x)
    } else {
      stop("Genomic regions not detected as character vector or matrix rows or column names.")
    }
    df = df |> 
      tidyr::separate_wider_delim(region_id, delim = stringr::regex("[:punct:]"), names = c("chr", "start", "end", "modality"), too_few = "align_start", cols_remove = FALSE) |> 
      dplyr::mutate(modality = ifelse(is.na(modality), .y, modality))
  }) |> 
    dplyr::bind_rows() |> 
    dplyr::mutate(unique_id = paste0(region_id, "_", modality)) |> 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |> 
    GenomeInfoDb::`seqlevelsStyle<-`("UCSC") |> 
    GenomicRanges::sort(ignore.strand = TRUE)
  names(peak_gr) = peak_gr$unique_id
  GenomeInfoDb::genome(peak_gr) = ucsc_genome
  GenomeInfoDb::seqinfo(peak_gr) = GenomeInfoDb::Seqinfo(genome = ucsc_genome)[as.character(GenomicRanges::seqnames(peak_gr)) |> unique(),]
  return(peak_gr)
}

#' Create Peak-Gene Links
#' Set up genomicranges objects for the data included in the correlation analysis. Data matrices must have *unique* feature IDs.
#' @param genes either a character vector of ensembl IDs or a matrix with the row or column names as the ensembl IDs
#' @param peaks a named list of character vectors containing renomic region. Genomic regions should begin with "chr". The list can also contain count matrices, where row names or column names are genomic regions.
#' @param promoter_upstream number of bases upstream of the TSS defining the promoter window. Default: 2000.
#' @param promoter_downstream number of bases downstream of the TSS defining the promoter window. Default: 1000.
#' @param require_promoter_peak if TRUE (default), genes with no overlapping promoter peaks are excluded. If FALSE, genes with no promoter peaks are retained as long as they have at least one distal peak; only distal_peak_to_gene links are generated for such genes.
#' @param promoter_modalities character vector of modality names whose peaks are eligible to be classified as promoter peaks based on overlap with the promoter window. Default: c("H3K27ac", "H3K4me3", "ATACSeq").
#' @param modality_max_dist named numeric vector specifying a maximum absolute distance from the TSS (in bp) for peaks of a given modality. Peaks exceeding the threshold for their modality are excluded from p2g_info. Names must match modality names in the peaks input. Default: c(H3K4me3 = 1e4).
#' @param locus_radius radius in bp around each TSS defining the distal peak search window (promoter region is excluded from this window). Default: 1e6.
#' @param verbose if TRUE (default), prints a summary of link and gene counts after construction. If fewer than 30% of queried genes produce any link, also prints the dist_to_tss percentile distribution to help diagnose window width.
#' @export
#' @examples
#' # genes = c("ENSRNOG00000000008", "ENSRNOG00000000082", "ENSRNOG00000001489")
#' # peaks = list(
#'      ATACSeq = c("chr2:101609880:101610274", "chr2:102047604:102048070", "chr3:79872103:79872545"),
#'      H3K4me3 = c("chr2:102549207:102550717", "chr3:79617231:79617692", "chr3:79610922:79613046")
#'      )
#' # createPeak2GeneObjects(genes, peaks, biomaRt::useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl", version = 109), "rn7")
#' 
createPeak2GeneObjects = function(genes, peaks, biomart_ensembl, ucsc_genome,
                                  promoter_upstream = 2000, promoter_downstream = 1000,
                                  require_promoter_peak = TRUE,
                                  promoter_modalities = c("H3K27ac", "H3K4me3", "ATACSeq"),
                                  modality_max_dist = c(H3K4me3 = 1e4),
                                  locus_radius = 1e6,
                                  verbose = TRUE) {
  gene_gr = createTSSGr(
    genes, biomart_ensembl, ucsc_genome,
    promoter_upstream = promoter_upstream, promoter_downstream = promoter_downstream
  )
  peak_gr = createPeakGr(peaks, ucsc_genome)
  
  promoter_options = subset(peak_gr, modality %in% promoter_modalities)
  promoter_locus = subset(gene_gr, annotation == "Promoter")
  promoter_olaps = GenomicRanges::findOverlaps(promoter_options, promoter_locus)
  promoter_peaks = promoter_options[promoter_olaps@from]
  promoter_peaks$ensembl_gene_id = promoter_locus$ensembl_gene_id[promoter_olaps@to]
  
  peak_gr$promoter_peak = peak_gr$unique_id %in% unique(promoter_peaks$unique_id)

  promoter_gr = subset(gene_gr, annotation == "Promoter") |>
    GenomicRanges::split(~ ensembl_gene_id)
  
  # tss_locus spans locus_radius +/- the TSS, with the promoter region removed
  tss_locus = subset(gene_gr, annotation == "TSS") |>
    BiocGenerics::unstrand()
  suppressWarnings(tss_locus <- tss_locus + locus_radius)
  tss_locus = GenomicRanges::trim(tss_locus) # trim to length of chromosome if you have genome info set in the granges object
  tss_locus = GenomicRanges::split(tss_locus, tss_locus$ensembl_gene_id)
  tss_locus = GenomicRanges::setdiff(tss_locus, promoter_gr[names(tss_locus)], ignore.strand = TRUE) |>
    unlist()
  tss_locus$ensembl_gene_id = names(tss_locus)
  names(tss_locus) = NULL
  
  distal_options = subset(peak_gr, !promoter_peak)
  distal_olaps = GenomicRanges::findOverlaps(distal_options, tss_locus)
  non_proximal_overlaps = distal_options[distal_olaps@from]
  non_proximal_overlaps$ensembl_gene_id = tss_locus$ensembl_gene_id[distal_olaps@to]
  
  dist_to_tss = calculateDirectedDistance(
    non_proximal_overlaps,
    subset(gene_gr, annotation == "TSS")[non_proximal_overlaps$ensembl_gene_id]
  )
  
  promoter_peaks_for_p2g_info = promoter_peaks |>
    `names<-`(NULL) |> 
    GenomicRanges::mcols() |>
    as.data.frame() |>
    dplyr::mutate(promoter_peak = TRUE) |>
    dplyr::mutate(dist_to_tss = calculateDirectedDistance(peak_gr[unique_id], subset(gene_gr, annotation == "TSS")[ensembl_gene_id]))
  
  p2g_info = non_proximal_overlaps |>
    `names<-`(NULL) |> 
    GenomicRanges::mcols() |>
    as.data.frame() |>
    dplyr::select(unique_id, region_id, modality, ensembl_gene_id) |>
    mutate(dist_to_tss = dist_to_tss) |>
    dplyr::mutate(promoter_peak = FALSE) |>
    dplyr::bind_rows(promoter_peaks_for_p2g_info) |>
    dplyr::arrange(ensembl_gene_id) |>
    dplyr::filter(is.na(modality_max_dist[modality]) | abs(dist_to_tss) <= modality_max_dist[modality])
  
  correlation_pairs = split(p2g_info, p2g_info$ensembl_gene_id) |>
    purrr::imap(function(.df, .gene) {
      promoter_peaks = dplyr::filter(.df, promoter_peak) |> dplyr::pull(unique_id)
      other_peaks = dplyr::filter(.df, !promoter_peak) |> dplyr::pull(unique_id)
      if (require_promoter_peak && length(promoter_peaks) == 0) return(NULL)
      if (length(promoter_peaks) == 0 && length(other_peaks) == 0) return(NULL)
      region_combinations = data.frame()
      if (length(promoter_peaks) > 0) {
        region_combinations = expand.grid(
          c(other_peaks, promoter_peaks), promoter_peaks
        ) |>
          dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
          dplyr::filter(!(as.character(Var1) == Var2)) |>
          dplyr::mutate(link_label = ifelse(
            Var1 %in% promoter_peaks,
            "promoter_peak_to_promoter_peak",
            "distal_peak_to_promoter_peak"
          )) |>
          dplyr::bind_rows(data.frame(
            Var1 = promoter_peaks,
            Var2 = paste0(.gene, "_RNASeq"),
            link_label = "promoter_peak_to_gene"
          ))
      }
      if (length(other_peaks) > 0) {
        region_combinations = dplyr::bind_rows(
          region_combinations,
          data.frame(
            Var1 = other_peaks,
            Var2 = paste0(.gene, "_RNASeq"),
            link_label = "distal_peak_to_gene"
          )
        )
      }
      region_combinations = region_combinations |>
        dplyr::rename("regulatory_element" = "Var1", "target_id" = "Var2")
      return(region_combinations)
    }) |>
    purrr::compact() |>
    dplyr::bind_rows(.id = "ensembl_gene_id") |>
    dplyr::mutate(chr = GenomeInfoDb::seqnames(subset(gene_gr, annotation == "TSS")[ensembl_gene_id]) |> as.character()) |>
    dplyr::mutate(chr = factor(chr, levels = gtools::mixedsort(unique(chr)))) |>
    dplyr::mutate(link_label = factor(link_label, levels = c("promoter_peak_to_gene", "distal_peak_to_gene", "distal_peak_to_promoter_peak", "promoter_peak_to_promoter_peak")))
  correlation_pairs = correlation_pairs |>
    left_join(p2g_info |> distinct(unique_id, re_modality = modality), by = c("regulatory_element" = "unique_id")) |> 
    left_join(p2g_info |> distinct(unique_id, t_modality = modality), by = c("target_id" = "unique_id")) |> 
    mutate(across(ends_with("modality"), ~ ifelse(is.na(.x), "RNASeq", .x))) |> 
    dplyr::mutate(modality_pair = ifelse(re_modality < t_modality, paste0(re_modality, "-", t_modality), paste0(t_modality, "-", re_modality)),
                  modality_pair = factor(modality_pair)) |>
    dplyr::select(-c(re_modality, t_modality))
  
  n_queried = length(unique(subset(gene_gr, annotation == "TSS")$ensembl_gene_id))
  n_any_peak = length(unique(p2g_info$ensembl_gene_id))
  n_promoter_peak = length(unique(
    dplyr::filter(p2g_info, promoter_peak)$ensembl_gene_id
  ))
  n_linked = length(unique(correlation_pairs$ensembl_gene_id))

  links_by_label = table(correlation_pairs$link_label)
  links_by_modality = table(correlation_pairs$modality_pair)

  stats = list(
    n_genes_queried        = n_queried,
    n_genes_any_peak       = n_any_peak,
    n_genes_promoter_peak  = n_promoter_peak,
    n_genes_linked         = n_linked,
    links_by_label         = links_by_label,
    links_by_modality_pair = links_by_modality
  )

  if (verbose) {
    message(sprintf(
      "Genes queried: %d | with any peak in locus: %d | with promoter peak: %d | linked: %d",
      n_queried, n_any_peak, n_promoter_peak, n_linked
    ))
    message("Links by type:")
    message(paste(sprintf("  %-40s %d", names(links_by_label), as.integer(links_by_label)), collapse = "\n"))
    message("Links by modality pair:")
    message(paste(sprintf("  %-40s %d", names(links_by_modality), as.integer(links_by_modality)), collapse = "\n"))

    if (n_linked / n_queried < 0.5) {
      message(sprintf(
        "\nWarning: only %.0f%% of queried genes produced links. Consider adjusting promoter_upstream/promoter_downstream/locus_radius.",
        100 * n_linked / n_queried
      ))
      message("dist_to_tss percentiles across all peaks in p2g_info:")
      pcts = quantile(abs(p2g_info$dist_to_tss), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
      message(paste(sprintf("  p%-2s: %d bp", names(pcts), as.integer(pcts)), collapse = "\n"))
      stats$dist_to_tss_percentiles = pcts
    }
  }

  result = list(
    gene_gr = gene_gr,
    peak_gr = peak_gr,
    p2g_info = p2g_info,
    correlation_pairs = correlation_pairs
  )
  attr(result, "stats") = stats
  return(result)
}

#' Use of matrix functions to adjust data via linear model
#' 
#' This is a wrapper function that adjusts for multiple additional experimental variables while protecting the contrast(s) of interest. 
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
