#' Signed genomic distance relative to a reference feature
#'
#' Computes the distance between `subject` and `target` genomic ranges and
#' returns a signed value: **negative** if `subject` is upstream of `target`,
#' **positive** if downstream. Strand is taken into account — for genes on the
#' minus strand, upstream corresponds to higher genomic coordinates.
#'
#' This is used in place of `GenomicRanges::distance()`, which returns only
#' unsigned distances, because peak-to-TSS interpretation requires knowing
#' whether a peak lies upstream or downstream. The `dist_to_tss` column in
#' `p2g_info` is populated using this function.
#'
#' @param subject A `GRanges` object representing the feature whose position is
#'   being evaluated (e.g., a peak).
#' @param target A `GRanges` object representing the reference feature
#'   (e.g., a gene TSS). Must be the same length as `subject`.
#'
#' @return A numeric vector of signed distances, the same length as `subject`.
#'   Negative values indicate that the subject is upstream of the target TSS;
#'   positive values indicate downstream.
#' @export
calculateDirectedDistance = function(subject, target) {
  dist = GenomicRanges::distance(subject, target)
  subject_upstream = GenomicRanges::start(subject) < GenomicRanges::start(target)
  dist = ifelse(as.logical(GenomicRanges::strand(target) == "+") & subject_upstream | (as.logical(GenomicRanges::strand(target) == "-") & !subject_upstream), dist * -1, dist)
  return(dist)
}

#' Build a GRanges object of canonical gene TSSs and their promoter windows
#'
#' Queries Ensembl via `biomaRt` for the canonical TSS of each requested gene
#' and returns a `GRanges` object containing two annotation types per gene: the
#' point-width TSS (`annotation == "TSS"`) and the flanking promoter window
#' (`annotation == "Promoter"`). Only autosomes (chromosomes 1–22) are
#' retained. Seqlevels are set to UCSC style (e.g. `"chr1"`).
#'
#' @param ensembl_ids A character vector of Ensembl gene IDs (e.g.
#'   `"ENSRNOG00000000008"`), or a matrix whose row names or column names are
#'   Ensembl gene IDs.
#' @param biomart_ensembl A `Mart` object returned by
#'   `biomaRt::useEnsembl()` pointing to the appropriate organism and Ensembl
#'   version (e.g. `biomaRt::useEnsembl(biomart = "genes", dataset =
#'   "rnorvegicus_gene_ensembl", version = 109)`).
#' @param ucsc_genome UCSC genome assembly name (e.g. `"rn7"`, `"hg38"`).
#'   Used to populate `seqinfo` on the returned object.
#' @param promoter_upstream Number of bases upstream of the TSS to include in
#'   the promoter window. Default: `2000`.
#' @param promoter_downstream Number of bases downstream of the TSS to include
#'   in the promoter window. Default: `1000`.
#'
#' @return A sorted `GRanges` object with one TSS entry and one Promoter entry
#'   per gene. Metadata columns: `ensembl_gene_id`, `annotation` (`"TSS"` or
#'   `"Promoter"`). Names are set to `ensembl_gene_id` for TSS entries.
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

#' Build a GRanges object from a named list of peak regions
#'
#' Parses genomic coordinates from a named list of character vectors or count
#' matrices and returns a single `GRanges` object with modality labels
#' attached. Peak IDs are expected in the format `chr:start:end` or
#' `chr:start-end`; the list name is used as the modality if not encoded in
#' the ID itself.
#'
#' @param peaks A named list where each element is either a character vector of
#'   genomic region strings (e.g. `"chr1:1000-2000"`) or a count matrix whose
#'   row names or column names are genomic regions. List names must be modality
#'   labels (e.g. `"ATACSeq"`, `"H3K27ac"`). **Modality names must match
#'   exactly** — see [createPeak2GeneObjects()] for the set of names recognised
#'   by internal filters.
#' @param ucsc_genome UCSC genome assembly name (e.g. `"rn7"`, `"hg38"`).
#'   Used to populate `seqinfo` on the returned object.
#'
#' @return A sorted `GRanges` object. Metadata columns: `region_id` (original
#'   coordinate string), `modality`, `unique_id` (region_id pasted with
#'   modality, used as a stable identifier throughout the pipeline).
#' @keywords internal
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

#' Set up peak-to-gene link objects for correlation analysis
#'
#' The primary setup function for the PeakGeneNet pipeline. For each gene,
#' peaks within the promoter window are designated as *promoter peaks* and
#' peaks within `locus_radius` of the TSS (outside the promoter window) are
#' designated as *distal peaks*. Four link types are enumerated:
#' `promoter_peak_to_gene`, `distal_peak_to_gene`,
#' `distal_peak_to_promoter_peak`, and `promoter_peak_to_promoter_peak`.
#'
#' The four link types form a topology that ensures every distal regulatory
#' hypothesis is anchored through a proximal element. When
#' `require_promoter_peak = FALSE` this topology is reduced for genes without
#' promoter peaks: only `distal_peak_to_gene` links are generated, and those
#' links lack a promoter-peak anchor — interpret them with caution.
#'
#' **Modality naming:** modality names in `peaks` must match the strings used
#' internally (`"ATACSeq"`, `"H3K27ac"`, `"H3K4me1"`, `"H3K4me3"`) exactly,
#' including capitalisation. Mismatches will silently bypass modality-specific
#' filters (e.g. `promoter_modalities`, `modality_max_dist`).
#'
#' @param genes A character vector of Ensembl gene IDs, or a count matrix
#'   whose row names or column names are Ensembl gene IDs.
#' @param peaks A named list of peak inputs. Each element may be a character
#'   vector of genomic region strings (e.g. `"chr1:1000-2000"`) or a count
#'   matrix with genomic regions as row or column names. List names are used as
#'   modality labels and must match the expected modality name strings exactly.
#' @param biomart_ensembl A `Mart` object from `biomaRt::useEnsembl()` for
#'   the target organism and Ensembl version (e.g.
#'   `biomaRt::useEnsembl(biomart = "genes", dataset =
#'   "rnorvegicus_gene_ensembl", version = 109)`).
#' @param ucsc_genome UCSC genome assembly name (e.g. `"rn7"`, `"hg38"`).
#' @param promoter_upstream Number of bases upstream of the TSS defining the
#'   promoter window. Default: `2000`.
#' @param promoter_downstream Number of bases downstream of the TSS defining
#'   the promoter window. Default: `1000`.
#' @param require_promoter_peak If `TRUE` (default), genes with no overlapping
#'   promoter peaks are excluded entirely. If `FALSE`, such genes are retained
#'   and receive only `distal_peak_to_gene` links; the standard four-link
#'   topology is not available for these genes because there is no proximal
#'   regulatory anchor.
#' @param promoter_modalities Character vector of modality names whose peaks
#'   are eligible to be classified as promoter peaks based on overlap with the
#'   promoter window. Default: `c("H3K27ac", "H3K4me3", "ATACSeq")`.
#' @param modality_max_dist Named numeric vector specifying a hard maximum
#'   absolute distance from the TSS (in bp) for peaks of a given modality.
#'   Peaks beyond this threshold are excluded from `p2g_info` even if they
#'   fall within `locus_radius`. Default: `c(H3K4me3 = 1e4)`. The H3K4me3
#'   default reflects the biology: H3K4me3 is a promoter-restricted histone
#'   mark and distal H3K4me3 signal is generally considered artifactual.
#' @param locus_radius Radius in bp around each TSS defining the distal peak
#'   search window (the promoter region is excluded from this window). Default:
#'   `1e6`.
#' @param verbose If `TRUE` (default), prints a summary of gene and link counts
#'   after construction. If fewer than 50% of queried genes produce any link,
#'   also prints `dist_to_tss` percentiles to help diagnose whether window
#'   parameters need adjustment.
#'
#' @return A named list with four elements and a `stats` attribute:
#'   \describe{
#'     \item{`gene_gr`}{`GRanges` of TSS and promoter window entries for all
#'       queried genes (output of [createTSSGr()]). Used internally for overlap
#'       detection; useful for visualisation.}
#'     \item{`peak_gr`}{`GRanges` of all input peaks across all modalities
#'       (output of [createPeakGr()]). The `promoter_peak` metadata column
#'       indicates whether each peak overlaps any gene's promoter window.}
#'     \item{`p2g_info`}{Data frame linking each peak to every gene whose locus
#'       it falls in. Columns include `unique_id`, `region_id`, `modality`,
#'       `ensembl_gene_id`, `dist_to_tss` (signed, see
#'       [calculateDirectedDistance()]), and `promoter_peak`. This table
#'       carries the distance and promoter-peak metadata used for interpreting
#'       correlation results.}
#'     \item{`correlation_pairs`}{Data frame enumerating every peak–target pair
#'       to be correlated, with columns `ensembl_gene_id`,
#'       `regulatory_element`, `target_id`, `link_label`, `chr`, and
#'       `modality_pair`. This is the primary input to
#'       [correlateByChromosome()].}
#'   }
#'   The `stats` attribute is a list of gene and link count summaries (and
#'   optionally `dist_to_tss` percentiles if the yield was low).
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

#' Adjust a count matrix for covariates while protecting contrasts of interest
#'
#' Fits a linear model across all features simultaneously and subtracts the
#' contribution of technical variables, leaving the biological contrast(s) of
#' interest intact. This is the recommended preprocessing step before
#' computing correlations with [correlateByChromosome()].
#'
#' Correlations in PeakGeneNet are computed on covariate-adjusted matrices,
#' not raw counts. Adjustment removes unwanted variation (batch, sex, etc.)
#' while preserving the contrast under study. VST-transformed counts are the
#' expected input for RNA-seq data; peak modalities should be similarly
#' variance-stabilised before input.
#'
#' @param counts A numeric matrix of transformed counts (VST, log, inverse
#'   rank normalised, etc.) with sample IDs as row names and feature IDs as
#'   column names.
#' @param covariate_df A data frame with all covariates, including both the
#'   contrast variable(s) and variables to be regressed out. Row names must be
#'   sample IDs in the same order as `counts`.
#' @param vars_to_protect A character vector of column names in `covariate_df`
#'   whose coefficients should **not** be subtracted (i.e. the biological
#'   contrasts to preserve). All other variables have their effects removed.
#'   Pass `""` if no variables should be protected.
#' @param return_coefficients Logical. If `TRUE`, the full coefficient matrix
#'   from `lm.fit()` is attached to the result as an attribute named
#'   `"coefficients"`. Default: `FALSE`.
#'
#' @return A numeric matrix of adjusted counts with the same dimensions as
#'   `counts`. If `return_coefficients = TRUE`, the matrix has a
#'   `"coefficients"` attribute containing the fitted coefficients.
#' @export
adjustCovariateMatrix = function(counts, covariate_df, vars_to_protect, return_coefficients = FALSE) {
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
