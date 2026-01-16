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
#' @examples
#' # createTSSObject(c("ENSR00129", "ENSR00139"), biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 105), "hg38")
#' @keywords internal
createTSSGr = function(ensembl_ids, biomart_ensembl, ucsc_genome) {
  # ensembl = biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 105)
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
  Seqinfo::seqinfo(gene_gr) = Seqinfo::Seqinfo(genome = ucsc_genome)[as.character(GenomicRanges::seqnames(gene_gr)) |> unique(),]
  
  promoter_gr = gene_gr |> 
    IRanges::promoters(upstream = 2000, downstream = 1000) |> 
    dplyr::mutate(annotation = "Promoter")
  
  gene_gr = plyranges::bind_ranges(gene_gr, promoter_gr) |> 
    GenomicRanges::sort(ignore.strand = TRUE)
  return(gene_gr)
}

createPeakGr = function(peak_counts, ucsc_genome) {
  peak_gr = purrr::imap(peak_counts, function(.x, .y) {
    if (inherits(.x, "matrix")) {
      df = data.frame(region_id = rownames(.x))
    } else {
      df = data.frame(region_id = .x)
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
  Seqinfo::seqinfo(peak_gr) = Seqinfo::Seqinfo(genome = ucsc_genome)[as.character(GenomicRanges::seqnames(peak_gr)) |> unique(),]
  return(peak_gr)
}

#' Set up genomicranges objects for the data included in the correlation analysis
#' @export
#' @examples
#' # createPeak2GeneObjects(gene_counts, peak_counts, biomaRt::useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl", version = 109), "rn7")
#' 
createPeak2GeneObjects = function(gene_counts, peak_counts, biomart_ensembl, ucsc_genome) {
  gene_gr = createTSSGr(rownames(gene_counts), biomart_ensembl, ucsc_genome)
  peak_gr = createPeakGr(peak_counts, ucsc_genome)
  
  promoter_peaks = peak_gr |>
    subset(modality %in% c("H3K27ac", "H3K4me3", "ATACSeq")) |>
    # IRanges::subsetByOverlaps(subset(gene_gr, annotation == "Promoter"))
    plyranges::join_overlap_inner(subset(gene_gr, annotation == "Promoter"))
  
  peak_gr = peak_gr |>
    dplyr::mutate(promoter_peak = unique_id %in% unique(promoter_peaks$unique_id))

  promoter_gr = subset(gene_gr, annotation == "Promoter") |>
    GenomicRanges::split(~ ensembl_gene_id)
  
  # the tss_locus object is the one that spans 1Mb +/- the gene TSS, its promoter regions removed from it
  tss_locus = subset(gene_gr, annotation == "TSS") |>
    BiocGenerics::unstrand() |>
    plyranges::stretch(2e6) |>
    GenomicRanges::trim() # trim to length of chromosome if you have genome info set in the granges object
  tss_locus = GenomicRanges::split(tss_locus, tss_locus$ensembl_gene_id)
  tss_locus = GenomicRanges::setdiff(tss_locus, promoter_gr[names(tss_locus)], ignore.strand = TRUE) |>
    unlist()
  tss_locus$ensembl_gene_id = names(tss_locus)
  names(tss_locus) = NULL
  
  non_proximal_overlaps = plyranges::join_overlap_inner(subset(peak_gr, !promoter_peak), tss_locus)
  
  dist_to_tss = calculateDirectedDistance(
    non_proximal_overlaps,
    subset(gene_gr, annotation == "TSS")[non_proximal_overlaps$ensembl_gene_id]
  )
  
  promoter_peaks_for_paired_df = promoter_peaks |>
    plyranges::remove_names() |>
    GenomicRanges::mcols() |>
    as.data.frame() |>
    dplyr::mutate(promoter_peak = TRUE) |>
    dplyr::select(-annotation) |>
    dplyr::mutate(dist_to_tss = calculateDirectedDistance(peak_gr[unique_id], subset(gene_gr, annotation == "TSS")[ensembl_gene_id]))
  
  paired_df = non_proximal_overlaps |>
    plyranges::remove_names() |>
    GenomicRanges::mcols() |>
    as.data.frame() |>
    dplyr::select(unique_id, region_id, modality, ensembl_gene_id) |>
    mutate(dist_to_tss = dist_to_tss) |>
    dplyr::mutate(promoter_peak = FALSE) |>
    dplyr::bind_rows(promoter_peaks_for_paired_df) |>
    dplyr::arrange(ensembl_gene_id) |>
    # remove H3K4me3 connections that are too far away from the promoter region and would not be expected biologically
    dplyr::filter(!((modality == "H3K4me3") & (abs(dist_to_tss) > 1e4)))
  
  all_pairs = split(paired_df, paired_df$ensembl_gene_id) |>
    purrr::imap(function(.df, .gene) {
      promoter_peaks = dplyr::filter(.df, promoter_peak) |> dplyr::pull(unique_id)
      if (length(promoter_peaks) == 0) return(NULL)
      other_peaks = dplyr::filter(.df, !promoter_peak) |> dplyr::pull(unique_id)
      region_combinations = expand.grid(c(other_peaks, promoter_peaks), promoter_peaks) |>
        dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
        dplyr::filter(!(as.character(Var1) == Var2)) |>
        dplyr::mutate(link_label = ifelse(Var1 %in% promoter_peaks, "promoter_peak_to_promoter_peak", "distal_peak_to_promoter_peak")) |>
        dplyr::bind_rows(data.frame(Var1 = promoter_peaks, Var2 = paste0(.gene, "_RNASeq"), link_label = "promoter_peak_to_gene"))
      if (length(other_peaks) > 0) {
        region_combinations = dplyr::bind_rows(region_combinations, data.frame(Var1 = c(other_peaks), Var2 = paste0(.gene, "_RNASeq"), link_label = "distal_peak_to_gene"))
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
  # all_pairs = all_pairs |>
  #   dplyr::mutate(re_modality = peak_gr[regulatory_element]$modality)
  # dplyr::mutate(re_modality = removeStrAroundCharacter(regulatory_element, "_", "before"),
  #               t_modality = removeStrAroundCharacter(target_id, "_", "before")) |>
  #   dplyr::mutate(modality_pair = ifelse(re_modality < t_modality, paste0(re_modality, "-", t_modality), paste0(t_modality, "-", re_modality)),
  #                 modality_pair = factor(modality_pair)) |>
  #   dplyr::select(-c(re_modality, t_modality))
  
  return(list(gene_gr = gene_gr, peak_gr = peak_gr, paired_df = paired_df, all_pairs = all_pairs))
}
