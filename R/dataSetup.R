#' Genomic distance with the direction 
#' 
#' This function assigns a negative value if the subject is upstream from the target and positive value if the subject is downstream from the target, while taking into account if the gene is located on the negative or positive strand
#' 
#' @return A numeric vector of signed distances
#' @export
#' @examples
#' # example code
#' @keywords internal
calculateDirectedDistance = function(subject, target) {
  dist = GenomicRanges::distance(subject, target)
  subject_upstream = start(subject) < start(target)
  dist = ifelse(as.logical(strand(target) == "+") & subject_upstream | (as.logical(strand(target) == "-") & !subject_upstream), dist * -1, dist)
  return(dist)
}

#' Setup gene TSS object
#' 
#' This creates a dataframe with canonical TSS for each gene that is included in the analysis
#' @examples
#' # createTSSObject(c("ENSR00129", "ENSR00139"), biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 105))
#' @keywords internal
createTSSObject = function(ensembl_ids, biomart_ensembl) {
  # ensembl = biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 105)
  gene_tss = biomaRt::getBM(attributes = c("transcription_start_site", "transcript_is_canonical", "ensembl_gene_id", "chromosome_name", "strand"), filters = "ensembl_gene_id", values = ensembl_ids, mart = biomart_ensembl) |> 
    rename("start" = "transcription_start_site", "chr" = "chromosome_name") |> 
    mutate(end = start) |> 
    drop_na(transcript_is_canonical) |> #remove TSS associated with any non-canonical transcripts
    select(-transcript_is_canonical) |> 
    filter(chr %in% as.character(1:22)) |> 
    mutate(strand = ifelse(strand == 1, "+", "-"),
           annotation = "TSS") |> 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |> 
    GenomicRanges::sort(ignore.strand = TRUE)
  names(gene_tss) = gene_tss$ensembl_gene_id
  genome(gene_tss) = "hg38"
  tmp_chr_info = getChromInfoFromNCBI("GrCh38")
  seqlengths(gene_tss) = tmp_chr_info |> filter(SequenceName %in% as.character(1:22)) |> select(SequenceName, SequenceLength) |> pull(SequenceLength)
  isCircular(gene_tss) = rep(FALSE, 22)
  return(gene_tss)
}