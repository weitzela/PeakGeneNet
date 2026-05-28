#' @keywords internal
#' @noRd
.customLogBreaks = function(pvals) {
  exponents = unique(round(quantile(-log10(pvals), c(0, 0.25, 0.75, 1))))
  breaks = 10^(-exponents)
  return(breaks)
}

#' @keywords internal
#' @noRd
.addScaleBar = function(data, scale_width, scale_y = 0.8, lwd = 0.25, indent_left_bp = 100) {
  if (scale_width < 1000) {
    scale_label = paste0(scale_width, " bp")
  } else {
    scale_label = paste0(scale_width / 1000, " kb")
  }
  
  # accepts either a dataframe (will search for the maximum x to anchor the scale bar)
  # or a single numeric value (used directly as the right anchor)
  if (inherits(data, "data.frame")) {
    if (!("start" %in% colnames(data))) {
      data = data |>
        select(where(is.numeric))
      start_col = grep("start", colnames(data), value = TRUE)[1]
      scale_start = max(data[[start_col]]) - scale_width - indent_left_bp
    } else {
      scale_start = max(data$start) - scale_width - indent_left_bp
    }
  } else if (inherits(data, "numeric")) {
    scale_start = data - scale_width - indent_left_bp
  }
  
  tick_end_height_npc = 0.02
  list(
    annotate("segment",
             x = c(scale_start, scale_start + scale_width, scale_start),
             xend = c(scale_start, scale_start + scale_width, scale_start + scale_width),
             y = I(c(scale_y + tick_end_height_npc, scale_y + tick_end_height_npc, scale_y)),
             yend = I(c(scale_y - tick_end_height_npc, scale_y - tick_end_height_npc, scale_y)),
             linewidth = lwd),
    annotate("text",
             x = scale_start + (scale_width / 2), y = I(scale_y),
             label = scale_label, vjust = 1.5, size = 5.5, size.unit = "pt")
  )
}

#' Build the long-format link dataframe used by [plotGeneLinks()]
#'
#' One row per peak per bezier control point (`reg_mid`, `link_mid`,
#' `target_mid`). The `y` column encodes the bezier apex height as
#' `-log10(P) * sign(peak_log2FC)` at `link_mid` and 0 at the anchors.
#' Indexes `peak_stats` with the filtered (non-unique) `regulatory_element`
#' vector so that duplicate links — one peak associated with multiple
#' targets — align row-for-row when `bind_cols`ing peak stats onto the
#' filter result. The `slice_min` collapse to one row per peak happens
#' after the join.
#' @keywords internal
#' @noRd
.buildLinkDf = function(pk2g_res, peak_stats, ensembl_to_plot, gene_tss, contrast_labels) {
  filtered = pk2g_res |>
    filter(ensembl_gene_id == ensembl_to_plot, region2gene_dir != "drop") |>
    select(ensembl_gene_id, regulatory_element, target_id, P, region2gene_dir)
  aligned_stats = peak_stats[filtered$regulatory_element]
  filtered |>
    bind_cols(
      data.frame(peak_log2FC = aligned_stats$log2FoldChange, peak_P = aligned_stats$pvalue)
    ) |>
    mutate(peak_dir = ifelse(peak_log2FC > 0, contrast_labels[1], contrast_labels[2])) |>
    separate_wider_delim(regulatory_element, regex("[:punct:]"),
                         names = c("reg_chr", "reg_start", "reg_end", "reg_modality"),
                         too_few = "align_start", cols_remove = FALSE) |>
    mutate(across(ends_with(c("start", "end")), ~ as.numeric(.x))) |>
    mutate(reg_mid = reg_start + ((reg_end - reg_start) / 2),
           target_mid = gene_tss) |>
    slice_min(P, n = 1, with_ties = FALSE, by = "regulatory_element") |>
    mutate(log10P = -log10(P),
           link_mid = (reg_mid + target_mid) / 2,
           region2gene_dir = factor(region2gene_dir, levels = c("pos", "neg"))) |>
    pivot_longer(c(reg_mid, link_mid, target_mid), names_to = "x_sequence", values_to = "x") |>
    mutate(y = ifelse(x_sequence == "link_mid", log10P * sign(peak_log2FC), 0)) |>
    relocate(regulatory_element, x_sequence, x, y)
}

#' Build the peak-triangle polygon dataframe used by [plotGeneLinks()]
#'
#' Apex height is offset per modality (default map: ATACSeq=1, H3K27ac=1.33,
#' H3K4me3=1.66, H3K4me1=2; modalities not in the map default to 1) so
#' overlapping peaks of different modalities remain visually distinguishable
#' even when colored uniformly. `reg_modality` is factored only when
#' `modality_colors` is supplied, using `names(modality_colors)` as the
#' level ordering for the fill scale.
#' @keywords internal
#' @noRd
.buildPeakShapeDf = function(link_df, region_gr, peak_height_npc, peak_width_npc, modality_colors) {
  peak_width_multiplier = round(GenomicRanges::width(range(region_gr)) * peak_width_npc)
  modality_height_default = c("ATACSeq" = 1, "H3K27ac" = 1.33, "H3K4me3" = 1.66, "H3K4me1" = 2)
  
  shape_df = link_df |>
    distinct(regulatory_element, reg_modality, reg_start, reg_end, peak_log2FC) |>
    mutate(reg_mid = (reg_start + reg_end) / 2,
           reg_start = reg_mid - peak_width_multiplier,
           reg_end = reg_mid + peak_width_multiplier) |>
    pivot_longer(c(reg_start, reg_mid, reg_end), names_to = "x_sequence", values_to = "x")
  
  height_lookup = modality_height_default[shape_df$reg_modality]
  height_lookup[is.na(height_lookup)] = 1
  shape_df = shape_df |>
    mutate(y = ifelse(x_sequence == "reg_mid",
                      0.5 + (peak_height_npc * height_lookup * sign(peak_log2FC)),
                      0.5))
  
  if (!is.null(modality_colors)) {
    shape_df = shape_df |>
      mutate(reg_modality = factor(reg_modality,
                                   levels = intersect(names(modality_colors), unique(link_df$reg_modality))))
  }
  shape_df
}

#' Build feature-track tick layers (e.g. SNPs, TF binding sites)
#' @keywords internal
#' @noRd
.buildFeatureTrackLayers = function(feature_track, locus_range, height_npc, color) {
  if (is.null(feature_track)) return(NULL)
  feature_overlap = IRanges::subsetByOverlaps(feature_track, locus_range) |>
    as.data.frame()
  if (nrow(feature_overlap) == 0) return(NULL)
  list(
    geom_segment(data = feature_overlap, inherit.aes = FALSE,
                 aes(x = start, xend = start,
                     y = I(0.5 - height_npc),
                     yend = I(0.5 + height_npc)),
                 linewidth = 0.1, color = "white"),
    geom_segment(data = feature_overlap, inherit.aes = FALSE,
                 aes(x = start, xend = start,
                     y = I(0.5 - height_npc),
                     yend = I(0.5 + height_npc)),
                 linewidth = 0.1, alpha = 0.5, color = color)
  )
}

#' Build the gene-body rectangle layer for genes other than the focal one
#' @keywords internal
#' @noRd
.buildOtherGeneBodiesLayer = function(gene_info, ensembl_to_plot, locus_range, tss_y_npc_adj) {
  xlims = c(GenomicRanges::start(locus_range), GenomicRanges::end(locus_range))
  gene_body_df = IRanges::subsetByOverlaps(gene_info, locus_range) |>
    as.data.frame() |>
    filter(ensembl_gene_id != ensembl_to_plot) |>
    mutate(start = ifelse(start < xlims[1], xlims[1], start),
           end = ifelse(end > xlims[2], xlims[2], end))
  if (nrow(gene_body_df) == 0) return(NULL)
  geom_rect(data = gene_body_df, inherit.aes = FALSE,
            aes(xmin = start, xmax = end,
                ymin = I(0.5 - (tss_y_npc_adj * 0.9)),
                ymax = I(0.5 + (tss_y_npc_adj * 0.9))),
            linewidth = 0, fill = "peachpuff4")
}

#' Build the TSS arrow path dataframe (three points: bottom of stem, top
#' of stem, arrow tip — direction from gene strand)
#' @keywords internal
#' @noRd
.buildTssArrowDf = function(gene_tss, gene_strand, locus_x_range, tss_y_npc_adj) {
  arrow_length_bp = round(diff(locus_x_range) * 0.025)
  arrow_x_end = if (gene_strand == "+") gene_tss + arrow_length_bp else gene_tss - arrow_length_bp
  data.frame(
    x = c(gene_tss, gene_tss, arrow_x_end),
    y = c(0.5 - tss_y_npc_adj, 0.5 + (tss_y_npc_adj * 3), 0.5 + (tss_y_npc_adj * 3)),
    lwd = 0.5,
    group = 1
  )
}

#' Plot peak-to-gene links for a single gene
#'
#' Visualises peak-to-gene links from PeakGeneNet output for a single gene,
#' drawing each link as a bezier curve from the peak to the gene's TSS. The
#' plot is split horizontally around the gene-body midline: peaks with
#' positive log2FC sit above (`contrast_labels[1]`), peaks with negative
#' log2FC sit below (`contrast_labels[2]`). Bezier apex height encodes the
#' peak-to-gene correlation p-value, color encodes correlation direction
#' (`link_dir_colors`), and line thickness encodes the peak's own DESeq2
#' p-value.
#'
#' Optional layers:
#' * `feature_track`: short vertical ticks across the midline marking
#'   arbitrary genomic features (SNPs, TF binding sites, etc.) overlapping
#'   the locus. Pre-filter the input GRanges to control which features appear.
#' * `modality_colors`: per-modality fill coloring of peak triangles, with
#'   per-modality height offsets so overlapping peaks remain distinguishable.
#' * `include_other_gene_bodies`: rectangles for other gene bodies
#'   overlapping the locus.
#'
#' @param ensembl_to_plot Character. Ensembl gene ID to plot. Must be in
#'   `names(gene_info)`.
#' @param pk2g_res Data frame, typically the output of [processCorrelations()].
#'   Must contain columns `ensembl_gene_id`, `regulatory_element`,
#'   `target_id`, `P`, and `region2gene_dir`.
#' @param peak_stats `GRanges` keyed by `regulatory_element` (i.e.
#'   `names(peak_stats)` matches values in `pk2g_res$regulatory_element`,
#'   which follow the `unique_id` convention from [createPeakGr()]:
#'   `"chr:start-end_modality"`). Must have `log2FoldChange` and `pvalue`
#'   metadata columns, e.g. attached from a DESeq2 analysis on the same peaks.
#' @param gene_info `GRanges` keyed by `ensembl_gene_id`, with one entry per
#'   gene whose `start`/`end` define the gene body and `strand` must be `"+"`
#'   or `"-"`. The TSS is derived as the strand-aware 5' end of the gene
#'   body (`start` for `"+"`, `end` for `"-"`). Typically obtained by
#'   subsetting the `gene_gr` output of [createPeak2GeneObjects()] to
#'   `annotation == "GeneBody"`.
#' @param contrast_labels Length-2 character vector. Labels for the two
#'   contrast halves (positive vs negative peak log2FC). Default:
#'   `c("Up", "Down")`.
#' @param contrast_colors Length-2 character vector. Background tint colors
#'   for the two contrast halves. Distinct from `link_dir_colors` because
#'   the two encode different things — peak direction vs link direction.
#'   Default: `c("#E78AC3", "#A6D854")`.
#' @param link_dir_colors Length-2 character vector. Colors for positive
#'   vs negative peak-to-gene correlation links (the bezier curves).
#'   Default: `c("#B2182B", "#2166AC")`.
#' @param modality_colors Optional named character vector mapping modality
#'   to fill color, e.g. `c(ATACSeq = "#abc", H3K27ac = "#def")`. When
#'   non-NULL, peaks are colored by modality and given a modality-specific
#'   vertical offset so overlapping peaks remain visually distinguishable.
#'   When `NULL` (default), all peaks are filled with `peak_color`.
#' @param peak_color Single fill color for peak triangles when
#'   `modality_colors = NULL`. Default: `"black"`.
#' @param peak_alpha Alpha for peak triangle fill. Default: `0.9`.
#' @param peak_height_npc Peak triangle apex height in NPC units. Default:
#'   `0.025`.
#' @param peak_width_npc Peak triangle width as a fraction of the locus
#'   width. Default: `0.01`.
#' @param bkg_lighten Lightening factor passed to [colorspace::lighten()]
#'   for the contrast background tints. Default: `0.8`.
#' @param feature_track Optional `GRanges` of genomic features to mark on
#'   the midline (e.g. SNPs, TF binding sites). Features overlapping the
#'   plotted locus are drawn as short vertical ticks at their start
#'   coordinates. Pre-filter the input GRanges to control which features
#'   appear. Default: `NULL`.
#' @param feature_track_height_npc Half-height of feature ticks in NPC
#'   units. Default: `0.05`.
#' @param feature_track_color Tick color for the feature track. Default:
#'   `"gray25"`.
#' @param include_other_gene_bodies If `TRUE` (default), other gene bodies
#'   overlapping the plotted locus are drawn as small rectangles on the
#'   midline.
#' @param scale_bar_width_bp Width of the scale bar in bp. Default: `2e5`.
#'
#' @return A `ggplot` object with a `"ylim"` attribute giving the symmetric
#'   y-limits used for `coord_cartesian`.
#' @import ggplot2
#' @export
plotGeneLinks = function(ensembl_to_plot,
                         pk2g_res, peak_stats, gene_info,
                         contrast_labels = c("Up", "Down"),
                         contrast_colors = c("#E78AC3", "#A6D854"),
                         link_dir_colors = c("#B2182B", "#2166AC"),
                         modality_colors = NULL,
                         peak_color = "black",
                         peak_alpha = 0.9,
                         peak_height_npc = 0.025,
                         peak_width_npc = 0.01,
                         bkg_lighten = 0.8,
                         feature_track = NULL,
                         feature_track_height_npc = 0.05,
                         feature_track_color = "gray25",
                         include_other_gene_bodies = TRUE,
                         scale_bar_width_bp = 2e5) {
  
  # ---- validate inputs ----
  if (!ensembl_to_plot %in% names(gene_info)) {
    stop("ensembl_to_plot ('", ensembl_to_plot, "') not found in names(gene_info).")
  }
  gene_strand = as.character(GenomicRanges::strand(gene_info[ensembl_to_plot]))
  if (!gene_strand %in% c("+", "-")) {
    stop("gene_info['", ensembl_to_plot, "'] must have strand '+' or '-' (got '",
         gene_strand, "'). The TSS arrow direction cannot be determined without an explicit strand.")
  }
  if (is.null(names(peak_stats))) {
    stop("peak_stats must be named by regulatory_element (set names(peak_stats) <- ...).")
  }
  needed_peaks = pk2g_res |>
    filter(ensembl_gene_id == ensembl_to_plot, region2gene_dir != "drop") |>
    pull(regulatory_element) |>
    unique()
  missing_peaks = setdiff(needed_peaks, names(peak_stats))
  if (length(missing_peaks) > 0) {
    stop("peak_stats is missing entries for ", length(missing_peaks),
         " regulatory element(s): ",
         paste0(head(missing_peaks, 5), collapse = ", "),
         if (length(missing_peaks) > 5) ", ..." else "")
  }
  
  # ---- build dataframes ----
  # .buildLinkDf indexes `peak_stats` with the non-unique filtered regulatory_element
  # vector internally so duplicate links align row-for-row before slice_min collapses
  # them. `region_gr` here is unique-subsetted for `range()`-based extent (collapses anyway).
  gene_tss = if (gene_strand == "+") {
    GenomicRanges::start(gene_info[ensembl_to_plot])
  } else {
    GenomicRanges::end(gene_info[ensembl_to_plot])
  }
  region_gr = peak_stats[needed_peaks]
  link_df = .buildLinkDf(pk2g_res, peak_stats, ensembl_to_plot, gene_tss, contrast_labels)
  peak_shape_df = .buildPeakShapeDf(link_df, region_gr, peak_height_npc, peak_width_npc, modality_colors)
  locus_range = range(region_gr)
  use_modality_colors = !is.null(modality_colors)
  
  # ---- assemble plot: background + contrast labels + feature track + bezier links ----
  p = link_df |>
    ggplot(aes(x = x, y = y * 2,  # *2 so the bezier cubic apex matches log10(p): https://ggforce.data-imaginist.com/reference/geom_bezier.html
               color = region2gene_dir,
               linewidth = peak_P)) +
    annotate("rect",
             xmin = c(-Inf, -Inf), xmax = c(Inf, Inf),
             ymin = c(0, -Inf), ymax = c(Inf, 0),
             fill = colorspace::lighten(contrast_colors, bkg_lighten),
             color = NA, lwd = 0) +
    annotate("text", x = I(rep(.01, 2)), y = I(c(0.99, 0.01)),
             label = contrast_labels, color = contrast_colors, fontface = "bold",
             size = 8, size.unit = "pt", hjust = "inward", vjust = "inward") +
    .buildFeatureTrackLayers(feature_track, locus_range, feature_track_height_npc, feature_track_color) +
    ggforce::geom_bezier(aes(group = regulatory_element), color = "white") +
    ggforce::geom_bezier(aes(group = regulatory_element), alpha = 0.8)
  
  # ---- peak triangles (fill mapping branches on whether modality coloring is on) ----
  if (use_modality_colors) {
    p = p +
      geom_polygon(data = peak_shape_df, inherit.aes = FALSE,
                   aes(x = x, y = I(y)), fill = "white") +
      geom_polygon(data = peak_shape_df, inherit.aes = FALSE,
                   aes(x = x, y = I(y), fill = reg_modality),
                   alpha = peak_alpha, key_glyph = draw_key_point) +
      scale_fill_manual(values = modality_colors)
  } else {
    p = p +
      geom_polygon(data = peak_shape_df, inherit.aes = FALSE,
                   aes(x = x, y = I(y)), fill = "white") +
      geom_polygon(data = peak_shape_df, inherit.aes = FALSE,
                   aes(x = x, y = I(y)),
                   fill = peak_color, alpha = peak_alpha, key_glyph = draw_key_point)
  }
  
  # ---- scales, guides, scale bar, axes, theme ----
  p = p +
    scale_color_manual(values = link_dir_colors, labels = c("Positive", "Negative")) +
    scale_linewidth_continuous(
      range = c(0.05, 1),
      limits = c(10^(-ceiling(max(-log10(link_df$peak_P)))), 1),
      breaks = .customLogBreaks(link_df$peak_P),
      trans = scales::trans_new("neg_log10",
                                transform = function(x) -log10(x),
                                inverse = function(x) 10^(-x)),
      labels = function(x) paste0("×10<sup>", log10(x), "</sup>")
    )
  
  guide_list = list(
    linewidth = guide_legend(title = "Peak DESeq(*p*)", ncol = 2, order = 3, byrow = TRUE, position = "inside"),
    color = guide_legend(title = "Pk2G Link", override.aes = list(shape = NA, lwd = 1, alpha = 1), nrow = 1, order = 1)
  )
  if (use_modality_colors) {
    present_modalities = intersect(names(modality_colors), unique(link_df$reg_modality))
    guide_list$fill = guide_legend(
      title = "Modality",
      override.aes = list(color = modality_colors[present_modalities],
                          alpha = 1, size = 2, shape = 17),
      nrow = 1, order = 2
    )
  }
  p = p +
    do.call(guides, guide_list) +
    .addScaleBar(max(link_df$x), scale_bar_width_bp, 0.95, indent_left_bp = 1e4) +
    scale_x_continuous(labels = scales::label_comma(prefix = paste0(unique(link_df$reg_chr), ":")),
                       n.breaks = 5,
                       expand = expansion(mult = c(0.008, 0.008))) +
    scale_y_continuous(labels = function(x) {
      # y-axis sign comes from sign(peak_log2FC), so ticks belong to the contrast halves
      cols = case_when(x > 0 ~ contrast_colors[1], x < 0 ~ contrast_colors[2], TRUE ~ "black")
      value = ifelse(x == 0, "1", paste0("×10<sup>", str_remove_all(10^(-abs(x)), "1e"), "</sup>"))
      paste0("<span style='color:", cols, "'>", value, "</span>")
    }) +
    ylab("Pk2G Link Significance<br><span style = 'font-size:6pt'>(*p*-value)</span>") +
    geom_hline(yintercept = 0, linewidth = 0.1) +
    theme(axis.title.y = ggtext::element_markdown(),
          axis.text.y = ggtext::element_markdown(),
          axis.line.x = element_blank(),
          legend.key.spacing = unit(0.1, "mm"),
          axis.title.x = element_blank(),
          legend.direction = "vertical",
          legend.text = ggtext::element_markdown(margin = margin(0, 2, 0, 0)),
          legend.title = ggtext::element_markdown(margin = margin()),
          legend.position.inside = c(1, 0),
          legend.justification.inside = c(1, 0),
          legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent"))
  
  # ---- gene-body and TSS layers (added last so coord_cartesian sees the full y-range) ----
  tss_y_npc_adj = peak_height_npc / 4
  
  if (include_other_gene_bodies) {
    p = p + .buildOtherGeneBodiesLayer(gene_info, ensembl_to_plot, locus_range, tss_y_npc_adj)
  }
  
  ymax = max(abs(layer_scales(p)$y$range$range))
  ylims = c(-ymax, ymax)
  
  tss_arrow_df = .buildTssArrowDf(gene_tss, gene_strand, range(link_df$x), tss_y_npc_adj)
  
  p = p +
    geom_path(data = tss_arrow_df, aes(x = x, y = I(y)),
              linewidth = 0.25, inherit.aes = FALSE,
              arrow = arrow(length = unit(0.1, "cm"), ends = "last", type = "closed"),
              position = position_nudge(x = -1),
              linejoin = "round", show.legend = FALSE) +
    annotate("rect",
             xmin = GenomicRanges::start(gene_info[ensembl_to_plot]),
             xmax = GenomicRanges::end(gene_info[ensembl_to_plot]),
             ymin = I(0.5 - tss_y_npc_adj), ymax = I(0.5 + tss_y_npc_adj),
             linewidth = 0, fill = "black") +
    coord_cartesian(ylim = ylims)
  
  p = p |> `attr<-`("ylim", ylims)
  return(p)
}
