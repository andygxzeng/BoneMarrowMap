#' Visualize Mapping error cutoff for QC purposes.
#'
#' Low quality cells with low RNA counts and low transcriptional diversity tend to erroneously map to orthochromatic erythroblasts.
#' This will help the user to set a mapping error QC threshold to remove these low quality cells
#'
#' @param query Query Seurat object after symphony mapping and mapping error QC calculation
#' @param mapQC_values Column in query object with mapping error QC status for each cell (either "Pass" or "Fail"). Can set as NULL if no filtering desired.
#' @param mapQC_class Column in query object with mapping error QC status for each cell (either "Pass" or "Fail"). Can set as NULL if no filtering desired.
#'
#' @import ggplot2
#' @importFrom ggpubr theme_pubr
#' @importFrom ggpubr stat_compare_means
#' @importFrom Seurat DimPlot
#' @importFrom Seurat PercentageFeatureSet
#' @return Lost of diagnostic QC plots to evaluate the current mapping error QC cutoff.
#' @export
plot_MappingErrorQC <- function(query, mapQC_values = 'mapping_error_score', mapQC_class = 'mapping_error_QC'){

  if (!mapQC_values %in% colnames(query@meta.data)) {
    stop('Label \"{mapQC_values}\" is not available in the reference metadata. Did you forget to run calculate_MappingError() ?')
  }
  if (!mapQC_class %in% colnames(query@meta.data)) {
    stop('Label \"{mapQC_class}\" is not available in the reference metadata. Did you forget to run calculate_MappingError() ?')
  }

  # Get Percentage Hemoglobin genes, this is to show that low quality cells mis-mapping to orthochromatic erythroblasts do not actually have Hb expression
  query <- Seurat::PercentageFeatureSet(query, features = base::intersect(c('HBB', 'HBA1', 'HBA2', 'HBD', 'HBM'), rownames(query[['RNA']])),
                                        assay = 'RNA', col.name = 'pct_Hb')

  # Plot Mapping error histogram to identify cutoff at tail of distribution
  p1 <- query@meta.data %>%
    ggplot2::ggplot(aes(x = get(mapQC_values), fill = get(mapQC_class))) +
    ggplot2::geom_histogram(bins = 200) + ggpubr::theme_pubr() +
    ggplot2::ylab('Frequency') + ggplot2::xlab(mapQC_values) + labs(fill = mapQC_class) +
    ggplot2::ggtitle('Mapping Error Cutoff')

  # Visualize RNA Counts by cell; MapQC Fail cells should have lower RNA counts
  p2 <- query@meta.data %>%
    ggplot2::ggplot(aes(x = get(mapQC_class), y = nCount_RNA, fill = get(mapQC_class))) +
    ggplot2::geom_boxplot(outlier.size = 0.3) + ggplot2::scale_y_log10() + ggpubr::theme_pubr(legend = 'none') +
    ggpubr::stat_compare_means(comparisons = list(c('Fail', 'Pass')), method = 'wilcox.test') +
    ggplot2::ylab('RNA Count') + ggplot2::xlab(mapQC_class) + labs(fill = mapQC_class) +
    ggplot2::ggtitle('MapQC vs RNA Counts')

  # Visualize mapping error QC on the UMAP
  p3 <- Seurat::DimPlot(query, reduction = 'umap_projected', group.by = c(mapQC_class), raster=FALSE)

  # Visualize percent Hemoglobin genes, true Poly/Orthochromatic Erythroblasts will have >50% Hb
  p4 <- query@meta.data %>%
    ggplot2::ggplot(aes(x = get(mapQC_class), y = pct_Hb, fill = get(mapQC_class))) +
    ggplot2::geom_boxplot(outlier.size = 0.3) + ggpubr::theme_pubr(legend = 'none') +
    ggplot2::ylab('Percent Hemoglobin') + ggplot2::xlab(mapQC_class) + ylim(c(0,100)) + labs(fill = mapQC_class) +
    ggplot2::ggtitle('MapQC vs % Hemoglobin')

  # Return as list of plots
  mapQC_plots <- list('MapQC_Histogram' = p1, 'MapQC_RNAcounts' = p2, 'MapQC_UMAP' = p3, 'MapQC_Hemoglobin' = p4)
  return(mapQC_plots)
}
