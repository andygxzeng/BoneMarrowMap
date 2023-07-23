#' Density plots for bone marrow projection results
#'
#' @param query_obj Query Seurat object as returned by Symphony map_Query
#' @param batch_key Batch key specifying condition to subset and display samples by
#' @param sample_name Specific sample to display density plot for
#' @param ref_obj Bone Marrow Reference object to pull background UMAP plot from
#' @param Hierarchy_only Whether to only show the Hematopoietic hierarchy and exclude T/NK/Plasma/Stromal cells from the plot. Default is FALSE
#' @param downsample_reference Whether to downsample the BM reference to reduce file size. Default is TRUE
#' @param downsample_frac If downsampling reference, what proportion of cells to sample from the BM reference. Default is 0.25
#' @param query_point_size Point size of individual query cells on plot. Default is 0.2
#' @param query_contour_size Contour size depicting query cell density on plot. Default is 0.3
#' @param saveplot Whether to save the plot. Default is TRUE
#' @param device Whether to save the plot as a pdf or png. Default is pdf
#' @param save_folder Folder to save plots in. If this folder does not exist, it will be created
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom ggpointdensity geom_pointdensity
#' @importFrom jcolors scale_color_jcolors_contin
#' @importFrom grDevices heat.colors
#' @return A ggplot object with projected coordinates of the query cells overlaid over a background of the reference cells.
#' @export
plot_Projection_BM <- function(query_obj, batch_key, sample_name, ref_obj, Hierarchy_only = FALSE, downsample_reference = TRUE, downsample_frac = 0.25,
                               query_point_size = 0.2, query_contour_size = 0.3, saveplot = TRUE, device = 'pdf', save_folder = 'projectionFigures/'){

  if (!batch_key %in% colnames(query_obj@meta.data)) {
    stop('Label \"{batch_key}\" is not available in the query metadata.')
  }
  if (!sample_name %in% unique(query_obj@meta.data[[batch_key]])) {
    stop('Label \"{sample_name}\" not found in the \"{batch_key}\" variable within the query metadata')
  }

  # Prepare query_obj data for that sample
  dat <- query_obj@meta.data %>% tibble::rownames_to_column('Cell') %>%
    dplyr::filter(.data[[batch_key]] == sample_name) %>%
    dplyr::left_join(query_obj@reductions$umap@cell.embeddings %>% data.frame() %>% tibble::rownames_to_column('Cell'), by = 'Cell') %>%
    dplyr::filter(mapping_error_QC == 'Pass') %>%
    dplyr::mutate(ref_query = 'query')
  if(Hierarchy_only){
    dat <- dat %>% dplyr::filter(umap_1 > -5)
  }

  # Prepare background from Reference UMAP
  background <- data.frame(ref_obj$umap$embedding) %>% dplyr::rename(umap_1 = X1, umap_2 = X2)
  if(Hierarchy_only){
    background <- background %>% dplyr::filter(umap_1 > -5)
  }
  if(downsample_reference){
    set.seed(123)
    background <- background %>% dplyr::sample_frac(downsample_frac)
  }

  ## Density plot overlaid on the background
  heatpalette <- grDevices::heat.colors(12)
  p <- dat %>%
    ggplot2::ggplot(aes(x = umap_1, y = umap_2)) +
    ggplot2::geom_point(data = background, color='#E3E3E3', size=0.05, alpha=0.5) +
    ggpointdensity::geom_pointdensity(size=query_point_size) +
    jcolors::scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 1.75) +
    ggplot2::geom_density_2d(alpha=0.4, color='black', h = 1.5, size=query_contour_size) +
    ggplot2::theme_void() + ggplot2::ggtitle(sample_name) +
    ggplot2::theme(strip.text.x = ggplot2::element_text(size=18), legend.position='none')

  if(saveplot){
    if(!file.exists(save_folder)){
      dir.create(file.path(paste0('./', save_folder)))
    }
    if(Hierarchy_only){
      if(device %in% c('pdf', 'PDF', 'Pdf')){
        ggsave(paste0(save_folder, 'density_', sample_name, '_projectedUMAP_HierarchyOnly.pdf'), height = 4, width = 4.5, device = 'pdf')
      } else if(device %in% c('png', 'PNG', 'Png')){
        ggsave(paste0(save_folder, 'density_', sample_name, '_projectedUMAP_HierarchyOnly.png'), height = 4, width = 4.5, device = 'png', dpi = 240)
      } else {
        stop('Please specify device as either a pdf or png.')
      }
    } else {
      if(device %in% c('pdf', 'PDF', 'Pdf')){
        ggsave(paste0(save_folder, 'density_', sample_name, '_projectedUMAP.pdf'), height = 4, width = 6, device = 'pdf')
      } else if(device %in% c('png', 'PNG', 'Png')){
        ggsave(paste0(save_folder, 'density_', sample_name, '_projectedUMAP.pdf'), height = 4, width = 6, device = 'png', dpi = 240)
      } else {
        stop('Please specify device as either a pdf or png.')
      }
    }
  }
  return(p)
}
