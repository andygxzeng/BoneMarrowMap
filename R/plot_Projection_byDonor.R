#' Density plots for bone marrow projection results
#'
#' @param query_obj Query Seurat object as returned by Symphony map_Query
#' @param batch_key Batch key specifying condition to subset and display samples by
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
#' @return A ggplot object with projected coordinates of the query cells overlaid over a background of the reference cells.
#' @export
plot_Projection_byDonor <- function(query_obj, batch_key, ref_obj, Hierarchy_only = FALSE, downsample_reference = TRUE, downsample_frac = 0.25,
                                    query_point_size = 0.2, query_contour_size = 0.3, saveplot = TRUE, device = 'pdf', save_folder = 'projectionFigures/'){

  if (!batch_key %in% colnames(query_obj@meta.data)) {
    stop('Label \"{batch_key}\" is not available in the query metadata.')
  }

  # Loop through donors / variables from batch key
  for(donor in unique(query_obj@meta.data[[batch_key]]) ){

    # add plot from each donor into list
    donor_plots[[donor]] <-
      plot_Projection_BM(
        query_obj = query_obj,
        batch_key = batch_key,
        sample_name = donor,
        ref_obj = ref_obj,
        Hierarchy_only = Hierarchy_only, # Whether to exclude T/NK/Plasma/Stromal cells
        downsample_reference = downsample_reference,
        downsample_frac = downsample_frac,   # down-sample reference cells to 25%; reduces figure file size
        query_point_size = query_point_size,   # adjust size of query cells based on # of cells
        query_contour_size = query_contour_size,
        saveplot = saveplot,
        device = device,
        save_folder = save_folder)
  }

  return(donor_plots)
}
