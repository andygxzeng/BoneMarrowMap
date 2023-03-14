#' KNN Prediction of Hematopoietic Cell Types
#'
#' Uses knnPredict_Seurat to classify celltypes based on nearest neighbours from reference
#' Provides initial estimates for all cells regardless of mapping QC filter
#' Provides final estimates for mapQC passing cells and assigns mapQC failing cells as NA.
#'
#' @param query_obj Query Seurat object as returned by Symphony map_Query
#' @param ref_obj Custom Symphony Reference object
#' @param ref_label Annotation column from Symphony Reference object to use for classification
#' @param k Number of nearest neighbours from Symphony Reference to use for classification
#' @param mapQC_col Column in query object with mapping error QC status for each cell (either "Pass" or "Fail"). Can set as NULL if no filtering desired.
#' @param initial_label Column name to assign initial cell type assignments to, prior to mapQC filtering
#' @param final_label Column name to assign final cell type assignments to, following mapQC filtering
#'
#' @importFrom dplyr %>%
#' @return A Seurat object, with predicted CellType labels stored in the 'final_label' column of the meta.data. All predicted labels including cells that failed QC can be found in 'initial_label'.
#' @export
predict_CellTypes <- function(query_obj, ref_obj, ref_label = 'CellType_Annotation', k = 30, mapQC_col = 'mapping_error_QC',
                              initial_label = 'initial_CellType', final_label = 'predicted_CellType'){

  if (!ref_label %in% colnames(ref_obj$meta_data)) {
    stop('Label \"{ref_label}\" is not available in the reference metadata.')
  }

  ## run KNN-based prediction
  query_obj <- knnPredict_Seurat(
    query_obj = query_obj,
    ref_obj = ref_obj,
    label_transfer = ref_label,
    col_name = initial_label,
    k = k
  )

  # Convert initial celltype to character
  query_obj@meta.data[[initial_label]] <- query_obj@meta.data[[initial_label]] %>% as.character()

  if(is.null(mapQC_col)){
    ## If no mapping error QC, final celltype same as initial.
    query_obj@meta.data[[final_label]] <- query_obj@meta.data[[initial_label]]
    query_obj@meta.data[[paste0(final_label, '_prob')]] <- query_obj@meta.data[[paste0(initial_label, '_prob')]]
  } else if(is.na(mapQC_col)) {
    ## If no mapping error QC, final celltype same as initial.
    query_obj@meta.data[[final_label]] <- query_obj@meta.data[[initial_label]]
    query_obj@meta.data[[paste0(final_label, '_prob')]] <- query_obj@meta.data[[paste0(initial_label, '_prob')]]
  } else {
    ## If mapping error QC, assign annotations of mapQC failing cells as NA.
    query_obj@meta.data[[final_label]] <- ifelse(query_obj@meta.data[[mapQC_col]] %in% c('Fail', 'fail'), NA,
                                                 query_obj@meta.data[[initial_label]])
    query_obj@meta.data[[paste0(final_label, '_prob')]] <- ifelse(query_obj@meta.data[[mapQC_col]] %in% c('Fail', 'fail'), NA,
                                                                  query_obj@meta.data[[paste0(initial_label, '_prob')]])
  }

  return(query_obj)
}
