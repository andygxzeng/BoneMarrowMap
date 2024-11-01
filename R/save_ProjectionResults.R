#' Function to save CellType and Pseudotime predictions and UMAP coordinates
#'
#' @param query_obj Query Seurat object after symphony mapping and celltype classification
#' @param file_name File name of csv to store projection annotations and UMAP coordinates
#' @param celltype_label Column in query metadata with projected CellType classification.
#' @param celltype_KNNprob_label Column in query metadata with KNN probabilities of the projected CellType classification.
#' @param pseudotime_label Column in query metadata with predicted Pseudotime values
#' @param save_AUCell_scores Boolean. Whether to save AUCell scores along with projection results. If TRUE, looks for metadata columns with suffix '_AUC' and saves those columns. This only works if AUCell scores are saved to the metadata. If you have them saved as an assay within the seurat object, please just save the seurat object with saveRDS.
#'
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom utils write.csv
#' @return Saves a csv file with projected annotations to file_name
#' @export
save_ProjectionResults = function(query_obj, file_name, celltype_label = 'predicted_CellType',
                                  celltype_KNNprob_label = 'predicted_CellType_prob', pseudotime_label = 'predicted_Pseudotime',
                                  save_AUCell_scores = FALSE){

  if(is.null(celltype_label)){
    warning('Label \"{celltype_label}\" set to NULL. This column will not be saved to annotation file')
  } else if (!celltype_label %in% colnames(query_obj@meta.data)) {
      stop('Label \"{celltype_label}\" is not available in the query metadata and is not NULL. Did you forget to run predict_CellTypes?
           \nIf you want to exclude this annotation from the output, set celltype_label to NULL')
  }

  if(is.null(pseudotime_label)){
    warning('Label \"{pseudotime_label}\" set to NULL. This column will not be saved to annotation file')
  } else if (!pseudotime_label %in% colnames(query_obj@meta.data)) {
      stop('Label \"{pseudotime_label}\" is not available in the query metadata and is not NULL. Did you forget to run predict_CellTypes?
           \nIf you want to exclude this annotation from the output, set pseudotime_label to NULL')
  }

  # Get projection results
  projection_results <- dplyr::bind_cols(
    # CellType and Pseudotime labels from metadata
    query_obj@meta.data %>% dplyr::select(mapping_error_score, mapping_error_QC),
    # umap coordinates
    data.frame(query_obj@reductions$umap@cell.embeddings) %>%
      dplyr::rename(UMAP1_projected = umap_1, UMAP2_projected = umap_2)
  )


  # Add celltype annotations if not null
  if(!is.null(celltype_label)){
    if (!celltype_KNNprob_label %in% colnames(query_obj@meta.data)) {
        warning('KNN prob for label \"{celltype_label}\" not found in the query metadata. This will not be saved to the annotation file.')
        projection_results <- projection_results %>% dplyr::bind_cols(query_obj@meta.data %>% dplyr::select(contains(celltype_label)))
    } else {
      projection_results <- projection_results %>%
        dplyr::bind_cols(query_obj@meta.data %>% dplyr::select(contains(celltype_label), all_of(celltype_KNNprob_label)))
    }
  }
  # Add pseudotime annotations if not null
  if(!is.null(pseudotime_label)){
    projection_results <- projection_results %>% dplyr::bind_cols(query_obj@meta.data %>% dplyr::select(all_of(pseudotime_label)))
  }

  # Save AUCell Scores if specified
  if(save_AUCell_scores == TRUE){
    projection_results <- projection_results %>% dplyr::bind_cols(query_obj@meta.data %>% dplyr::select(contains('_AUC')))
  }

  # Save as csv
  projection_results %>% tibble::rownames_to_column('Cell') %>% utils::write.csv(file_name, row.names = FALSE)
}


