#' Annotate Hematopoietic Hierarchy Pseudotime metric
#'
#' Assigns each query cell the median pseudotime of k-nearest neighbours by UMAP from the reference map.
#' Provides initial estimates for all cells regardless of mapping QC filter
#' Provides final estimates for mapQC passing cells and assigns mapQC failing cells as NA.
#'
#' @param query_obj Query Seurat object as returned by Symphony mapQuery
#' @param ref_obj Custom Symphony Reference object
#' @param pseudotime_label Pseudotime annotation column from Symphony Reference object
#' @param k Number of nearest neighbours from Symphony Reference to use for pseudotime annotation
#' @param mapQC_class Column in query object with mapping error QC status for each cell (either "Pass" or "Fail"). Can set as NULL if no filtering desired.
#' @param initial_label Column name to assign initial pseudotime annotations to, prior to mapQC filtering
#' @param final_label Column name to assign final pseudotime annotations to, following mapQC filtering
#'
#' @importFrom dplyr %>%
#' @importFrom RANN nn2
#' @importFrom Seurat Embeddings
#' @return A Seurat object, with annotated pseudotime labels stored in the 'final_label' column of the meta.data.
#' @export
#'
predict_Pseudotime <- function(query_obj, ref_obj, pseudotime_label = 'Pseudotime', k = 30, mapQC_class = 'mapping_error_QC',
                               initial_label = 'initial_Pseudotime', final_label = 'predicted_Pseudotime'){

  if (!pseudotime_label %in% colnames(ref_obj$meta_data)) {
    stop('Label \"{pseudotime_label}\" is not available in the reference metadata.')
  }

  # get pseudotime values from reference
  ref_pseudotime <- ref_obj$meta_data[pseudotime_label] %>% data.matrix()

  # Assign pseudotime as median of k-nearest UMAP neighbours
  query_nn <- RANN::nn2(data = ref_obj$umap$embedding,
                        query = query_obj@reductions$umap@cell.embeddings, k = k, eps = 0)
  query_obj@meta.data[[initial_label]] <- apply(query_nn$nn.idx, 1, function(x) {median(ref_pseudotime[x,])})

  if(is.null(mapQC_class)){
    ## If no mapping error QC, final celltype same as initial.
    query_obj@meta.data[[final_label]] <- query_obj@meta.data[[initial_label]]
  } else if(is.na(mapQC_class)) {
    ## If no mapping error QC, final celltype same as initial.
    query_obj@meta.data[[final_label]] <- query_obj@meta.data[[initial_label]]
  } else {
    ## If mapping error QC, assign annotations of mapQC failing cells as NA.
    query_obj@meta.data[[final_label]] <- ifelse(query_obj@meta.data[[mapQC_class]] %in% c('Fail', 'fail'), NA,
                                                 query_obj@meta.data[[initial_label]])
  }

  return(query_obj)
}
