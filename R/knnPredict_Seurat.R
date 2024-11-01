#' KNN Prediction of annotation label from reference to query
#'
#' This is modified from the original Seurat_Utils.R script within the Symphony github page
#'
#' @param query_obj Query Seurat object as returned by map_Query
#' @param ref_obj Custom Symphony Reference object
#' @param label_transfer Annotation column from Symphony Reference object to use for classification
#' @param col_name Name of column in query object to store annotations
#' @param k Number of nearest neighbours from Symphony Reference to use for classification
#' @param confidence Whether to return KNN confidence scores (proportion of neighbours voting for the predicted annotation)
#' @param seed random seed (k-NN has some stochasticity in the case of ties)
#'
#' @importFrom class knn
#' @importFrom Seurat Embeddings
#' @return A Seurat object, with predicted reference labels stored in the 'col_name' column of the meta.data
#' @export
knnPredict_Seurat <- function(query_obj, ref_obj, label_transfer, col_name, k = 5, confidence = TRUE, seed = 0)
{
  set.seed(seed)
  if (!label_transfer %in% colnames(ref_obj$meta_data)) {
    stop('Label \"{label_transfer}\" is not available in the reference metadata.')
  }

  if (confidence) {
    knn_pred <- class::knn(t(ref_obj$Z_corr), Seurat::Embeddings(query_obj, 'harmony_projected'),
                           ref_obj$meta_data[[label_transfer]], k = k, prob = TRUE)
    knn_prob = attributes(knn_pred)$prob
    query_obj@meta.data[[col_name]] <- knn_pred
    query_obj@meta.data[paste0(col_name, '_prob')] = knn_prob
  } else {
    knn_pred <- class::knn(t(ref_obj$Z_corr), Seurat::Embeddings(query_obj, 'harmony_projected'),
                           ref_obj$meta_data[[col_name]], k = k, prob = FALSE)
    query_obj@meta.data[[col_name]] <- knn_pred
  }
  return(query_obj)
}
