#' Per-cell Confidence Score:
#'
#' Calculates the weighted Mahalanobis distance for the query cells to reference clusters. Returns a vector
#' of distance scores, one per query cell. Higher distance metric indicates less confidence.
#'
#' @param reference Custom Symphony Reference object with covariance matrix for centroids
#' @param query Query Seurat object as returned by map_Query()
#' @param MAD.threshold Median absolute deviation cutoff to identify query cells with low-quality mapping
#'
#' @importFrom stats mahalanobis
#' @importFrom stats median
#' @importFrom stats mad
#' @return A vector of per-cell mapping metric scores for each cell.
#' @export
#'
calculate_MappingError = function(query, reference, MAD.threshold = 2) {

  query_pca = t(query@reductions$pca@cell.embeddings)
  query_R = query@reductions$harmony@misc$R

  # Calculate the Mahalanobis distance from each query cell to all centroids
  mah_dist_ks = matrix(rep(0, len = ncol(query_pca) * ncol(reference$centroids)), nrow = ncol(query_pca))
  for (k in 1:ncol(reference$centroids)) {
    mah_dist_ks[, k] = sqrt(stats::mahalanobis(x = t(query_pca), center = reference$center_ks[, k], cov = reference$cov_ks[[k]]))
  }

  # Return the per-cell score, which is the average of the distances weighted by the clusters the cell belongs to
  maha = rowSums(mah_dist_ks * t(query_R))
  query$mapping_error_score <- maha
  query$mapping_error_QC <- ifelse(query$mapping_error_score < (stats::median(query$mapping_error_score) + MAD.threshold*stats::mad(query$mapping_error_score)), 'Pass', 'Fail')
  return(query)
}
