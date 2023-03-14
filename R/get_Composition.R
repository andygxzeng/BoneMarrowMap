#' Get Composition for projected CellTypes
#'
#' @param query_obj Query Seurat object as returned by Symphony map_Query
#' @param donor_key Metadata column specifying donor to tally celltype abundance by
#' @param celltype_label Column in query metadata with projected CellType classification.
#' @param mapQC_col Column in query object with mapping error QC status for each cell (either "Pass" or "Fail"). Can set as NULL if no filtering desired.
#' @param knn_prob_cutoff CellType classification KNN probability threshold for composition analysis. Default is NULL.
#' @param return_type Format to return composition data, either matrix with counts, matrix with proportions, or in a long format.
#'
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @return A data frame with celltype abundance data for each donor.
#' @export
get_Composition = function(query_obj, donor_key, celltype_label = 'predicted_CellType', mapQC_col = 'mapping_error_QC',
                           knn_prob_cutoff = NULL, return_type = c('long', 'count', 'proportion')){

  if (!donor_key %in% colnames(query_obj@meta.data)) {
    stop('Label \"{donor_key}\" is not available in the query metadata.')
  }
  if (!celltype_label %in% colnames(query_obj@meta.data)) {
    stop('Label \"{celltype_label}\" is not available in the query metadata.')
  }
  if (!paste0(celltype_label,'_prob') %in% colnames(query_obj@meta.data)) {
    stop('Label \"{celltype_label}\" KNN prob is not available in the query metadata.')
  }
  if (is.null(knn_prob_cutoff)) {
    knn_prob_cutoff <- 0
  } else if (is.na(knn_prob_cutoff)) {
    knn_prob_cutoff <- 0
  }

  # if mapQC col is not NA or NULL
  if(!is.null(mapQC_col)){
    if(!is.na(mapQC_col)) {
      query_composition <- query_obj@meta.data %>%
        dplyr::filter(.data[[paste0(celltype_label,'_prob')]] > knn_prob_cutoff) %>%
        dplyr::group_by(.data[[donor_key]], .data[[celltype_label]]) %>%
        dplyr::summarise(count = n()) %>%
        dplyr::group_by(.data[[donor_key]]) %>%
        dplyr::mutate(proportion = count / sum(count)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(-proportion) %>% dplyr::arrange(.data[[donor_key]])
    }
  } else {
    # check mapQC column
    if (!mapQC_col %in% colnames(query_obj@meta.data)) {
      stop('Label \"{mapQC_col}\" is not available in the query metadata.')
    }
    # Filter mapQC == 'Pass' before getting composition
    query_composition <- query_obj@meta.data %>%
      dplyr::filter(.data[[mapQC_col]] %in% c('Pass', 'pass', 'PASS')) %>%
      dplyr::filter(.data[[paste0(celltype_label,'_prob')]] > knn_prob_cutoff) %>%
      dplyr::group_by(.data[[donor_key]], .data[[celltype_label]]) %>%
      dplyr::summarise(count = n()) %>%
      dplyr::group_by(.data[[donor_key]]) %>%
      dplyr::mutate(proportion = count / sum(count)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(-proportion) %>% dplyr::arrange(.data[[donor_key]])
  }

  # If return as count or proportion matrix
  if(return_type %in% c('count', 'proportion')){
    query_composition <- query_composition %>%
      tidyr::pivot_wider(id_cols = donor_key, names_from = celltype_label, values_from = return_type) %>% replace(is.na(.), 0)
  }

  return(query_composition)
}
