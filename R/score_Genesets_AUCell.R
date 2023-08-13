#' Scoring Genesets on Seurat object by AUCell:
#'
#' Split Cell x Gene expression matrix into batches and run AUCell to reduce working memory load.
#'
#' @param scrna Seurat object to be scored by AUCell. "counts" slot within RNA assay will be used so this must not be empty.
#' @param genesets List of lists, which each component list being the genes that belong to an individual geneset to score by AUCell.
#' @param nbatches Number of batches to split the Cell x Gene expression matrix into. These will be processed sequentially to save memory.
#' @param ncores Number of cores to use for parallel computation of AUCell scores.
#' @param output Preferred output format of AUCell scores. This can either be "metadata", "assay", or "dataframe". "metadata" returns a seurat object with AUCell scores in the metadata. "assay" returns a seurat object with AUCell scores as a separate assay. "dataframe" returns a dataframe with AUCell scores for each cell.
#' @param assay_name If preferred output is "assay", specify the name of the assay to store the results in. Default is "AUCell".
#'
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat AddMetaData
#' @importFrom Seurat CreateAssayObject
#'
#' @return Depending on output parameter, either a seurat object with AUCell scores within metadata or a separate assay, or AUCell scores as a separate dataframe
#' @export
#'
score_Genesets_AUCell <- function(scrna, genesets, nbatches = 10, ncores = 1, output = c('metadata', 'assay', 'dataframe'), assay_name = 'AUCell'){

  AUCell_scores <- AUCell_batch(Seurat::GetAssayData(scrna, assay = 'RNA', slot = 'counts'), genesets=genesets, num_batches=nbatches, num_cores=ncores)
  colnames(AUCell_scores) <- paste0(colnames(AUCell_scores), "_AUC")

  if(output == 'metadata'){
    scrna <- Seurat::AddMetaData(scrna, as.data.frame(AUCell_scores))
    out <- scrna
  }else if(output == 'assay'){
    scrna[[assay_name]] <- Seurat::CreateAssayObject(t(AUCell_scores))
    out <- scrna
  }else if(output == 'dataframe'){
    AUCell_scores <- as.data.frame(AUCell_scores)
    out <- AUCell_scores
  }else{
    stop('output parameter \"{output}\" is not a valid option. Please specify output as either "metadata", "assay", or "dataframe".')
  }
  return(out)
}
