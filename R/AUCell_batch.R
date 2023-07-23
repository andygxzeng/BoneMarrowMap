#' Run AUCell Scoring in Batches:
#'
#' Split Cell x Gene expression matrix into batches and run AUCell to reduce working memory load.
#'
#' @param inp_data Cell x Gene expression matrix with cells as columns and genes as rows. Raw counts are acceptable.
#' @param genesets List of lists, which each component list being the genes that belong to an individual geneset to score by AUCell.
#' @param num_batches Number of batches to split the Cell x Gene expression matrix into. These will be processed sequentially to save memory.
#' @param num_cores Number of cores to use for parallel computation of AUCell scores.
#'
#' @importFrom AUCell AUCell_buildRankings
#' @importFrom AUCell AUCell_calcAUC
#' @importFrom SummarizedExperiment assay
#' @import doMC
#' @return A matrix of AUCell scores for each geneset for each cell in the input dataset
#'
AUCell_batch <- function(inp_data, genesets, num_batches = 10, num_cores = 1) {
  num_cells <- ncol(inp_data)
  batch_size <- ceiling(num_cells/num_batches)
  score_mat <- c()
  print('Running AUCell scoring')
  for (i in 1:num_batches) {
    print(paste('batch', i))
    ind1 <- (i-1)*batch_size + 1
    ind2 <- i*batch_size
    if (ind2 > num_cells) {
      ind2 <- num_cells
    }
    gene_rankings <- AUCell::AUCell_buildRankings(inp_data[,ind1:ind2], plotStats = FALSE, nCores = num_cores) #splitByBlocks=TRUE
    score_mat_i <- AUCell::AUCell_calcAUC(geneSets = genesets, rankings = gene_rankings, nCores = num_cores)
    score_mat_i <- t(SummarizedExperiment::assay(score_mat_i, 'AUC'))
    score_mat <- rbind(score_mat, score_mat_i)
    gc(full = TRUE, verbose = TRUE)
  }
  print('Finished Scoring')
  return(score_mat)
}
