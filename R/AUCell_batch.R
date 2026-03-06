#' Run AUCell Scoring in Batches:
#'
#' Split Cell x Gene expression matrix into batches and run AUCell to reduce working memory load.
#'
#' @param inp_data Cell x Gene expression matrix with cells as columns and genes as rows. Raw counts are acceptable.
#' @param genesets List of lists, which each component list being the genes that belong to an individual geneset to score by AUCell.
#' @param num_batches Number of batches to split the Cell x Gene expression matrix into. These will be processed sequentially to save memory.
#' @param ncores Number of cores to use for parallel computation. Default is 1 (serial).
#'   For newer AUCell versions, BoneMarrowMap uses AUCell_run (and BiocParallel if available).
#'   For older AUCell versions, it automatically falls back to AUCell_buildRankings + AUCell_calcAUC with nCores.
#'
#' @importFrom AUCell AUCell_buildRankings
#' @importFrom AUCell AUCell_calcAUC
#' @importFrom SummarizedExperiment assay
#' @return A matrix of AUCell scores for each geneset for each cell in the input dataset
#'
AUCell_batch <- function(inp_data, genesets, num_batches = 10, ncores = 1) {
  has_aucell_run <- "AUCell_run" %in% getNamespaceExports("AUCell")
  has_biocparallel <- requireNamespace("BiocParallel", quietly = TRUE)

  # Build BPPARAM only when AUCell_run path is available
  BPPARAM <- NULL
  if (has_aucell_run && has_biocparallel) {
    if (ncores > 1) {
      if (Sys.info()[["sysname"]] == "Windows") {
        BPPARAM <- BiocParallel::SnowParam(workers = ncores)
      } else {
        BPPARAM <- BiocParallel::MulticoreParam(workers = ncores)
      }
    } else {
      BPPARAM <- BiocParallel::SerialParam()
    }
  }

  if (has_aucell_run) {
    message("Using AUCell::AUCell_run backend")
  } else {
    message("Using legacy AUCell backend (AUCell_buildRankings + AUCell_calcAUC)")
  }

  num_cells <- ncol(inp_data)
  batch_size <- ceiling(num_cells/num_batches)
  score_mat <- c()
  message("Running AUCell scoring")
  for (i in 1:num_batches) {
    message(paste("batch", i))
    ind1 <- (i-1)*batch_size + 1
    ind2 <- i*batch_size
    if (ind2 > num_cells) {
      ind2 <- num_cells
    }

    if (has_aucell_run) {
      if (is.null(BPPARAM)) {
        score_mat_i <- AUCell::AUCell_run(inp_data[, ind1:ind2], geneSets = genesets)
      } else {
        score_mat_i <- AUCell::AUCell_run(inp_data[, ind1:ind2], geneSets = genesets, BPPARAM = BPPARAM)
      }
      score_mat_i <- t(SummarizedExperiment::assay(score_mat_i, "AUC"))
    } else {
      ncores_legacy <- ncores
      if (Sys.info()[["sysname"]] == "Windows") {
        ncores_legacy <- 1
      }
      gene_rankings <- AUCell::AUCell_buildRankings(inp_data[, ind1:ind2], plotStats = FALSE, nCores = ncores_legacy)
      score_mat_i <- AUCell::AUCell_calcAUC(geneSets = genesets, rankings = gene_rankings, nCores = ncores_legacy)
      score_mat_i <- t(SummarizedExperiment::assay(score_mat_i, "AUC"))
    }

    score_mat <- rbind(score_mat, score_mat_i)
    gc(full = TRUE, verbose = TRUE)
  }
  message("Finished Scoring")
  return(score_mat)
}
