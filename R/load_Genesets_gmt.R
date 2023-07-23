#' Load Geneset:
#'
#' Loads genesets stored in a gmt file into R as a list of lists.
#'
#' @param path path to geneset file, should end with .gmt
#' @param ignore_cols Number of columns to skip from gmt file for extracting genes. First two columns are typically geneset name and description, so these will be skipped when extracting a list of genes for each geneset.
#'
#' @return A list of lists, with each component list being the genes belonging to an individual geneset.
#' @export
#'
load_Genesets_gmt <- function(path, ignore_cols = 2){
  x <- scan(path, what="", sep="\n")
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  for(i in 1:ignore_cols){
    y <- lapply(y, `[`, -1)
  }
  return(y)
}
