#' Create Seurat Object from Symphony Reference. No Gene Expression data but can visualize annotations
#'
#' @param ref_obj Custom Symphony Reference object
#'
#' @importFrom dplyr rename
#' @importFrom tibble column_to_rownames
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat CreateDimReducObject
#' @return A Seurat object with embeddings and metadata from the reference object. To save space, no raw count data is present in this object.
#' @export
#'
create_ReferenceObject = function(ref_obj){
  ReferenceSeuratObj <- Seurat::CreateSeuratObject(data.frame(Cell = ref_obj$meta_data$Cell, placeholder = 0) %>% tibble::column_to_rownames('Cell') %>% data.matrix() %>% t(),
                                                   meta.data = ref_obj$meta_data %>% tibble::column_to_rownames('Cell'))

  refUMAP <- data.frame(ref_obj$umap$embedding) %>% dplyr::rename(umap_1 = X1, umap_2 = X2) %>% data.matrix()
  rownames(refUMAP) <- ref_obj$meta_data$Cell
  ReferenceSeuratObj@reductions[['umap']] <- Seurat::CreateDimReducObject(refUMAP, key='umap')

  return(ReferenceSeuratObj)
}
