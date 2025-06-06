% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_CellTypes.R
\name{predict_CellTypes}
\alias{predict_CellTypes}
\title{KNN Prediction of Hematopoietic Cell Types}
\usage{
predict_CellTypes(
  query_obj,
  ref_obj,
  ref_label = "CellType_Annotation",
  k = 30,
  mapQC_col = "mapping_error_QC",
  initial_label = "initial_CellType",
  final_label = "predicted_CellType",
  include_broad = TRUE
)
}
\arguments{
\item{query_obj}{Query Seurat object as returned by Symphony map_Query}

\item{ref_obj}{Custom Symphony Reference object}

\item{ref_label}{Annotation column from Symphony Reference object to use for classification}

\item{k}{Number of nearest neighbours from Symphony Reference to use for classification}

\item{mapQC_col}{Column in query object with mapping error QC status for each cell (either "Pass" or "Fail"). Can set as NULL if no filtering desired.}

\item{initial_label}{Column name to assign initial cell type assignments to, prior to mapQC filtering}

\item{final_label}{Column name to assign final cell type assignments to, following mapQC filtering}

\item{include_broad}{Whether to include broad celltype assignments. If so, ref_label must be assigned as "CellType_Annotation"}
}
\value{
A Seurat object, with predicted CellType labels stored in the 'final_label' column of the meta.data. All predicted labels including cells that failed QC can be found in 'initial_label'.
}
\description{
Uses knnPredict_Seurat to classify celltypes based on nearest neighbours from reference
Provides initial estimates for all cells regardless of mapping QC filter
Provides final estimates for mapQC passing cells and assigns mapQC failing cells as NA.
}
