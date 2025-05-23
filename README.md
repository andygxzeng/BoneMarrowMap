# BoneMarrowMap

### Rapid and accurate projection of scRNA-seq data onto a curated atlas of Bone Marrow hematopoiesis with balanced representation of HSPCs and differentiated cells.

We provide a large-scale reference atlas of human hematopoiesis with balanced representation of CD34+ stem and progenitor cells together with terminally differentiated populations.
Our reference map of 263,159 hematopoietic cells provides us with deep portraits of lineage transitions in both early and late stages of human hematopoiesis. 

HSPC annotations within this atlas were carefully validated to maximize concordance with immunophenotypically pure populations
established through decades of functional studies. In particular, our transcriptional HSC state is highly concordant with functional LT-HSCs.

We provide functions to rapidly and accurately map query scRNA-seq profiles of either normal or leukemic hematopoietic cells onto our bone marrow atlas. 
Mapping, QC Filtering, CellType prediction, Pseudotime prediction, and Composition analysis can be performed from raw count matrices within minutes
(<10 min total for ~100,000 single cells on a personal laptop).

To cite this package or learn more, please refer to the publication in Blood Cancer Discovery: 
https://aacrjournals.org/bloodcancerdiscov/article/doi/10.1158/2643-3230.BCD-24-0342/761938/Single-cell-Transcriptional-Atlas-of-Human

To explore this dataset interactively please use our [CellxGene portal](https://cellxgene.cziscience.com/collections/f6c50495-3361-40ed-a819-fb9644396ed9), or download annotated single-cell objects can also be downloaded in scanpy format [(h5ad file)](https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_Annotated_Dataset_expandedFeatures.h5ad) or in Seurat format [(RDS file)](https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_Annotated_Dataset_expandedFeatures.rds).

![BoneMarrowMap](https://raw.githubusercontent.com/andygxzeng/BoneMarrowMap_Extras/main/BoneMarrow_ReferenceMap.png)


## Installation

```
# install dependencies
if(!require(BiocManager, quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AUCell", "doMC", "BiocNeighbors"))
if(!require(devtools, quietly = TRUE)) install.packages("devtools")
devtools::install_github("jaredhuling/jcolors") # dependency that is no longer on CRAN

# install BoneMarrowMap package
devtools::install_github('andygxzeng/BoneMarrowMap')
```

## [In-Depth Tutorial](https://htmlpreview.github.io/?https://github.com/andygxzeng/BoneMarrowMap_Extras/blob/main/BoneMarrowMap_Tutorial.nb.html)

We provide an [in-depth tutorial notebook](https://htmlpreview.github.io/?https://github.com/andygxzeng/BoneMarrowMap_Extras/blob/main/BoneMarrowMap_Tutorial.nb.html) for downloading the reference object and an example query dataset (Roy et al, Cell Rep 2021) and performing reference mapping, celltype classification, pseudotime prediction, and composition analysis on the example query dataset.

You can substitute the example dataset for your own scRNA-seq data. All you need is a raw or ambient RNA-corrected count matrix and cell annotations, ideally within a seurat object. Your data should be pre-filtered by n_counts, n_features, and percent_mito before mapping.

Query mapping is effective across sequencing technologies and hematopoietic tissues. Limitations include lack of coverage of mature Neutrophils (due to limitations in 10x technology) and T-cell precursors (due to tissue source).

To run the tutorial within R markdown, you can download the Rmd file directly from your R console and get started on the tutorial:
```
# Download R Notebook with mapping tutorial 
download.file('https://raw.githubusercontent.com/andygxzeng/BoneMarrowMap_Extras/main/BoneMarrowMap_Tutorial.Rmd', 
              destfile = 'BoneMarrowMap_Tutorial.Rmd')
```

## [Leukemia Projection Tutorial](https://htmlpreview.github.io/?https://github.com/andygxzeng/BoneMarrowMap_Extras/blob/main/LeukemiaProjection_Tutorial.nb.html)

We have also used this tool extensively for the mapping of leukemic blasts spanning all various acute and chronic leukemias, except T-ALL. Classification of leukemia depends on lineage-specific features identified by morphology and immunophenotyping. Projection of individual leukemic blasts based on thousands of transcriptomic lineage markers enables precise mapping of cellular states beyond what can be achieved in the clinic today. 

After completing the standard tutorial, please see the [leukemia projection tutorial](https://htmlpreview.github.io/?https://github.com/andygxzeng/BoneMarrowMap_Extras/blob/main/LeukemiaProjection_Tutorial.nb.html) to download an example dataset of three diverse AMLs and perform mapping, classification, pseudotime prediction, and composition analysis. This tutorial covers special considerations in mapping QC for leukemia samples and provides additional functions to score leukemia cells for functionally validated LSC signatures. 

To run the leukemia mapping tutorial within R markdown, you can download the Rmd file directly from your R console:
```
# Download R Notebook with mapping tutorial for leukemic samples
download.file('https://raw.githubusercontent.com/andygxzeng/BoneMarrowMap_Extras/main/LeukemiaProjection_Tutorial.Rmd', 
              destfile = 'LeukemiaProjection_Tutorial.Rmd')
```


## [B-cell Development Projection Tutorial](https://htmlpreview.github.io/?https://github.com/andygxzeng/b_development_map/blob/main/BALL_Projection_Example.nb.html)

For those interested in B-cell development and B-cell acute lymphoblastic leukemia, we have a second reference map of human B-cell development derived from fetal and post-natal tissue sources. This can be achieved through a second round of projection onto our B-cell development atlas. Please see our [B-cell Development projection tutorial](https://htmlpreview.github.io/?https://github.com/andygxzeng/b_development_map/blob/main/BALL_Projection_Example.nb.html) to download an example dataset of two phenotypically diverse B-ALL samples and perform mapping and classification using both reference atlases. To cite or learn more about this B-cell differentiation atlas and how we applied it to B-ALL, please see this pre-print: https://www.biorxiv.org/content/10.1101/2023.12.04.569954

To run the B-cell development mapping tutorial within R markdown, you can download the Rmd file directly from your R console:
```
# Download R Notebook with mapping tutorial for leukemic samples
download.file('https://raw.githubusercontent.com/andygxzeng/b_development_map/main/BALL_Projection_Example.Rmd', 
              destfile = 'BALL_Projection_Example.Rmd')
```


## Quick Start Guide

Through a few chunks of code, we will: 
  - Download and visualize the annotated Bone Marrow Reference Atlas
  - Download and load a query dataset of CD34+ HSPCs from Roy et al 2021
  - Map the query dataset and filter out low-quality cells 
  - Predict Hematopoietic Cell Type and Pseudotime score for each query cell
  - Visualize projection results and save the mapping results 
  
```
library(Seurat)
library(tidyverse)
library(symphony)
library(ggpubr)
library(patchwork)
library(BoneMarrowMap)
```

### Download reference object and UMAP model 

```
# Set directory to store projection reference files
projection_path = './'

# Download Bone Marrow Reference - 344 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_SymphonyReference.rds', 
                    destfile = paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
# Download uwot model file - 221 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_uwot_model.uwot', 
                    destfile = paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot'))

# Load Symphony reference
ref <- readRDS(paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
# Set uwot path for UMAP projection
ref$save_uwot_path <- paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot')

# Visualize Bone Marrow reference
ReferenceSeuratObj <- create_ReferenceObject(ref)
DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CellType_Annotation_formatted', raster=FALSE, label=TRUE, label.size = 4)
```

### Download example query data and load 
You can use your own query dataset, suggested format is a Seurat object with raw counts
```
# Load example data from Roy et al - 141 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/ExampleQuery_Roy2021.rds',
                    destfile = paste0(projection_path, 'ExampleQuery_Roy2021.rds'))

query <- readRDS(paste0(projection_path, 'ExampleQuery_Roy2021.rds'))
query
```

### Map Query data and evaluate mapping QC metrics (< 1min)
```
# batch variable to correct in the query data, set as NULL if no batches in query
batchvar <- 'sampleID'

# Map query dataset using Symphony 
query <- map_Query(
    exp_query = query[['RNA']]@counts, 
    metadata_query = query@meta.data,
    ref_obj = ref,
    vars = batchvar
)

# Run QC based on mapping error score, flag cells with mapping error >= 2.5 MADs above median
query <- query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5) 
plot_MappingErrorQC(query)
```

### Predict Cell Types from Hematopoietic Reference Map (10,000 cells in < 1min)
```
# Predict Hematopoietic Cell Types by KNN classification
query <- predict_CellTypes(
  query_obj = query, 
  ref_obj = ref, 
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

DimPlot(subset(query, mapping_error_QC == 'Pass'), group.by = c('predicted_CellType'), label=TRUE, label.size = 4)
```

### Predict Pseudotime along the Hematopoietic Hierarchy (< 1min)
```
# Predict Pseudotime values by KNN
query <- predict_Pseudotime(
  query_obj = query, 
  ref_obj = ref, 
  final_label = 'predicted_Pseudotime'
)

FeaturePlot(subset(query, mapping_error_QC == 'Pass'), features = c('predicted_Pseudotime'))
```

### Plot projection results, split by donor 
````
# Set batch/condition to be visualized individually
batch_key <- 'sampleID'

# returns a list of plots for each donor from a pre-specified batch variable
projection_plots <- plot_Projection_byDonor(
  query_obj = query, 
  batch_key = batch_key, 
  ref_obj = ref, 
  save_folder = 'projectionFigures/'
)

# show plots together with patchwork
patchwork::wrap_plots(projection_plots)
````

### Save annotations for each cell
```
save_ProjectionResults(
  query_obj = query, 
  file_name = 'querydata_projected_labeled.csv')
```

