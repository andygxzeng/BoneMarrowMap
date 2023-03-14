# BoneMarrowMap

### Rapid and accurate projection of scRNA-seq data onto a curated atlas of Bone Marrow hematopoiesis with balanced representation of HSPCs and differentiated cells.

Our goal in developing this resource is to provide a large-scale and comprehensive reference atlas of human hematopoiesis in the bone marrow, 
and critically to ensure balanced representation of CD34+ stem and progenitor cells together with CD34- terminally differentiated populations.
This provides us with deep portraits of lineage transitions in both early and late stages of human hematopoiesis. 

HSPC annotations within this atlas were carefully validated to maximize concordance with immunophenotypically pure populations
established through decades of functional studies. In particular, our transcriptional HSC state is highly concordant with functional LT-HSCs.

We provide functions to rapidly and accurately map query scRNA-seq profiles of either normal or leukemic hematopoietic cells onto our bone marrow atlas. 
Mapping, QC Filtering, CellType prediction, Pseudotime prediction, and Composition analysis can be performed from raw count matrices within minutes
(~10min total for ~10,000 cells on a personal laptop).





## Quick Start

Through a few chunks of code, we will: 
  - Download and visualize the annotated Bone Marrow Reference Atlas
  - Download and load a query dataset of CD34+ HSPCs from Roy et al 2021
  - Map the query dataset and filter out low-quality cells 
  - Predict Hematopoietic Cell Type and Pseudotime score for each query cell
  - Visualize projection results and save the mapping results 
  

### Download reference object and UMAP model 

```
# Set directory to store projection reference files
projection_path = './'

# Download Bone Marrow Reference - 344 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_SymphonyRef.rds', 
                    destfile = paste0(projection_path, 'BoneMarrow_RefMap_SymphonyRef.rds'))
# Download uwot model file - 221 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_uwot_model.uwot', 
                    destfile = paste0(projection_path, 'BoneMarrow_RefMap_uwot_model.uwot'))

# Load Symphony reference
ref <- readRDS(paste0(projection_path, 'BoneMarrow_RefMap_SymphonyRef.rds'))
# Set uwot path for UMAP projection
ref$save_uwot_path <- paste0(projection_path, 'BoneMarrow_RefMap_uwot_model.uwot')

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
    exp_query = query@assays$RNA@counts, 
    metadata_query = query@meta.data,
    ref_obj = ref,
    vars = batchvar
)

# Run QC based on mapping error score, flag cells with mapping error >= 2.5 MADs above median
query <- query %>% calculate_MappingError(., reference = ref, MAD.threshold = 2.5) 
plot_MappingErrorQC(query)
```

### Predict Cell Types from Hematopoietic Reference Map (10,000 cells in < 10min)
```
# Predict Hematopoietic Cell Types by KNN classification
query <- predict_CellTypes(
  query_obj = query, 
  ref_obj = ref, 
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap', group.by = c('predicted_CellType'), raster=FALSE, label=TRUE, label.size = 4)
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

