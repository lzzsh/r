# Install the CARD and MuSiC packages from GitHub
devtools::install_github('YingMa0107/CARD')  # Install CARD
devtools::install_github('xuranw/MuSiC')     # Install the dependency package

### Load the CARD, MuSiC, Seurat, patchwork, and tidyverse libraries
library(CARD)
library(MuSiC)
library(Seurat)
library(patchwork)
library(tidyverse)

# Load the GBM4 spatial transcriptomics data saved in a previous session
load('GBM4.rdata')

### Get the counts expression matrix from the spatial transcriptomics data
spatial_count <- GBM4@assays$Spatial@counts
print(spatial_count[1:4,1:4])  # Display a small portion of the count matrix

### Get the spatial location matrix
spatial_loca <- GBM4@images$GBM4@coordinates
spatial_location <- spatial_loca[,2:3]
# The column names must be 'x' and 'y', or else CARD_deconvolution will throw an error
colnames(spatial_location) <- c("x", "y")
print(spatial_location[1:4,])  # Display a small portion of the location matrix

# Load the GBM-scRNA single-cell RNA sequencing data saved in a previous session
load('scRNA.rdata')

# Get the counts matrix from the single-cell RNA-seq data
sc_count <- scRNA@assays$RNA@counts

# Get the cell annotation matrix from the single-cell RNA-seq data
sc_meta <- scRNA@meta.data %>%
  rownames_to_column("cellID") %>%
  dplyr::select(cellID, orig.ident, celltype) %>%
  mutate(CB = cellID) %>%
  column_to_rownames("CB")
print(head(sc_meta))  # Display the first few rows of the meta data

## Create a CARD object
CARD_obj <- createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "celltype",
  ct.select = unique(sc_meta$celltype),  # Column name for cell types
  sample.varname = "orig.ident"
)

## Perform quality control on both the scRNASeq and spatially-resolved datasets
## ...

# Deconvolute using CARD
CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)

# Visualize the cell type distribution as pie charts at each spot
colors <- c("#4DAF4A","#F0027F","#377EB8","#FDC086","#A6761D","#FFFF00","#BEAED4",
            "#BF5B17","#666666")
p1 <- CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD,
                         spatial_location = CARD_obj@spatial_location,
                         colors = colors)

# Plot the spatial regions (e.g., normal vs. tumor)
p2 <- SpatialPlot(GBM4, group.by = 'Region', cols = c('Normal'='#007799','Tumor'='#AA0000'))

# Combine the two plots
p1 + p2

# Visualize the proportion of selected cell types across spots
ct.visualize <- c("OPC like", "AC like", "MES like")
p3 <- CARD.visualize.prop(proportion = CARD_obj@Proportion_CARD,
                          spatial_location = CARD_obj@spatial_location,
                          ct.visualize = ct.visualize,
                          colors = c("lightblue", "lightyellow", "red"),
                          NumCols = 3, pointSize = 1)  # Size of the points in the plot
print(p3)

# Visualize the original gene expression levels
p4 <- CARD.visualize.gene(
  spatial_expression = CARD_obj@spatial_countMat,
  spatial_location = CARD_obj@spatial_location,
  gene.visualize = c("SNAP25", "SYT1", "JUNB"),
  colors = NULL,
  NumCols = 3
)
print(p4)

# Save the CARD deconvolution results
save(CARD_obj, file = 'CARD_obj.rdata')