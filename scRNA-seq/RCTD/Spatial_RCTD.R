library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spacexr)

# E9.5
# cell_position <- as.matrix(cbind(sobj_9.5@images$D02165C2@coordinates$imagerow,sobj_9.5@images$D02165C2@coordinates$imagecol))
# colnames(cell_position) <- c("Spatial_1","Spatial_2")
# rownames(cell_position) <- colnames(sobj_9.5)
# sobj_9.5[["Spatial_location"]] <- CreateDimReducObject(embeddings = cell_position, key = "Spatial_", assay = "RNA")
# p <- SpatialDimPlot(sobj_9.5, stroke = 0, pt.size = 50)
# ggsave("Spatial_DimPlot_E9_5.pdf",p,device = "pdf",width=12,height=8)
# 
# 
# cell_position <- as.matrix(cbind(sobj_10.5@images$D02165C3@coordinates$imagerow,sobj_10.5@images$D02165C3@coordinates$imagecol))
# colnames(cell_position) <- c("Spatial_1","Spatial_2") 
# rownames(cell_position) <- colnames(sobj_10.5)
# sobj_10.5[["Spatial_location"]] <- CreateDimReducObject(embeddings = cell_position, key = "Spatial_", assay = "RNA")
# p <- SpatialDimPlot(sobj_10.5, stroke = 0, pt.size = 60)
# ggsave("Spatial_DimPlot_E10_5.pdf",p,device = "pdf",width=12,height=8)

# sobj_9.5 <- SCTransform(sobj_9.5, assay = "Spatial", ncells = 3000, verbose = FALSE)
# sobj_9.5 <- RunPCA(sobj_9.5)
# sobj_9.5 <- RunUMAP(sobj_9.5, dims = 1:30)
# sobj_9.5 <- FindNeighbors(sobj_9.5, dims = 1:30)
# sobj_9.5 <- FindClusters(sobj_9.5, resolution = 0.3, verbose = FALSE)

# Read and update the Seurat object for the reference
ref <- readRDS("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/integration/sobj_x9_allcelltype_mnn_v4.rds")
ref <- UpdateSeuratObject(ref)
Idents(ref) <- "celltype"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$celltype)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)


# Create the reference object
reference <- Reference(counts, cluster, nUMI)

# Define the list of sample names
sample_names <- c("9.5","10.5","11.5","12.5","13.5","14.5","15.5","16.5","18.5")

# Loop over each sample
for (sample in sample_names) {
  # Process the current sample
  sobj_sample <- readRDS(paste0("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/spatial/rawdata/E", sample, ".rds"))
  query_counts <- sobj_sample[["Spatial"]]$counts
  coords <- GetTissueCoordinates(sobj_sample)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  
  # Create the query object
  query <- SpatialRNA(coords, query_counts, colSums(query_counts))
  
  # Run RCTD analysis
  RCTD <- create.RCTD(query, reference, max_cores = 4)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  # Add metadata and save the results
  sobj_sample <- AddMetaData(sobj_sample, metadata = RCTD@results$results_df)
  
  # Save rds files
  saveRDS(RCTD, file = paste0("Spatial_", sample, "_RCTD_Score.Rds"))
  saveRDS(sobj_sample, file = paste0("Spatial_", sample, "_RCTD.Rds"))
}

