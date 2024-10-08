library(Seurat)
library(tidyverse)
library(openxlsx)
library(viridis)
job_id <- "855"

sobj <- readRDS("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/GC_subtype/sobj_x9_TB_mnn_v5.rds")
md <- sobj@meta.data
# sobj <- subset(sobj, subset = stage %in% "WT")
# md <- sobj@meta.data
sobj$celltype <- sobj$celltype_v5
unique(sobj$celltype)
# selected_cluster <- c("Macrophage-D", "Macrophage-F")
# setdiff(selected_cluster, sobj$celltype)
# sobj <- subset(sobj, subset = celltype %in% selected_cluster)
Idents(sobj) <- "celltype"
DimPlot(sobj, label = TRUE)

# md <- sobj@meta.data
# Idents(sobj) <- "celltype"
output_file <- "DEG_TB"

# DefaultAssay(sobj) <- "RNA"
# md <- sobj@meta.data

output_xlsx <- paste0( output_file, ".xlsx")

# sobj_file <- paste0("../data/sobj_BGI.rds")
# sobj <- readRDS(sobj_file)
# DefaultAssay(sobj) <- "Spatial"
# DimPlot(sobj, label = TRUE, raster = FALSE)
# Subset the Seurat object by each stage
# sobj <- subset(sobj, subset = stage == current_stage)
# print(paste("Dimensions of sobj for stage", current_stage, ":", paste(dim(sobj), collapse = " x ")))

clusters <- Idents(sobj) %>% unique() %>% sort()
cluster_sizes <- table(Idents(sobj))
cluster_sizes
# levels(sobj)
# clusters <- c("Decidual-center", "Decidual-top", "Decidual-surrounding")

markers_list <- list()
all_markers <- NULL

for(cluster_id in clusters){
  # Skip cluster if it contains less than 3 cells
  if(cluster_sizes[[cluster_id]] < 3){
    print(paste("Skipping cluster ", cluster_id, " because it has fewer than 3 cells"))
    next
  }
  
  print(paste("Processing cluster: ", cluster_id))
  # Find markers for a specific group
  markers <- FindMarkers(sobj, ident.1 = cluster_id, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers <- cbind(row.names(markers), markers)
  colnames(markers)[1] <- "gene"
  markers$cluster <- cluster_id
  markers <- markers[order(-markers$avg_log2FC), ]
  row.names(markers) <- NULL
  markers_list[[cluster_id]] <- markers
  all_markers <- rbind(all_markers, markers)
}

saveRDS(all_markers, file = all_markers_rds)
write.csv(all_markers, file = all_markers_csv, row.names = FALSE)
markers_list[["all_markers"]] <- all_markers

# Define a function to create the Excel file
create_zoomed_excel <- function(data_list, sheet_names, filename, zoom = 200) {
  wb <- createWorkbook()
  
  for(i in seq_along(data_list)) {
    addWorksheet(wb, sheet_names[i], zoom = zoom)
    markers <- data_list[[i]]
    writeData(wb, sheet_names[i], data_list[[i]])
  }
  
  saveWorkbook(wb, filename, overwrite = TRUE)
}

# Now you can create an Excel file with multiple sheets, all with 200% zoom in one line:
create_zoomed_excel(markers_list, names(markers_list), output_xlsx)
