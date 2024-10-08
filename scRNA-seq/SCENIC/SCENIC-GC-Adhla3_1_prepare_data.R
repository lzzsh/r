library(Seurat)
library(Matrix)
library(tidyverse)
library(openxlsx)

sobj <- readRDS("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/GC_subtype/sobj_x9_TB_mnn_v5.rds")
sobj$celltype <- sobj$celltype_v5
md <- sobj@meta.data
annotation <- md$celltype_v5 %>% unique()

filtered_cluster <- c(
  "JZP",
  "LaTP",
  "GC precursor",
  "SynTI",
  "S-TGC precursor",
  "SpT precursor",
  "S-TGC",
  "P-TGC",
  "SynTII",
  "GC-Prl7b1",
  "LaTP2",
  "SynTII precursor",
  "SynTI precursor",
  "SpT-outer",
  "SpT-inner",
  "GC-Aldh1a3"
)
# class_df <- read.csv("../../2023-12-27_Subtypes_with_RCTD/data/cell_types.csv", row.names = 1)
# filtered_cluster <- class_df %>% filter(class == "TB") %>% row.names()

sobj <- subset(sobj, subset = celltype %in% filtered_cluster)

# Extract barcodes from the original Seurat object
all_barcodes <- colnames(sobj)
set.seed(123) # Setting a seed for reproducibility
selected_barcodes <- sample(all_barcodes, 10000)
sobj <- subset(sobj, cells = selected_barcodes)

Idents(sobj) <- "celltype"
DimPlot(sobj, label = TRUE)

# sobj$celltype <- sobj$integrated_snn_res.0.5
# sobj[["RNA"]] <- CreateAssayObject(counts = sobj@assays$Spatial@counts)
# sobj[["Spatial"]] <- NULL
# sobj[["SCT"]] <- NULL
# DefaultAssay(sobj) <- "RNA"
# sobj[["integrated"]] <- NULL
# saveRDS(sobj, file = "sobj_E18.5.rds")

mat <- sobj@assays[["RNA"]]@layers[["counts"]]
row.names(mat) <- row.names(sobj)
colnames(mat) <- colnames(sobj)
md <- sobj@meta.data
cellInfo <- md[ , c("orig.ident", "celltype")]
colnames(cellInfo)[2] <- "CellType"

markers <- read.xlsx("DEG_TB.xlsx", sheet = "GC-Aldh1a3")
# filtered_markers <- markers %>% filter(cluster == "GC-Aldh1a3")
# genes <- filtered_markers$gene
genes <- markers$gene
mat <- mat[genes, ]

saveRDS(mat, file = "exprMat_TB.Rds")
saveRDS(cellInfo, file = "cellInfo_TB.Rds")
