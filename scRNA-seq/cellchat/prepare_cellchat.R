library(Seurat)

setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/GC_subtype/")

sobj <- readRDS("sobj_x9_TB_mnn_v5.rds")
DimPlot(sobj, label = TRUE, raster=FALSE)

# sobj <- subset(sobj, subset = celltype %in% c("Secondary P-TGC",
#                                               "_GC",
#                                               "SpA-TGC",
#                                               "_SpT",
#                                               "Decidual-center",
#                                               "Decidual-surrounding",
#                                               "Endothelial",
#                                               "_SynTI",
#                                               "_SynTII",
#                                               "_S-TGC"))
Idents(sobj) <- "celltype"
DimPlot(sobj, label = TRUE, raster=FALSE)

data.input <- GetAssayData(sobj, assay = "RNA", slot = "counts")
meta <- sobj@meta.data
meta$labels <- meta$celltype_v5
meta$labels <- as.factor(meta$labels)
unique(meta$labels)

saveRDS(data.input, file = "./cellchat/cellchat_data.input.rds")
saveRDS(meta, file = "./cellchat/cellchat_meta.rds")
