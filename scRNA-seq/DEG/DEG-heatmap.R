library(Seurat)
library(tidyverse)
library(viridis)
library(RColorBrewer)

sobj <- readRDS("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/GC_subtype/sobj_x9_TB_mnn_v5.rds")
md <- sobj@meta.data
selected_cluster <- c("GC-Prl7b1", "GC-Aldh1a3")
setdiff(selected_cluster, sobj$celltype_v5)
sobj <- subset(sobj, subset = celltype_v5 %in% selected_cluster)
sobj$celltype <- sobj$celltype_v5
Idents(sobj) <- "celltype"
DimPlot(sobj)
md <- sobj@meta.data
str(md)

new_order <- c("GC-Aldh1a3", "GC-Prl7b1")
# new_order <- rev(new_order)

annotation <- unique(sobj$celltype)
setdiff(new_order, annotation)
setdiff(annotation, new_order)

# genes <- rev(genes)
all_genes <- row.names(sobj)
setdiff(genes, all_genes)

sobj$celltype <- factor(sobj$celltype, levels = new_order)
Idents(sobj) <- "celltype"
DimPlot(sobj, label = TRUE)

# GC-DEG
both <- c("Gbe1", "Tbc1d1", "Pygl", "Lrrc8d", "Mcu", "Tiam1", "Stat3", "Ppp3ca", "Ogt", "Zbtb20", "Hk1", "Raf1", "Map4k4")
GC.P <- c("Ccbe1", "Pik3cb", "Tgfbr3", "Btg1", "Hpgd", "Prdm1", "Zmiz1", "Col3a1", "Prl2a1", "Prl7b1", "Prl7c1", "Prl8a6", "Slc12a2", "Slc4a8", "Bcl2l1", "H2-D1", "Bst2", "H2-K1")
GC.A <- c("Tbc1d8", "Ralgps2", "Tiam2", "Actn1",  "Vstm4", "Tns1")
genes <- c(both, GC.A, GC.P)

colors <- rev(brewer.pal(n = 11, name = "RdBu"))
sobj <- ScaleData(sobj, features = genes)
p <- DoHeatmap(object = sobj, features = genes, draw.lines = FALSE, raster = FALSE, group.colors = c("#CCCC00", "#FF66FF")) + 
  scale_fill_gradientn(limits = c(-2, 2), colors = colors)
p
ggsave("DEG_heatmap.pdf", plot = p, width = 8, height = 10)