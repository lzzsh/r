library(CytoTRACE2)
library(RColorBrewer)
library(tidyverse)

seurat_obj <- readRDS("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/GC_subtype/sobj_x9_TB_mnn_v5.rds")

# UMAP
umap_data <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
umap_data$celltype <- seurat_obj$celltype_v5

colors <- c(
  "#669933",    "#339966",    "#339999",    "#85C247",    "#A3D175",    "#BDC25B",
  "#AC80D0",    "#D0AC80",    "#D080CC",    "#D08380",    "#D3D3D3",    "#D080A4",
  "#A6DBDE",    "#E633D1",    "#EC7969",    "#E63377",    "#E64933",    "#47A3FF",
  "#06D6A0",    "#FFFF0A",    "#FF8667",    "#E7C19C",    "#EEB1AA",    "#F5934E",
  "#DB8B0A",    "#A30059",    "#4AC170",    "#7BBC54",    "#FFB3FF",    "#BDBD2E",
  "#FFA64D",    "#CDCD32",    "#CD8032",    "#BD752E",    "#DDA873",    "#F990EE",
  "#D2CE2E")

cluster_names <- c("S-TGC precursor",			 
                   "LaTP",					 
                   "LaTP2",					 
                   "SynTI precursor",		 
                   "SynTII precursor",		 
                   "Proliferating EC",		 
                   "S-TGC",					 
                   "SynTII",					 
                   "SynTI",					 
                   "Angiogenic EC",			 
                   "ExE endoderm",			 
                   "Erythrocyte",				 
                   "Parietal endoderm",		 
                   "SpT precursor",			 
                   "JZP",						 
                   "GC precursor",			 
                   "P-TGC",					 
                   "SpT",						 
                   "GC",						 
                   "Nourishing DSC",			 
                   "DSC precursor",			 
                   "Angiogenic DSC",			 
                   "Endometrium stromal cell", 
                   "Epithelial",				 
                   "Venous EC",				 
                   "Mesenchyme",				 
                   "Megakaryocyte",			 
                   "Lymphatic EC",			 
                   "B cell",					 
                   "Monocyte",				 
                   "Neutrophil",				 
                   "DC",						 
                   "T cell",					 
                   "Macrophage",				 
                   "NK",
                   "GC-Prl7b1",
                   "GC-Aldh1a3")
umap_data$celltype <- factor(seurat_obj$celltype_v5,cluster_names)
names(colors) <- cluster_names

p <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = celltype)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = colors) +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "black", color = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "white"),
    legend.background = element_rect(fill = "black", color = "black"),
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # Increase legend dot size
p

# ggsave("TB_UMAP.pdf",p,device = "pdf",width=12,height=8)


# extract JZP,GC precursor,GC-1,GC-2
cell_types_to_extract <- c("JZP", "GC precursor", "GC-Prl7b1", "GC-Aldh1a3")
Idents(seurat_obj) <- seurat_obj$celltype_v5
seurat_obj_selected <- subset(seurat_obj, idents = cell_types_to_extract)

# rerun UMAP 
seurat_obj_selected <- RunUMAP(seurat_obj_selected, reduction = "integrated.mnn", dims = 1:30, reduction.name = "selected.umap.mnn")
umap_data <- as.data.frame(Embeddings(seurat_obj_selected, reduction = "selected.umap.mnn"))
umap_data$celltype <- seurat_obj_selected$celltype_v5

# rotate the graph
umap_data$selectedumapmnn_1 <- -(umap_data$selectedumapmnn_1)
umap_data$selectedumapmnn_2 <- -(umap_data$selectedumapmnn_2)

# UMAP grouped by cell type
p <- ggplot(umap_data, aes(x = selectedumapmnn_1, y = selectedumapmnn_2, color = celltype)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = colors) +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "black", color = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "white"),
    legend.background = element_rect(fill = "black", color = "black"),
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # Increase legend dot size
p

# ggsave("GC_lineage_UMAP.pdf",p,device = "pdf",width=12,height=8)

# UMAP grouped by stage
umap_data$stage <- seurat_obj_selected$stage
colors <- colorRampPalette(brewer.pal(9, "Blues"))(9)
p <- ggplot(umap_data, aes(x = selectedumapmnn_1, y = selectedumapmnn_2, color = stage)) +
  geom_point(size = 0.3) +
  coord_fixed() +
  theme_minimal() +
  scale_color_manual(values = colors) + 
  theme(
    panel.background = element_rect(fill = "black", color = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "white"),
    legend.background = element_rect(fill = "black", color = "black"),
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # Increase legend dot size

ggsave("GC_lineage_UMAP_grouped_stage.pdf",p,device = "pdf",width=22,height=15)

# running CytoTRACE 2 main function - cytotrace2
# setting is_seurat = TRUE as input is a Seurat object
# setting slot_type = "counts" as the gene expression data is stored in the "counts" slot of the Seurat object
seurat_obj_selected[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap_data[,-3]), key = "selectedumapmnn_", assay = "RNA")
cytotrace2_result <- cytotrace2(seurat_obj_selected, 
                                is_seurat = TRUE, 
                                slot_type = "counts", 
                                batch_size = 10000, 
                                smooth_batch_size = 1000)
# saveRDS(cytotrace2_result, file = "../data/sobj_cytotrace2_result.rds")

md <- cytotrace2_result@meta.data

DimPlot(cytotrace2_result, group.by = "celltype_v5", label = TRUE)

metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(cytotrace2_result, features = "CytoTRACE2_Relative")
color <- colorRampPalette(brewer.pal(11, "RdYlGn"))(100) %>% rev()
FeaturePlot(cytotrace2_result, features = "CytoTRACE2_Relative") + scale_color_gradientn(colors = color)

p <- FeaturePlot(cytotrace2_result, features = "CytoTRACE2_Relative") + 
  scale_color_gradientn(
    colors = color,
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("0.0 (More diff.)", "1.0 (Less diff.)")
  ) + coord_fixed()

# ggsave(paste0("../plots/CytoTRACE2_TB_", job_id, ".pdf"), p, width = 7, height = 7, limitsize = FALSE)


# making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame(phenotype = seurat_obj_selected@meta.data$celltype_v5) %>% set_rownames(., colnames(seurat_obj_selected))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation, 
                  is_seurat = TRUE)

plots$CytoTRACE2_UMAP
plots$CytoTRACE2_Potency_UMAP
plots$CytoTRACE2_Relative_UMAP
plots$Phenotype_UMAP
plots$CytoTRACE2_Boxplot_byPheno

# UMAP splited by stage
seurat_obj_selected$celltype_v5 <- factor(seurat_obj_selected$celltype_v5, levels = c("JZP", "GC precursor", "GC-Aldh1a3", "GC-Prl7b1"))
p <- DimPlot(seurat_obj_selected, reduction = "selected.umap.mnn", group.by = "celltype_v5", split.by = "stage", cols = colors) + scale_y_reverse()
ggsave("GC_lineage_UMAP_splited_stage.pdf",p,device = "pdf",width=12,height=8)


# umap_data$stage <- seurat_obj_selected$stage
# p <- ggplot(umap_data, aes(x = selectedumapmnn_1, y = selectedumapmnn_2, color = celltype)) +
#   geom_point(size = 0.3) +
#   scale_color_manual(values = colors) +
#   coord_fixed() +
#   theme_minimal() +
#   theme(
#     panel.background = element_rect(fill = "black", color = "white"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks = element_blank(),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 14, color = "white"),
#     legend.background = element_rect(fill = "black", color = "black"),
#   ) +
#   guides(color = guide_legend(override.aes = list(size = 5)))+
#   facet_wrap(~ stage, ncol = 9)
# p


