library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

sobj <- readRDS("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/integration/E9_5_to_18_5.Rds")
sobj@meta.data <- sobj@meta.data %>%
  mutate(doublets = case_when(sobj$nFeature_RNA > 3200 ~ "doublet",
                              TRUE ~ "normal"))

p11 <- DimPlot(
  sobj,
  reduction = "umap.rpca",
  group.by = "doublets",
  combine = FALSE, label.size = 2,
  cols = c("firebrick", "lightgrey")
)

p22 <- DimPlot(
  sobj,
  reduction = "umap.har",
  group.by = "doublets",
  combine = FALSE, label.size = 2,
  cols = c("firebrick", "lightgrey")
)


p33 <- DimPlot(
  sobj,
  reduction = "umap.mnn",
  group.by = "doublets",
  combine = FALSE, label.size = 2,
  cols = c("firebrick", "lightgrey")
)
pp <- wrap_plots(c(p11, p22, p33), ncol = 3, byrow = F)

ggsave("clus_annotation_before_rmDoublets_nFeature.pdf",pp,device = "pdf",width=35,height=10)

# remove doublets
sobj<-subset(sobj,subset=nFeature_RNA < 3200)
sobj[["RNA"]] <- split(sobj[["RNA"]], f = sobj$sample_name)

sobj <- IntegrateLayers(
  object = sobj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

sobj <- IntegrateLayers(
  object = sobj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

sobj <- IntegrateLayers(
  object = sobj, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

sobj <- FindNeighbors(sobj, reduction = "integrated.rpca", dims = 1:30)
sobj <- FindClusters(sobj, resolution = 0.1, cluster.name = "rpca_clusters")
sobj <- RunUMAP(sobj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:30)
sobj <- FindClusters(sobj, resolution = 0.1, cluster.name = "har_clusters")
sobj <- RunUMAP(sobj, reduction = "harmony", dims = 1:30, reduction.name = "umap.har")


sobj <- FindNeighbors(sobj, reduction = "integrated.mnn", dims = 1:30)
sobj <- FindClusters(sobj, resolution = 0.1, cluster.name = "mnn_clusters")
sobj <- RunUMAP(sobj, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")

p11 <- DimPlot(
  sobj,
  reduction = "umap.rpca",
  group.by = "celltype",
  combine = FALSE, label.size = 2
)

p22 <- DimPlot(
  sobj,
  reduction = "umap.har",
  group.by = "celltype",
  combine = FALSE, label.size = 2
)


p33 <- DimPlot(
  sobj,
  reduction = "umap.mnn",
  group.by = "celltype",
  combine = FALSE, label.size = 2
)
pp <- wrap_plots(c(p11, p22, p33), ncol = 3, byrow = F)
pp
ggsave("clus_annotation_after_rmDoublets_nFeature.pdf",pp,device = "pdf",width=35,height=10)
