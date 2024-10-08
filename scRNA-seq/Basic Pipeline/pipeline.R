###########cellranger aggr#############
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(batchelor)
library(SeuratWrappers)

options(future.globals.maxSize = 20 * 1024^3)  

# Load the data
sobj.data <- Read10X(data.dir = "/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/aggr_result/AGG_9_to_18/outs/count/filtered_feature_bc_matrix")
sobj <- CreateSeuratObject(counts = sobj.data,project = "sobj3k",min.cells = 3,min.features = 200)
# 138202 cells

# Calculate the percentage of mitochondrial
sobj[['percent.mt']] <- PercentageFeatureSet(sobj,pattern = "^mt-")
VlnPlot(sobj,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)

plot1 <- FeatureScatter(sobj,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(sobj,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1 + plot2

# sobj$sample <- unlist(lapply(rownames(sobj@meta.data),function(x) unlist(strsplit(x,"-"))[2]))
# sobj@meta.data <- sobj@meta.data %>%
#   mutate(sample_name = case_when(sobj$sample == "1" ~ "E10_5",
#                                  sobj$sample == "2" ~ "E11_5",
#                                  sobj$sample == "3" ~ "E12_5",
#                                  sobj$sample == "4" ~ "E14_5",
#                                  sobj$sample == "5" ~ "E15_5",
#                                  sobj$sample == "6" ~ "E16_5",
#                                  sobj$sample == "7" ~ "E9_5",
#                                  sobj$sample == "8" ~ "E13_5",
#                                  sobj$sample == "9" ~ "E18_5"))

# sobj$sample_name <- factor(sobj$sample_name,levels = c("E9_5",
#                                                        "E10_5",
#                                                        "E11_5",
#                                                        "E12_5",
#                                                        "E13_5",
#                                                        "E14_5",
#                                                        "E15_5",
#                                                        "E16_5",
#                                                        "E18_5"))

# Quality control
p1 <- sobj@meta.data %>% ggplot() +
  ggridges::geom_density_ridges(aes(x = nFeature_RNA, y = sample_name, fill = sample_name)) +
  theme_classic() +
  labs(title = "Density Plot of nFeature_RNA") + 
  theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_x_continuous(limits = c(0, 3000)) +
  geom_vline(xintercept = 500, size = 0.7) +
  scale_fill_manual(values = morandi_colors) +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(t = 20, b = 20, l = 0, r = 0)))
p1
ggsave("DensityPlot_nFeature.pdf",p1,device = "pdf",width=15,height=10)

p2 <- sobj@meta.data %>% ggplot() +
  ggridges::geom_density_ridges(aes(x = nCount_RNA, y = sample_name, fill = sample_name)) +
  theme_classic() +
  labs(title = "Density Plot of nCount_RNA") + 
  theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_x_continuous(limits = c(0, 5000)) +
  geom_vline(xintercept = 500, size = 0.7) +
  scale_fill_manual(values = morandi_colors) +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(t = 20, b = 20, l = 0, r = 0)))
p2
ggsave("DensityPlot_nCount.pdf",p2,device = "pdf",width=15,height=10)

# Set the cutoff of nFeature and percent.mt
sobj <- subset(sobj,subset=nFeature_RNA>500 & percent.mt<10)
# 121752 cells

sobj <- NormalizeData(sobj)

sobj <- FindVariableFeatures(sobj,selection.method = "vst",nfeatures = 2000)
top10 <- head(VariableFeatures(sobj),10)
plot1 <- VariableFeaturePlot(sobj)
plot2 <- LabelPoints(plot = plot1,points = top10,repel = TRUE)
plot2

all.genes <- rownames(sobj)
sobj <- ScaleData(sobj,features = all.genes)


# Reduction
sobj <- RunPCA(sobj,features = VariableFeatures(object=sobj))
print(sobj[["pca"]],dims = 1:5,nfeatures = 5)

# VizDimLoadings(sobj,dims = 1:2,reduction = "pca")
# DimHeatmap(sobj,dims = 1,cells = 500,balanced = TRUE)

ElbowPlot(sobj)


# Find neighbors and clusters
sobj <- FindNeighbors(sobj,dims = 1:30)
sobj <- FindClusters(sobj,resolution = seq(0.1, 0.2, 0.1))
head(Idents(sobj),5)


# Run UMAP
sobj <- RunUMAP(sobj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(sobj, reduction = "umap.unintegrated", group.by = c("sample_name", "seurat_clusters"))

# Integration
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
sobj <- FindClusters(sobj, resolution = 2, cluster.name = "rpca_clusters")
sobj <- RunUMAP(sobj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:30)
sobj <- FindClusters(sobj, resolution = 2, cluster.name = "har_clusters")
sobj <- RunUMAP(sobj, reduction = "harmony", dims = 1:30, reduction.name = "umap.har")

sobj <- FindNeighbors(sobj, reduction = "integrated.mnn", dims = 1:30)
sobj <- FindClusters(sobj, resolution = 2, cluster.name = "mnn_clusters")
sobj <- RunUMAP(sobj, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")

p1 <- DimPlot(
  sobj,
  reduction = "umap.rpca",
  group.by = c("sample_name", "seurat_clusters"),
  combine = FALSE, label.size = 2
)

p2 <- DimPlot(
  sobj,
  reduction = "umap.har",
  group.by = c("sample_name", "seurat_clusters"),
  combine = FALSE, label.size = 2
)

p3 <- DimPlot(
  sobj,
  reduction = "umap.mnn",
  group.by = c("sample_name", "seurat_clusters"),
  combine = FALSE, label.size = 2
)

p <- wrap_plots(c(p1, p2, p3), ncol = 3, byrow = F)
# ggsave("p.pdf",p,device = "pdf",width=20,height=10)

# rejoin the layers
sobj <- JoinLayers(sobj)

# SingleR
library(SingleR)
library(scater)
ref <- as.SingleCellExperiment(sobj_reference)
ref <- logNormCounts(ref)

test <- as.SingleCellExperiment(sobj)
test <- logNormCounts(test)
pred <- SingleR(test = test, 
               ref = ref, 
               labels = ref$celltype,
               de.method="wilcox")
sobj$celltype <- pred$labels


# Find markers for each cluster
sobj.markers <- FindAllMarkers(sobj,only.pos = TRUE,min.pct = 0.25,logfc.threshold =0.25 )
sobj.markers %>%
  group_by(cluster) %>%
  slice_max(n=2,order_by = avg_log2FC)
# VlnPlot(sobj, features = c("Gpc6", "Gata4"))
# p1 <- FeaturePlot(sobj, features = c("Gpc6", "Gata4"))
# saveRDS(sobj, file="E9_5_to_18_5.Rds")


######################remove doublets#################
# split the RNA measurements into several layers
sobj[["RNA"]] <- split(sobj[["RNA"]], f = sobj$sample_name)

sobj <- subset(sobj,subset=nFeature_RNA<3200)
# cells

sobj <- NormalizeData(sobj)
sobj <- FindVariableFeatures(sobj,selection.method = "vst",nfeatures = 2000)
sobj <- ScaleData(sobj)
sobj <- RunPCA(sobj,features = VariableFeatures(object=sobj))

sobj <- FindNeighbors(sobj,dims = 1:30,reduction = "pca")
sobj <- FindClusters(sobj,resolution = 0.1,cluster.name = "unintegrated_clusters")
sobj <- RunUMAP(sobj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


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
sobj <- FindClusters(sobj, resolution = 2, cluster.name = "rpca_clusters")
sobj <- RunUMAP(sobj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:30)
sobj <- FindClusters(sobj, resolution = 2, cluster.name = "har_clusters")
sobj <- RunUMAP(sobj, reduction = "harmony", dims = 1:30, reduction.name = "umap.har")


sobj <- FindNeighbors(sobj, reduction = "integrated.mnn", dims = 1:30)
sobj <- FindClusters(sobj, resolution = 2, cluster.name = "mnn_clusters")
sobj <- RunUMAP(sobj, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")


saveRDS(sobj, file="E9_5_to_18_5_rmDoublets_nfeature.Rds")