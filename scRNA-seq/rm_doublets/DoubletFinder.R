library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(batchelor)
library(SeuratWrappers)

options(future.globals.maxSize = 20 * 1024^3)  


sobj.data<-Read10X(data.dir = "/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/aggr_result/AGG_9_to_18/outs/count/filtered_feature_bc_matrix")
sobj<-CreateSeuratObject(counts = sobj.data,project = "sobj3k",min.cells = 3,min.features = 200)
# 138202 cells

sobj[['percent.mt']]<-PercentageFeatureSet(sobj,pattern = "^mt-")

# sobj<-subset(sobj,subset=nFeature_RNA>500 & percent.mt<10)
# 121752 cells

sobj<-subset(sobj,subset=nFeature_RNA>500 & percent.mt<10)
# cells

sobj<-NormalizeData(sobj)
sobj<-FindVariableFeatures(sobj,selection.method = "vst",nfeatures = 2000)
sobj<-ScaleData(sobj)
sobj<-RunPCA(sobj,features = VariableFeatures(object=sobj))
sobj<-FindNeighbors(sobj,dims = 1:30,reduction = "pca")
sobj<-FindClusters(sobj,resolution = 0.1,cluster.name = "unintegrated_clusters")
sobj<-RunUMAP(sobj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

## (2)pK Identification ----------------------------------------------------------
sweep.res.list_sobj <- paramSweep(sobj, PCs = 1:30, sct = FALSE)
#head(sweep.res.list_sobj)
sweep.stats_sobj <- summarizeSweep(sweep.res.list_sobj, GT = FALSE)
bcmvn_sobj <- find.pK(sweep.stats_sobj) 
mpK<-as.numeric(as.vector(bcmvn_sobj$pK[which.max(bcmvn_sobj$BCmetric)]))


## (3) Homotypic Doublet Proportion Estimate -------------------------------------
annotations <- sobj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = 0.05

nExp_poi <- round(DoubletRate*length(sobj$seurat_clusters)) 

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## (4) Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
sobj <- doubletFinder(sobj, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
sobj <- doubletFinder(sobj, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)

## Plot results ---------------------------------------------------------------------------
saveRDS(sobj, file = "/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/integration/E9_5_to_18_5_rmDoublets_DoubleFinder1.Rds")

sobj@meta.data$DF_hi.lo <- sobj@meta.data$DF.classifications_0.25_0.23_118588
sobj@meta.data$DF_hi.lo[which(sobj@meta.data$DF_hi.lo == "Doublet" & sobj@meta.data$DF.classifications_0.25_0.23_109060 == "Singlet")] <- "Doublet-Low Confidience"
sobj@meta.data$DF_hi.lo[which(sobj@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
table(sobj@meta.data$DF_hi.lo)

## Visualization
DimPlot(sobj, reduction = "umap", group.by ="DF_hi.lo",cols =c("black","red","gold"))

## remove doublets
sobj<-subset(sobj,subset=DF_hi.lo == "Singlet")

saveRDS(sobj, file = "/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/integration/E9_5_to_18_5_rmDoublets_DoubleFinder2.Rds")

######### re-run UMAP ###########
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(batchelor)
library(SeuratWrappers)
sobj <- readRDS("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/integration/E9_5_to_18_5_rmDoublets_DoubleFinder2.Rds")
# split the RNA measurements into several layers
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

saveRDS(sobj, file = "/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/integration/E9_5_to_18_5_rmDoublets_DoubleFinder3.Rds")