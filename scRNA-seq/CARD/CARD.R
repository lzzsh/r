devtools::install_github('YingMa0107/CARD') #安装CARD
devtools::install_github('xuranw/MuSiC') #安装依赖包
 
###加载CARD包
library(CARD)
library(MuSiC)
library(Seurat)
library(patchwork)
library(tidyverse)
 
#载入实战2保存的 GBM4 空转数据
load('GBM4.rdata')
 
###获取空转的counts表达矩阵
spatial_count <-  GBM4@assays$Spatial@counts
spatial_count[1:4,1:4]
 
###获取空转的空间位置矩阵
spatial_loca <- GBM4@images$GBM4@coordinates
spatial_location <- spatial_loca[,2:3]
#名字必须是x y ，不然后面CARD_deconvolution会报错
colnames(spatial_location) <- c("x","y")
spatial_location[1:4,]
 
#                   x   y
#AAACAAGTATCTCCCA-1 50 102
#AAACACCAATAACTGC-1 59  19
#AAACAGAGCGACTCCT-1 14  94
#AAACAGGGTCTATATT-1 47  13
 
 
#载入实战3保存的 GBM-scRNA 单细胞数据
load('scRNA.rdata')
 
#获取单细胞counts矩阵
sc_count <- scRNA@assays$RNA@counts
 
#获取单细胞细胞注释矩阵
sc_meta <- scRNA@meta.data %>% 
  rownames_to_column("cellID") %>%
  dplyr::select(cellID,orig.ident,celltype) %>% 
  mutate(CB = cellID) %>% 
  column_to_rownames("CB")
head(sc_meta)
 
head(sc_meta)
#                                                 cellID orig.ident   celltype
#GSM4119531_AAACCTGAGTCAAGGC GSM4119531_AAACCTGAGTCAAGGC GSM4119531   MES like
#GSM4119531_AAACCTGTCAGGCAAG GSM4119531_AAACCTGTCAGGCAAG GSM4119531 Macrophage
#GSM4119531_AAACCTGTCCTGCCAT GSM4119531_AAACCTGTCCTGCCAT GSM4119531   OPC like
 
 
##构建CARD对象
CARD_obj = createCARDObject( 
  sc_count = sc_count, 
  sc_meta = sc_meta, 
  spatial_count = spatial_count, 
  spatial_location = spatial_location, 
  ct.varname = "celltype", 
  ct.select = unique(sc_meta$celltype), #细胞类型列名
  sample.varname = "orig.ident")
 
## QC on scRNASeq dataset! ...
## QC on spatially-resolved dataset! ...
 
#CARD 解卷积
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
## create reference matrix from scRNASeq...
## Select Informative Genes! ...
## Deconvolution Starts! ...
## Deconvolution Finish! ...
 
#CARD-spot 可视化spot的细胞类型分布饼图
colors = c("#4DAF4A","#F0027F","#377EB8","#FDC086","#A6761D","#FFFF00","#BEAED4",
           "#BF5B17","#666666")
p1<-CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD,
                         spatial_location = CARD_obj@spatial_location, 
                         colors = colors)
 
p2<-SpatialPlot(GBM4,group.by = 'Region',cols = c('Normal'='#007799','Tumor'='#AA0000'))
p1+p2


#选择一些感兴趣的细胞类型分别进行可视化
ct.visualize = c("OPC like","AC like","MES like")
 
p3 <- CARD.visualize.prop(proportion = CARD_obj@Proportion_CARD,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = ct.visualize,                
  colors = c("lightblue","lightyellow","red"), 
  NumCols = 3,pointSize = 1)#图中spot大小
p3

#原始的基因表达
p4 <- CARD.visualize.gene(
  spatial_expression = CARD_obj@spatial_countMat,
  spatial_location = CARD_obj@spatial_location,
  gene.visualize = c("SNAP25","SYT1","JUNB"),
  colors = NULL,
  NumCols =3)
p4
 
#保存 CARD 反卷积结果
save(CARD_obj,file = 'CARD_obj.rdata')
