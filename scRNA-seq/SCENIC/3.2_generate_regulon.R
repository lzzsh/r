library(RColorBrewer)
library(SCENIC)
library(doParallel)
library(dplyr)

setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/scenic/")
scenicOptions <- readRDS("./int/scenicOptions.Rds")

GRNBoost_output <- read.delim("adj.tsv")
colnames(GRNBoost_output) <- c("TF", "Target", "weight")
saveRDS(GRNBoost_output, file="./int/1.4_GENIE3_linkList.Rds")

exprMat_filtered <- readRDS("./int/exprMat_filtered.rds")

runCorrelation(exprMat_filtered, scenicOptions)

runSCENIC_1_coexNetwork2modules(scenicOptions)

data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9

runSCENIC_2_createRegulons(scenicOptions)

scenicOptions@settings$nCores <- 1

exprMat <- readRDS(file = "./exprMat_TB.Rds")
exprMat <- as.matrix(exprMat)
exprMat_log <- log2(exprMat + 1)
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions)

cellInfo <- readRDS("./int/cellInfo.Rds")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]

regulonSelection <- loadInt(scenicOptions, "aucell_regulonSelection", ifNotExists="null", verbose=FALSE)
cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
cellInfo <- data.frame(cellInfo)
colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")

binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_full")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

pdf(file = "SCENIC_regulonActivity_byCellType_Scaled.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 48) # The height of the plot in inches

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")

dev.off()
