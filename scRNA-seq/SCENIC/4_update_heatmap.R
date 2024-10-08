library(RColorBrewer)
library(SCENIC)
library(doParallel)
library(dplyr)
library(AUCell)
library(ComplexHeatmap)

setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/scenic/Prl7b1/")
scenicOptions <- readRDS("int/scenicOptions.Rds")

cellInfo <- readRDS("int/cellInfo.Rds")

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

new_order <- c("JZP", "GC precursor", "GC-Aldh1a3", "GC-Prl7b1")
setdiff(new_order, colnames(regulonActivity_byCellType_Scaled))
setdiff(colnames(regulonActivity_byCellType_Scaled), new_order)
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[ , new_order]

pdf(file = "SCENIC_GC-Aldh1a3_reorder.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 48) # The height of the plot in inches

Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
        cluster_columns = FALSE)

dev.off()

row.names(regulonActivity_byCellType_Scaled)

filter_rows <- c( "Prdm1 (16g)", "Klf6 (24g)", "Creb5_extended (50g)", "Fosl2_extended (34g)", "Hmgb1 (15g)", "Nr3c1 (195g)", "Foxp1 (151g)", "Foxo4 (88g)")

regulonActivity_byCellType_Scaled_filter <- regulonActivity_byCellType_Scaled[filter_rows, ]

pdf(file = "SCENIC-GC-Aldh1a3.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 10) # The height of the plot in inches

Heatmap(regulonActivity_byCellType_Scaled_filter, name="Regulon activity",
        cluster_columns = FALSE)

dev.off()

regulon <- readRDS("./int/2.6_regulons_asGeneSet.Rds")