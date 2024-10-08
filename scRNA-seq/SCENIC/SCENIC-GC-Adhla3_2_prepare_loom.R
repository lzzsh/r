library(RColorBrewer)
library(SCENIC)
library(doParallel)
library(dplyr)

cores <- 15

exprMat <- readRDS(file = "~/exprMat_TB.Rds")
cellInfo <- readRDS(file = "~/cellInfo_TB.Rds")
dim(exprMat)
# Cell info/phenodata
head(cellInfo)
str(cellInfo)
unique(cellInfo$CellType)
cellInfo$CellType <- as.factor(cellInfo$CellType)
cellInfo$CellType <- droplevels(cellInfo$CellType)
unique(cellInfo$CellType)

dir.create("int")
saveRDS(cellInfo, file="~/int/cellInfo.Rds")

colVars <- list(CellType = c(
  "LaTP" = "#FF6347",              # Tomato
  "SynTI" = "#4682B4",             # Steel Blue
  "S-TGC precursor" = "#32CD32",   # Lime Green
  "JZP" = "#FFD700",               # Gold
  "P-TGC" = "#DA70D6",             # Orchid
  "SpT precursor" = "#F08080",     # Light Coral
  "SynTII" = "#8A2BE2",            # Blue Violet
  "SynTII precursor" = "#20B2AA",  # Light Sea Green
  "SynTI precursor" = "#FFA07A",   # Light Salmon
  "LaTP2" = "#87CEFA",             # Light Sky Blue
  "GC precursor" = "#B22222",      # Firebrick
  "S-TGC" = "#FF69B4",             # Hot Pink
  "GC-Prl7b1" = "#FFA500",         # Orange
  "SpT-outer" = "#008080",         # Teal
  "GC-Aldh1a3" = "#808000",        # Olive
  "SpT-inner" = "#6A5ACD"          # Slate Blue
))

saveRDS(colVars, file="~/int/colVars.Rds")

org <- "mgi" # or hgnc, or dmel
dbDir <- "/storage/liuxiaodongLab/liaozizhuo/data/cisTarget_databases/mm10"
myDatasetTitle <- "SCENIC" # choose a name for your analysis
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
data(defaultDbNames)
defaultDbNames$mgi[1] <- "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
defaultDbNames$mgi[2] <- "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
defaultDbNames
scenicOptions <- initializeScenic(org = org,
                                  dbDir = dbDir,
                                  dbs = defaultDbNames[["mgi"]],
                                  datasetTitle = myDatasetTitle,
                                  nCores = cores)
scenicOptions@settings$dbs <- defaultDbNames[["mgi"]]
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$seed <- 123

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "~/int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "~/int/colVars.Rds"
# Save to use at a later time...
saveRDS(scenicOptions, file="~/int/scenicOptions.Rds")

### Co-expression network
exprMat <- as.matrix(exprMat)
genesKept <- geneFiltering(exprMat,
                           scenicOptions = scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01
)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
saveRDS(exprMat_filtered, file = "~/int/exprMat_filtered.rds")

loom <- build_loom("~/exprMat_filtered_TB.loom", dgem = exprMat_filtered)
close_loom(loom)

# go out and run GRNBoost
