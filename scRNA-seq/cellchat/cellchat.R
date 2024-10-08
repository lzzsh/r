library(CellChat)
library(patchwork)
setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/GC_subtype/")
options(future.globals.maxSize = 3 * 1024^3)  # Set limit to 2 GiB

data.input <- readRDS("./cellchat//cellchat_data.input.rds")
meta <- readRDS("./cellchat//cellchat_meta.rds")
unique(meta$labels)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("./plots/net_circle.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_circle(cellchat@net$count, 
                 idents.use = c(1:2),
                 vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:15), remove.isolate = FALSE)

mat <- cellchat@net$weight
# Loop through the first 29 row names of the matrix
for (i in 1:15) {
  # Get the ith row name for the PDF file name
  pdf_file_name <- rownames(mat)[i]
  
  # Construct the full filename with the .pdf extension
  full_pdf_file_name <- paste0("./plots_v2/", pdf_file_name, ".pdf")
  
  # Create the bubble plot for the ith source
  plot <- netVisual_bubble(cellchat, sources.use = i, targets.use = 1:15, remove.isolate = FALSE)
  
  # Close the PDF device
  ggsave(filename = full_pdf_file_name, plot = plot, device = "pdf", width = 10, height = 18)
}

a <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,3,10,11), signaling = c("IGF"), remove.isolate = FALSE)
a <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:15), signaling = c("IGF","CCL","SEMA3"), remove.isolate = FALSE)
ggsave(filename = "cellchat_GC2_bubble.pdf", plot = a, device = "pdf", width = 10, height = 12)

netVisual_chord_gene(cellchat, sources.use = 3, targets.use = c(1:15), lab.cex = 0.01,legend.pos.y = 120)
