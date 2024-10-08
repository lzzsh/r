library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(spacexr)

# Load rds file
sobj <- readRDS("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/spatial/Spatial_16.5_RCTD.Rds")
RCTD_Score <- readRDS("/storage/liuxiaodongLab/liaozizhuo/Projects/YutingFu/Placenta_Project/result/spatial/Spatial_16.5_RCTD_Score.Rds")

# Rotate images
x_locations <- sobj@images[[1]]@coordinates$imagerow
y_locations <- sobj@images[[1]]@coordinates$imagecol
temp_x <- -(x_locations)
temp_y <- y_locations

# Add reduction object
cell_position <- as.matrix(cbind(temp_x,temp_y))
colnames(cell_position) <- c("Spatial_1","Spatial_2")
rownames(cell_position) <- colnames(sobj)
sobj[["Spatial_location"]] <- CreateDimReducObject(embeddings = cell_position, key = "Spatial_", assay = "RNA")

# 
umap_data <- as.data.frame(Embeddings(sobj, reduction = "Spatial_location"))
angle <-  -190 * pi / 180  
rotation_matrix  <- matrix(c(cos(angle), -sin(angle),
                            sin(angle), cos(angle)), nrow = 2, byrow = TRUE)


umap_data  <- as.data.frame(as.matrix(umap_data) %*% rotation_matrix)
colnames(umap_data) <- c("Spatial_1","Spatial_2")
umap_data$celltype <- sobj$first_type


# Set colors for all cell types
colors <- c(
    "#669933",    "#339966",    "#339999",    "#85C247",    "#A3D175",    "#BDC25B",
    "#AC80D0",    "#D0AC80",    "#D080CC",    "#D08380",    "#D3D3D3",    "#D080A4",
    "#A6DBDE",    "#E633D1",    "#EC7969",    "#E63377",    "#E64933",    "#47A3FF",
    "#06D6A0",    "#FFFF0A",    "#FF8667",    "#E7C19C",    "#EEB1AA",    "#F5934E",
    "#DB8B0A",    "#A30059",    "#4AC170",    "#7BBC54",    "#FFB3FF",    "#BDBD2E",
    "#FFA64D",    "#CDCD32",    "#CD8032",    "#BD752E",    "#DDA873")

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
                   "NK")
umap_data$celltype <- factor(umap_data$celltype,cluster_names)
names(colors) <- cluster_names

p <- ggplot(umap_data, aes(x = Spatial_1, y = Spatial_2, color = celltype)) +
  geom_point(size = 0.5) +
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
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 16)  # Increase legend title size
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # Increase legend dot size
p

# 
# highlight_celltype <- "Angiogenic DSC" #
# umap_data$color <- ifelse(umap_data$celltype == highlight_celltype, 
#                           colors[match(umap_data$celltype, names(colors))], 
#                           "#484846")

# Plot grouped by cell types
# p <- ggplot(umap_data, aes(x = Spatial_1, y = Spatial_2, color = color)) +
#   geom_point(size = 0.5) +
#   scale_color_identity() + 
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
#     legend.position = "none" 
#   )
# p
# 
# 
# ggsave("Spatial_DimPlot_RCTD_Angiogenic DSC_E9_5.pdf",p,device = "pdf",width=12,height=8)
# 
# 
# library(qpdf)
# 
# # 假设你有三个PDF文件需要合并
# pdf_files <- c("Spatial_DimPlot_RCTD_Angiogenic DSC_E9_5.pdf", "Spatial_DimPlot_RCTD_Angiogenic DSC_E10_5.pdf",
#                "Spatial_DimPlot_RCTD_Angiogenic DSC_E11_5.pdf","Spatial_DimPlot_RCTD_Angiogenic DSC_E12_5.pdf",
#                "Spatial_DimPlot_RCTD_Angiogenic DSC_E13_5.pdf","Spatial_DimPlot_RCTD_Angiogenic DSC_E14_5.pdf",
#                "Spatial_DimPlot_RCTD_Angiogenic DSC_E15_5.pdf","Spatial_DimPlot_RCTD_Angiogenic DSC_E16_5.pdf",
#                "Spatial_DimPlot_RCTD_Angiogenic DSC_E18_5.pdf")
# 
# # 创建一个新的PDF文件来保存合并后的结果
# output_file <- "merged_output.pdf"
# 
# # 合并PDF文件
# pdf_combine(pdf_files, output = output_file)


# Make plots
results <- RCTD_Score@results
norm_weights <- normalize_weights(results$weights)
cell_type_names <- RCTD_Score@cell_type_info$info[[2]]
spatialRNA <- RCTD_Score@spatialRNA
spatialRNA@coords$x <- umap_data$Spatial_1
spatialRNA@coords$y <- umap_data$Spatial_2

plot_puck_wrapper <- function(puck, plot_val, cell_type = NULL, minUMI = 0, maxUMI = 200000, min_val = NULL, max_val = NULL, title = NULL, my_cond = NULL) {
  UMI_filter = (puck@nUMI >= minUMI) & (puck@nUMI < maxUMI)
  ylimit = NULL
  if(!is.null(my_cond))
    my_cond = UMI_filter & my_cond
  else
    my_cond = UMI_filter
  if(!is.null(cell_type))
    my_cond = my_cond & (puck@cell_labels == cell_type)
  if(!is.null(min_val))
    my_cond = my_cond & (plot_val > min_val)
  if(!is.null(max_val)) {
    epsilon = 0.00001
    plot_val[plot_val >= max_val - epsilon] = max_val - epsilon
    if(!is.null(min_val))
      ylimit = c(min_val, max_val)
  }
  
p <- plot_puck_continuous(puck, names(which(my_cond)), plot_val, title = title, ylimit = ylimit)

# Add the color scale
p <- p +
  ggplot2::scale_color_gradientn(colors = c("blue", "yellow", "red"), limits = c(0, 1))

# Add theme to change the background color to black
p <- p +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "black"),
    plot.background = ggplot2::element_rect(fill = "black"),
    plot.title = ggplot2::element_text(color = "white", size = 36)
    # text = ggplot2::element_text(color = "white"),
    # axis.title.x=element_blank(),
    # axis.title.y=element_blank(),
    # axis.ticks=element_blank(),
    # axis.text.x=element_blank(),
    # axis.text.y=element_blank(),
    # legend.position="none"
  )

return(p)
}

plots <- vector(mode = "list", length = length(cell_type_names))
weights <- norm_weights
puck <- spatialRNA

for (i in 1:length(cell_type_names)) {
  cell_type = cell_type_names[i]
  plot_var <- weights[,cell_type]; names(plot_var) = rownames(weights)
  if(sum(weights[,cell_type]) > 0)
    plots[[i]] <- plot_puck_wrapper(puck, plot_var, NULL, minUMI = 0, maxUMI = 200000, min_val = 0, max_val = 1, title = cell_type)
}

pdf("/storage/liuxiaodongLab/liaozizhuo/plots/RCTD_E16.5.pdf", width = 20, height = 16)
invisible(lapply(plots, print))
dev.off()
