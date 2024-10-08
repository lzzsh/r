library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stats)
library(data.table)
library(dplyr)
library(openxlsx)

all_markers <- read.xlsx("DEG_TB.xlsx",sheet = "GC-Aldh1a3")

# Loop through each unique cluster in the data
unique_clusters <- unique(all_markers$cluster)

for (cluster_name in unique_clusters) {
DEG_data <- all_markers %>% filter(cluster == cluster_name)

# Gene name to GeneID conversion
gene_df <- bitr(DEG_data$gene, fromType = "SYMBOL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db)

colnames(gene_df)[1] <- "gene"
DEG_data1 <- left_join(gene_df, DEG_data)

# GO enrichment
GO_all <- enrichGO(gene = DEG_data1$ENTREZID,
                   keyType = "ENTREZID",
                   OrgDb = org.Mm.eg.db,
                   ont = "ALL",
                   pvalueCutoff = 0.01,
                   pAdjustMethod = "fdr",
                   minGSSize = 10,
                   maxGSSize = 500,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

GO_result <- data.frame(GO_all)

# Write GO_result to Excel
write.xlsx(GO_result, paste0("output/", cluster_name, "_GO_result.xlsx"))

# Top 12 significant GO terms for plotting
go_enrichment_pathway <- GO_result %>% group_by(ONTOLOGY) %>% top_n(n = 12, wt = -p.adjust)

# Create ggplot
p <- ggplot(go_enrichment_pathway, aes(x=reorder(Description, -Count), y=Count, fill=ONTOLOGY)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.6) +
  theme_minimal() +
  labs(x="GO Term", y="Gene_Number", title="Top 12 Enriched GO Terms") +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y') +
  coord_flip() +
  scale_fill_manual(values=c("CC"="skyblue","BP"="pink","MF"="lightgreen")) + 
  theme_bw()

# Save ggplot to PDF
pdf_file <- paste0("plots/", cluster_name, "_GO_plot.pdf")
pdf(pdf_file)
print(p)
dev.off()
}

# Selected GO Terms
GC.A <- read.xlsx("./output/GC-Aldh1a3_GO_result.xlsx")
GC.P <- read.xlsx("./output/GC-Prl7b1_GO_result.xlsx")
selected_description <- c("glycogen metabolic process",
                          "glucose homeostasis",
                          "glycogen biosynthetic process",
                          "actin filament organization",
                          "small GTPase mediated signal transduction",
                          "regulation of GTPase activity",
                          "chromatin remodeling",
                          "regulation of vasculature development",
                          "artery development",
                          "artery morphogenesis",
                          "epithelial cell migration")

GC.A <- GC.A %>%
  filter(Description %in% selected_description) %>%
  mutate(Gene = "GC-1")
GC.P <- GC.P %>%
  filter(Description %in% selected_description) %>%
  mutate(Gene = "GC-2")
plot <- rbind(GC.P,GC.A)
plot$p.adjust <- -10*(log(10,plot$p.adjust))

# GO plot
p <- ggplot(plot, aes(x = Description, y = Gene, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, stroke = 1) +  # shape 21 是一个有边框的圆
  theme_minimal() +
  coord_flip() +
  theme_bw() + 
  scale_fill_gradient(low = "#591877", high = "#F1E575") + 
  scale_size(range = c(2, 10)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )
p
ggsave("DEG-GO-bubble.pdf",p,device = "pdf",width=8,height=8)

