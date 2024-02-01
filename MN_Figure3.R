# The following script is for generating Figure 3 for Masako Nakanishi Gastro paper

# Set working directory
setwd("/users/michaelmartinez/Desktop/MN_Figure_3/")

# Set relative path to files
filePath <- file.path("/users/michaelmartinez/Desktop/AJB6_Ulceration/")

# Load relevant libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(pheatmap)
library(circlize)
library(clusterProfiler)
library(DOSE)
library(AnnotationDbi)
library(enrichplot)
library(org.Mm.eg.db)
library(gridExtra)
library(RColorBrewer)
library(stringr)

# -------------------------------------------------------------------------------- #
#                               Figure 3A                                          #
# -------------------------------------------------------------------------------- #

# -----  B6D:KO relative to A/D:KO top 50 bottom 50 heat-map
# First, read in the KO vs KO dds results and get results and counts
KOvsKO <- readRDS(paste(filePath, "Renamed_Samples_Figures/Renamed_RDS_files/ADKO_vs_B6DKO.rds", sep = ""))
dds <- as.data.frame(results(KOvsKO))
counts <- as.data.frame(counts(KOvsKO, normalized = TRUE))
res <- merge(dds, counts, by = 0, all = TRUE)
res <- na.omit(res)

# Take significant results with alpha = 0.05 and |log2FC| > 2 and order by decreasing value
res.top <- res[abs(res$log2FoldChange) > 2 & res$padj < 0.05,]
res.ordered <- res[order(res$log2FoldChange, decreasing = TRUE),]

# Get symbols, only keep unique ones
res.ordered$Symbols <- gsub("^[^-]+-(.*)$", "\\1", res.ordered$Row.names)
res.ordered <- res.ordered[!duplicated(res.ordered$Symbols),]
rownames(res.ordered) <- res.ordered$Symbols
res.ordered$Row.names <- NULL
res.ordered$Symbols <- NULL

# Take the top 50 and top 100 (top/bottom 25 and 50 respectively). This is so that we have options
total_rows <- nrow(res.ordered)
top100 <- c(1:50, (total_rows - 49):total_rows)
top50 <- c(1:25, (total_rows - 24):total_rows)

# Subset the res.ordered df into two new dataframes for top 100 and top 50
top100_subset <- res.ordered[top100,]
top50_subset <- res.ordered[top50,]

# Create a design data frame
KOdesign <- data.frame(Genotype = rep(c("A/D:KO", "B6D:KO"),
                                      c(12,11)))

# Factor the levels so they appear in the correct order
KOdesign$Genotype <- factor(KOdesign$Genotype, levels = c("A/D:KO", "B6D:KO"))

# Plot the top100 subset
top100.degs.subset <- top100_subset[order(top100_subset$log2FoldChange, decreasing = TRUE),]

# Pull the baseMean and Log2FC column as a matrix
top100.FC <- as.matrix(top100.degs.subset$log2FoldChange)
colnames(top100.FC) <- "Log2FC"
top100.BM <- as.matrix(top100.degs.subset$baseMean)
colnames(top100.BM) <- "Base Mean"

# Extract just the first column and the counts columns
X100 <- top100.degs.subset[,c(7:ncol(top100.degs.subset))]

# Transpose, center, and scale the normalized counts
X100.scaled <- t(apply(X100, 1, scale))
colnames(X100.scaled) <- colnames(X100)

# Create a named color vector for the heat-map
#-----GROUPS: A/D:KO = red, A/D:WT = forestgreen, B6D:KO = blue, B6D:WT = purple
color_vector <- c("red", "blue")
names(color_vector) <- c("A/D:KO", "B6D:KO")

# Setting "4" for fontface specifies bold and italic
hmAnno <- HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  annotation_height = unit(4, "cm"),
  show_annotation_name = FALSE,
  Group = anno_block(gp = gpar(fill = color_vector), labels = c("A/D:KO", "B6D:KO"),
                     labels_gp = gpar(col = "white", fontface = 4, fontsize = 60)),
  col = list(Group = color_vector))

# Set heatmap splitting pattern
hmSplit <- rep(1:2, c(12, 11))

# Map colors to values (`library(circlize)` required)
X100.FC.colors <- colorRamp2(c(min(top100.FC),
                               0,
                               max(top100.FC)),
                             c("blue", "white", "red"))
X100.BM.colors <- colorRamp2(c(quantile(top100.BM)[1],
                               quantile(top100.BM)[4]),
                             c("white", "red"))

# Define the number of slices to show in the heatmap
slices <- 12+11

# Heatmap for top 100 scaled data
X100.hmScaled <- Heatmap(X100.scaled,
                         column_labels = colnames(X100.scaled),
                         name = "Z-Score",
                         cluster_rows = FALSE,
                         cluster_columns = TRUE,
                         column_dend_height = unit(4.5, "cm"),
                         column_dend_gp = gpar(lwd = c(4)),
                         column_order = order(factor(colnames(X100.scaled))),
                         top_annotation = hmAnno,
                         column_split = hmSplit,
                         show_column_names = FALSE,
                         heatmap_legend_param = list(labels_gp = gpar(fonsize = 16),
                                                     title_gp = gpar(fontsize = 14),
                                                     legend_direction = "horizontal",
                                                     legend_width = unit(3, "cm"),
                                                     legend_height = unit(4, "cm")))

# Heatmap for log2FC values
X100.hml2fc <- Heatmap(top100.FC,
                       row_labels = rownames(X100.scaled),
                       cluster_rows = FALSE,
                       name = "log2 FC",
                       column_names_rot = 90,
                       col = X100.FC.colors,
                       cell_fun = function(j,i, x, y, w, h, col) {
                         grid.text(round(top100.FC[i, j],1), x, y, 
                                   gp = gpar(fontsize = 10, 
                                             col = "black"))},
                       heatmap_legend_param = list(labels_gp = gpar(fonsize = 16),
                                                   title_gp = gpar(fontsize = 14),
                                                   legend_direction = "horizontal",
                                                   legend_width = unit(3, "cm"),
                                                   legend_height = unit(4, "cm")))

# Heatmap for average expression
X100.hmbm <- Heatmap(top100.BM,
                     row_labels = rownames(X100.scaled),
                     row_names_gp = gpar(fontsize = 10),
                     cluster_rows = FALSE,
                     name = "Avg Expression",
                     column_names_rot = 90,
                     col = X100.BM.colors,
                     cell_fun = function(j, i, x, y, w, h, col) {
                       grid.text(round(top100.BM[i, j],1), x, y,
                                 gp = gpar(fontsize = 10,
                                           col = "black"))},
                     heatmap_legend_param = list(labels_gp = gpar(fonsize = 16),
                                                 title_gp = gpar(fontsize = 14),
                                                 legend_direction = "horizontal",
                                                 legend_width = unit(3, "cm"),
                                                 legend_height = unit(4, "cm")))
# Draw the heatmap
tiff("Fig3A_KOvsKO_Top100_W24_H16.tiff", width = 12, height = 20 , pointsize = 20, units = "in", res=300)
X100.Heatmap <- X100.hmScaled + X100.hml2fc + X100.hmbm
draw(X100.Heatmap, heatmap_legend_side = "bottom", legend_gap = unit(1, "cm"))
dev.off()


#-----Repeat with top 50
# Plot the top50 subset
top50.degs.subset <- top50_subset[order(top50_subset$log2FoldChange, decreasing = TRUE),]

# Pull the baseMean and Log2FC column as a matrix
top50.FC <- as.matrix(top50.degs.subset$log2FoldChange)
colnames(top50.FC) <- "Log2FC"
top50.BM <- as.matrix(top50.degs.subset$baseMean)
colnames(top50.BM) <- "Base Mean"

# Extract just the first column and the counts columns
X50 <- top50.degs.subset[,c(7:ncol(top50.degs.subset))]

# Transpose, center, and scale the normalized counts
X50.scaled <- t(apply(X50, 1, scale))
colnames(X50.scaled) <- colnames(X50)

# Create a named color vector for the heat-map
#-----GROUPS: A/D:KO = red, A/D:WT = forestgreen, B6D:KO = blue, B6D:WT = purple
color_vector <- c("red", "blue")
names(color_vector) <- c("A/D:KO", "B6D:KO")

# Setting "4" for fontface specifies bold and italic
hmAnno <- HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  annotation_height = unit(4, "cm"),
  show_annotation_name = FALSE,
  Group = anno_block(gp = gpar(fill = color_vector), labels = c("A/D:KO", "B6D:KO"),
                     labels_gp = gpar(col = "white", fontface = 4, fontsize = 60)),
  col = list(Group = color_vector))

# Set heatmap splitting pattern
hmSplit <- rep(1:2, c(12, 11))

# Map colors to values (`library(circlize)` required)
X50.FC.colors <- colorRamp2(c(min(top50.FC),
                               0,
                               max(top50.FC)),
                             c("blue", "white", "red"))
X50.BM.colors <- colorRamp2(c(quantile(top50.BM)[1],
                               quantile(top50.BM)[4]),
                             c("white", "red"))

# Define the number of slices to show in the heatmap
slices <- 12+11

# Heatmap for top 100 scaled data
X50.hmScaled <- Heatmap(X50.scaled,
                         column_labels = colnames(X50.scaled),
                         name = "Z-Score",
                         cluster_rows = FALSE,
                         cluster_columns = TRUE,
                         column_dend_height = unit(4.5, "cm"),
                         column_dend_gp = gpar(lwd = c(4)),
                         column_order = order(factor(colnames(X50.scaled))),
                         top_annotation = hmAnno,
                         column_split = hmSplit,
                         show_column_names = FALSE,
                         heatmap_legend_param = list(labels_gp = gpar(fonsize = 16),
                                                     title_gp = gpar(fontsize = 14),
                                                     legend_direction = "horizontal",
                                                     legend_width = unit(3, "cm"),
                                                     legend_height = unit(4, "cm")))

# Heatmap for log2FC values
X50.hml2fc <- Heatmap(top50.FC,
                       row_labels = rownames(X50.scaled),
                       cluster_rows = FALSE,
                       name = "log2 FC",
                       column_names_rot = 90,
                       col = X50.FC.colors,
                       cell_fun = function(j,i, x, y, w, h, col) {
                         grid.text(round(top50.FC[i, j],1), x, y, 
                                   gp = gpar(fontsize = 14, 
                                             col = "black"))},
                       heatmap_legend_param = list(labels_gp = gpar(fonsize = 16),
                                                   title_gp = gpar(fontsize = 14),
                                                   legend_direction = "horizontal",
                                                   legend_width = unit(3, "cm"),
                                                   legend_height = unit(4, "cm")))

# Heatmap for average expression
X50.hmbm <- Heatmap(top50.BM,
                     row_labels = rownames(X50.scaled),
                     row_names_gp = gpar(fontsize = 16),
                     cluster_rows = FALSE,
                     name = "Avg Expression",
                     column_names_rot = 90,
                     col = X50.BM.colors,
                     cell_fun = function(j, i, x, y, w, h, col) {
                       grid.text(round(top50.BM[i, j],1), x, y,
                                 gp = gpar(fontsize = 12,
                                           col = "black"))},
                     heatmap_legend_param = list(labels_gp = gpar(fonsize = 16),
                                                 title_gp = gpar(fontsize = 14),
                                                 legend_direction = "horizontal",
                                                 legend_width = unit(3, "cm"),
                                                 legend_height = unit(4, "cm")))
# Draw the heatmap
tiff("Fig3A_KOvsKO_Top50_W24_H16.tiff", width = 12, height = 20 , pointsize = 18, units = "in", res=300)
X50.Heatmap <- X50.hmScaled + X50.hml2fc + X50.hmbm
draw(X50.Heatmap, heatmap_legend_side = "bottom", legend_gap = unit(1, "cm"))
dev.off()

# -------------------------------------------------------------------------------- #
#                               Figure 3B                                          #
# -------------------------------------------------------------------------------- #

# ----- Figure 3B: IG Genes heat-map
# From the ordered DEG list, extract the immunoglobulin genes 
res.ordered$Symbols <- rownames(res.ordered)
IGsubset <- res.ordered[grepl("^ Ighg*|^ Igha*|^ Ighd* |^ Ighe* |^ Ighm* |^ Ighv.*", res.ordered$Symbols),]
rownames(IGsubset) <- IGsubset$Symbols
IGsubset$Row.names <- NULL
IGsubset$Symbols <- NULL

# Run the consistent heatmap features above (annotation, split, slices, etc...)

# Pull the baseMean and Log2FC column as a matrix
IG.FC <- as.matrix(IGsubset$log2FoldChange)
colnames(IG.FC) <- "Log2FC"
IG.BM <- as.matrix(IGsubset$baseMean)
colnames(IG.BM) <- "Base Mean"

# Extract just the first column and the counts columns
XIG <- IGsubset[,c(7:ncol(IGsubset))]

# Transpose, center, and scale the normalized counts
XIG.scaled <- t(apply(XIG, 1, scale))
colnames(XIG.scaled) <- colnames(XIG)

# Map colors to values (`library(circlize)` required)
XIG.FC.colors <- colorRamp2(c(min(IG.FC),
                              0,
                              max(IG.FC)),
                            c("blue", "white", "red"))
XIG.BM.colors <- colorRamp2(c(quantile(IG.BM)[1],
                              quantile(IG.BM)[4]),
                            c("white", "red"))

# Heatmap for top 100 scaled data
XIG.hmScaled <- Heatmap(XIG.scaled,
                        column_labels = colnames(XIG.scaled),
                        name = "Z-Score",
                        cluster_rows = FALSE,
                        cluster_columns = TRUE,
                        column_dend_height = unit(4.5, "cm"),
                        column_dend_gp = gpar(lwd = c(4)),
                        column_order = order(factor(colnames(XIG.scaled))),
                        top_annotation = hmAnno,
                        column_split = hmSplit,
                        show_column_names = FALSE,
                        heatmap_legend_param = list(labels_gp = gpar(fonsize = 16),
                                                    title_gp = gpar(fontsize = 14),
                                                    legend_direction = "horizontal",
                                                    legend_width = unit(3, "cm"),
                                                    legend_height = unit(4, "cm")))

# Heatmap for log2FC values
XIG.hml2fc <- Heatmap(IG.FC,
                      row_labels = rownames(XIG.scaled),
                      cluster_rows = FALSE,
                      name = "log2 FC",
                      column_names_rot = 90,
                      col = XIG.FC.colors,
                      cell_fun = function(j,i, x, y, w, h, col) {
                        grid.text(round(IG.FC[i, j],1), x, y, 
                                  gp = gpar(fontsize = 14, 
                                            col = "black"))},
                      heatmap_legend_param = list(labels_gp = gpar(fonsize = 16),
                                                  title_gp = gpar(fontsize = 14),
                                                  legend_direction = "horizontal",
                                                  legend_width = unit(3, "cm"),
                                                  legend_height = unit(4, "cm")))

# Heatmap for average expression
XIG.hmbm <- Heatmap(IG.BM,
                    row_labels = rownames(XIG.scaled),
                    row_names_gp = gpar(fontsize = 16),
                    cluster_rows = FALSE,
                    name = "Avg Expression",
                    column_names_rot = 90,
                    col = XIG.BM.colors,
                    cell_fun = function(j, i, x, y, w, h, col) {
                      grid.text(round(IG.BM[i, j],1), x, y,
                                gp = gpar(fontsize = 12,
                                          col = "black"))},
                    heatmap_legend_param = list(labels_gp = gpar(fonsize = 16),
                                                title_gp = gpar(fontsize = 14),
                                                legend_direction = "horizontal",
                                                legend_width = unit(3, "cm"),
                                                legend_height = unit(4, "cm")))
# Draw the heatmap
tiff("Fig3B_KO_IGs_W24_H16.tiff", width = 12, height = 20 , pointsize = 18, units = "in", res=300)
XIG.Heatmap <- XIG.hmScaled + XIG.hml2fc + XIG.hmbm
draw(XIG.Heatmap, heatmap_legend_side = "bottom", legend_gap = unit(1, "cm"))
dev.off()

# -------------------------------------------------------------------------------- #
#                               Figures 3E-F                                       #
# -------------------------------------------------------------------------------- #


# ----- Biological Network plots
# These results are from the FINAL_GSEA_Object_KOvsKO.rds

# Function for extracting a legend (yingtools2)
gg.legend <- function(a.gplot) {
  #extract legend, so it can be used with grid.arrange or whatever
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Get the Ensembl and Entrez IDs
res.ordered$Ensembl <- gsub("^(.*?) - .*", "\\1", res.ordered$Row.names)
res.ordered$Entrez <- mapIds(org.Mm.eg.db, key = res.ordered$Ensembl,
                             column = "ENTREZID", keytype = "ENSEMBL",
                             multiVals = "first")
res.ordered <- res.ordered[order(res.ordered$log2FoldChange, decreasing = TRUE),]
res.ordered.genes <- res.ordered$log2FoldChange

#Assign Entrez IDs as names for the genes
names(res.ordered.genes) <- res.ordered$Entrez

#Remove duplicated Entrez IDs and their corresponding values
unique_entrez_genes <- names(res.ordered.genes[!duplicated(names(res.ordered.genes))])
unique_genes <- res.ordered.genes[unique_entrez_genes]
unique_genes <- sort(unique_genes, decreasing = TRUE)

# Load in the GSEA results
GSEA <- readRDS("/users/michaelmartinez/Desktop/FINAL_GSEA_Object_KOvsKO.rds")
 
# Select only core enrichment genes with an |log2FC| > 2
core.genes <- str_split(as.data.frame(GSEA)[,"core_enrichment"], "/")
my.selected.genes <- names(unique_genes[abs(unique_genes) > 2.5])

# For each GO category, only keep the core enriched genes that are in that category AND have been selected
filtered.core.genes <- sapply(lapply(core.genes, function(x) x[x %in% my.selected.genes]),
                              paste, collapse = "/")
GSEA@result$core_enrichment <- filtered.core.genes
gse.filtered <- setReadable(GSEA, "org.Mm.eg.db", "ENTREZID")


# Plot CNET for A/D:KO categories
ADCNET <- cnetplot(gse.filtered, node_label = "gene", foldChange = unique_genes, cex_label_gene = 2.3, 
                   colorEdge = TRUE, layout = "kk", circular = TRUE,
                   showCategory = c("antimicrobial humoral immune response mediated by antimicrobial peptide",
                                                     "antimicrobial humoral response",
                                                     "defense response to Gram-positive bacterium",
                                                     "immunoglobulin production",
                                                     "leukocyte migration",
                                                     "leukocyte proliferation",
                                                     "mononuclear cell proliferation",
                                                     "production of molecular mediator of immune response",
                                                     "regulation of leukocyte proliferation",
                                                     "response to molecule of bacterial origin")) +
  theme(legend.position = "right",
        #legend.direction = "horizontal",
        legend.key.size = unit(1, "cm"), 
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 20)) +
  scale_color_gradient(low = "darkblue", high = "lightblue", name = "log2 FC")

# # Remove the legend
# AD_new_plot <- ADCNET + theme(legend.position = "none")
# 
# # Extract the legend from the original plot and combine to bottom
# AD_leg <- gg.legend(ADCNET)
# Final_ADCNET <- grid.arrange(AD_new_plot, bottom = AD_leg)

ggsave("Figure3E_AD_KOvsKO_CNET.tiff", ADCNET, width = 25, height = 24, dpi = 300)



# Reset the GSEA object and pick B6D categories now
GSEA <- readRDS("/users/michaelmartinez/Desktop/FINAL_GSEA_Object_KOvsKO.rds")

# For each GO category, only keep the core enriched genes that are in that category AND have been selected
core.genes <- str_split(as.data.frame(GSEA)[,"core_enrichment"], "/")
my.selected.genes <- names(unique_genes[abs(unique_genes) > 2.0])

#For each GO cateogry, only keep the core enriched genes that are in that category AND have been selected
filtered.core.genes <- sapply(lapply(core.genes, function(x) x[x %in% my.selected.genes]),
                              paste, collapse = "/")
GSEA@result$core_enrichment <- filtered.core.genes
gse.filtered <- setReadable(GSEA, "org.Mm.eg.db", "ENTREZID")

# Plot
B6DCNET <- cnetplot(gse.filtered, node_label = "gene", foldChange = unique_genes, cex_label_category = 0.5, cex_label_gene = 4.5, 
                    colorEdge = TRUE, 
                    circular = TRUE, showCategory = c("asymmetric synapse",
                                                      "channel activity",
                                                      "glutamatergic synapse",
                                                      "integral component of synaptic membrane",
                                                      "intrinsic component of synaptic membrane",
                                                      "ion channel activity",
                                                      "neuron to neuron synapse",
                                                      "passive transmembrane transporter activity",
                                                      "suppression of viral release by host",
                                                      "Z disc")) +
  theme(legend.position = "right",
        #legend.direction = "horizontal",
        legend.key.size = unit(1, "cm"), 
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 20)) +
  scale_color_gradient(low = "lightblue", high = "red", name = "log2 FC")

# Remove the legend 
# B6D_new_plot <- B6DCNET + theme(legend.position = "none")
# 
# # Extract the legend from the original plot and combine to bottom
# B6D_leg <- gg.legend(B6DCNET)
# Final_B6DCNET <- grid.arrange(B6D_new_plot, bottom = B6D_leg)
ggsave("Figure3F_B6D_KOvsKO_CNET.tiff", B6DCNET, width = 25, height = 24, dpi = 300)


# Lastly, read in immunological dataset
C7 <- readRDS("/users/michaelmartinez/Desktop/AJB6_Ulceration/ImmunologicalSignatures_GSEA.rds")
C7.readable <- setReadable(C7, org.Mm.eg.db, "ENTREZID")

core.genes <- str_split(as.data.frame(C7)[,"core_enrichment"], "/")
my.selected.genes <- names(unique_genes[abs(unique_genes) > 1.0])

#For each GO cateogry, only keep the core enriched genes that are in that category AND have been selected
filtered.core.genes <- sapply(lapply(core.genes, function(x) x[x %in% my.selected.genes]),
                              paste, collapse = "/")
C7@result$core_enrichment <- filtered.core.genes
C7.filtered <- setReadable(C7, "org.Mm.eg.db", "ENTREZID")
C7.df <- as.data.frame(C7.filtered)

# Plot
C7CNET <- cnetplot(C7.filtered, node_label = "gene", foldChange = unique_genes, cex_label_category = 0.5, cex_label_gene = 3.0, 
                    colorEdge = TRUE, layout = "kk", 
                    circular = FALSE, showCategory = c("GSE10325_CD4_TCELL_VS_BCELL_DN",
                                                       "FLETCHER_PBMC_BCG_10W_INFANT_BCG_STIMULATED_VS_UNSTIMULATED_10W_UP",
                                                       "FLETCHER_PBMC_BCG_10W_INFANT_PPD_STIMULATED_VS_UNSTIMULATED_10W_UP",
                                                       "GSE10325_LUPUS_CD4_TCELL_VS_LUPUS_BCELL_DN",
                                                       "GSE18281_SUBCAPSULAR_VS_CENTRAL_CORTICAL_REGION_OF_THYMUS_DN",
                                                       "GSE21360_SECONDARY_VS_QUATERNARY_MEMORY_CD8_TCELL_UP",
                                                       "GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY3_DN",
                                                       "GSE3039_NKT_CELL_VS_ALPHAALPHA_CD8_TCELL_UP",
                                                       "GSE45365_NK_CELL_VS_CD11B_DC_DN",
                                                       "MATSUMIYA_PBMC_MODIFIED_VACCINIA_ANKARA_VACCINE_AGE_18_55YO_2DY_UP")) +
  theme(legend.position = "right",
        legend.direction = "horizontal",
        legend.key.size = unit(1, "cm"), 
        legend.text = element_text(size = 8.5),
        legend.title = element_text(size = 12)) +
  scale_color_gradient(low = "darkblue", high = "lightblue", name = "log2 FC")

# Remove the legend 
# C7_new_plot <- C7CNET + theme(legend.position = "none")
# 
# # Extract the legend from the original plot and combine to bottom
# C7_leg <- gg.legend(C7CNET)
# Final_C7CNET <- grid.arrange(C7_new_plot, bottom = C7_leg)
ggsave("Figure3H_B6D_KOvsKO_CNET.tiff", C7CNET, width = 25, height = 24, dpi = 300)


# -------------------------------------------------------------------------------- #
#                               Figures 3C                                         #
# -------------------------------------------------------------------------------- #


# KEGG dotplot
keggdf <- read.csv("/users/michaelmartinez/Desktop/keggdf.csv", header = TRUE, sep = ",")
keggdf$Sign <- ifelse(keggdf$NES < 0, "Enriched in A/D:KO", "Enriched in B6D:KO")
keggdf$Description <- str_wrap(keggdf$Description, width = 30)
keggdf <- keggdf[order(keggdf$NES, decreasing = FALSE),]
keggdf$Description <- factor(keggdf$Description, levels = keggdf$Description)

# Plot dotplot
KEGGDot <- ggplot(keggdf, aes(x = NES, y = Description, color = p.adjust)) +
  geom_point(aes(size = setSize),alpha = 0.8) +
  theme_bw() +
  scale_size(range = c(3.5,10)) +
  scale_color_continuous(low = "red", high = "blue", name = "P adj") +
  facet_grid(~Sign, scales = "free") +
  xlab("NES") +
  theme(axis.text.y = element_text(size = 9)) +
  theme(strip.text = element_text(face = "bold", size = 30),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24))
KEGGDot
ggsave("Keggdotplot.tiff", KEGGDot, width = 16, height = 24, dpi = 300)

# -------------------------------------------------------------------------------- #
#                               Figures 3D                                         #
# -------------------------------------------------------------------------------- #


# ontology plot
gsea <- readRDS("/Users/michaelmartinez/Desktop/FINAL_GSEA_Object_KOvsKO.rds")
KOgsea <- as.data.frame(setReadable(gsea, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"))
gse.readable <- KOgsea

#Group the results based on BP ontology and separate into positive and negative, arrange by decreasing NES
BP_terms <- gse.readable[gse.readable$ONTOLOGY == "BP", ]
posBP_terms <- BP_terms[BP_terms$NES > 0,]
posBP_terms.orderedNES <- posBP_terms[order(posBP_terms$NES, decreasing = TRUE),]

negBP_terms <- BP_terms[BP_terms$NES < 0,]
negBP_terms.orderedNES <- negBP_terms[order(negBP_terms$NES, decreasing = FALSE),]
topnegBP_terms.orderedNES <- negBP_terms.orderedNES[1:15,]

#Combine the posBP_terms.orderedNES and topnegBP_terms.orderedNES dataframes
BP <- rbind(posBP_terms.orderedNES, topnegBP_terms.orderedNES)
BP$Sign <- ifelse(BP$NES > 0, "Enriched in B6D:KO", "Enriched in A/D:KO")

#Set the GO description to a factor so it is ordrede the same way in the plot
BP$Description <- factor(BP$Description, levels = BP$Description)

#Group the results based on MF ontology and separate into positive and negative, arrange by decreasing NES
MF_terms <- gse.readable[gse.readable$ONTOLOGY == "MF", ]
posMF_terms <- MF_terms[MF_terms$NES > 0, ]
posMF_terms.orderedNES <- posMF_terms[order(posMF_terms$NES, decreasing = TRUE),]

negMF_terms <- MF_terms[MF_terms$NES < 0,]
negMF_terms.orderedNES <- negMF_terms[order(negMF_terms$NES, decreasing = FALSE),]
topnegMF_terms.orderedNES <- negMF_terms.orderedNES[1:15,]

#Combine the posBP_terms.orderedNES and topnegBP_terms.orderedNES dataframes
MF <- rbind(posMF_terms.orderedNES, topnegMF_terms.orderedNES)
MF$Sign <- ifelse(MF$NES > 0, "Enriched in B6D:KO", "Enriched in A/D:KO")

#Set the GO description to a factor so it is ordrede the same way in the plot
MF$Description <- factor(MF$Description, levels = MF$Description)

#Group the results based on CC ontology and separate into positive and negative, arrange by decreasing NES
CC_terms <- gse.readable[gse.readable$ONTOLOGY == "CC", ]
posCC_terms <- CC_terms[CC_terms$NES > 0, ]
posCC_terms.orderedNES <- posCC_terms[order(posCC_terms$NES, decreasing = TRUE),]

negCC_terms <- CC_terms[CC_terms$NES < 0,]
negCC_terms.orderedNES <- negCC_terms[order(negCC_terms$NES, decreasing = FALSE),]
topnegCC_terms.orderedNES <- negCC_terms.orderedNES[1:2,]

#Combine the posBP_terms.orderedNES and topnegBP_terms.orderedNES dataframes
CC <- rbind(posCC_terms.orderedNES, topnegCC_terms.orderedNES)
CC$Sign <- ifelse(CC$NES > 0, "Enriched in B6D:KO", "Enriched in A/D:KO")

#Set the GO description to a factor so it is ordrede the same way in the plot
CC$Description <- factor(CC$Description, levels = CC$Description)

#Rbind BP and MF
BPMF <- rbind(BP, MF)

#Rbind BPMF and CC
all <- rbind(BPMF, CC)

#Order all GO terms by decreasing NES and set description as factor'
all$Description <- str_wrap(all$Description, width = 60)
all <- all[order(all$NES, decreasing = TRUE),]
all$Description <- factor(all$Description, levels = all$Description)

# Plot
OntologyPlot <- ggplot(all, aes(x = NES, y = Description, fill = ONTOLOGY, label = setSize)) +
  geom_bar(stat = "identity") +
  geom_text(size = 7, color = "black") +
  facet_wrap(~Sign, scales = "free_x") +
  scale_x_continuous(breaks = seq(-2,2,1)) +
  theme(axis.text.y = element_text(size = 9)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 26),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 20),
        strip.text = element_text(face = "bold",size = 26),
        title = element_text(size = 26),
        legend.text = element_text(size = 26),
        legend.position = "bottom",
        legend.title = element_text(size = 0)) 
OntologyPlot
ggsave("GO_ontology_customBarplot.tiff", OntologyPlot, width = 14, height = 24, dpi = 300)




