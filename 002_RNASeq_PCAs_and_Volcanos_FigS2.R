setwd("/Users/mikemartinez/Desktop/Renamed_Samples_Figures/")


#Load libraries
suppressWarnings({
  #Data manipulation
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(plyr))
  
  
  #Differential gene expression
  suppressPackageStartupMessages(library(DESeq2))
  suppressPackageStartupMessages(library(ashr))
  suppressPackageStartupMessages(library(RNAseqQC))
  suppressPackageStartupMessages(library(vsn))
  
  #Graphics and visualizations
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggh4x))
  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(patchwork))
  
  
  #GSEA analysis
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  suppressPackageStartupMessages(library(AnnotationDbi))
  suppressPackageStartupMessages(library(msigdbr))
  suppressPackageStartupMessages(library(enrichplot))
})

raw <- read.csv("/users/mikemartinez/Desktop/Renamed_Samples_Figures/New_Name_Counts/final_counts.csv", header = TRUE, sep = ",")

#Obtain gene symbols
raw$Symbol <- mapIds(org.Mm.eg.db, key = raw$Gene, column = "SYMBOL",
                     keytype = "ENSEMBL", multiVals = "first")

#Omit any gene with no symbol annotation or unknown genes
raw <- raw[!is.na(raw$Symbol),]

#Omit any gene with no MGI annotation or unknown genes and Rikens
raw <- raw[!grepl("^Gm\\d+$", raw$Symbol),]
raw <- raw[!grepl("Rik$", raw$Symbol),]
raw <- raw[!grepl("Rik\\d+$", raw$Symbol),]
raw <- raw[!grepl("^LOC", raw$Symbol),]
raw <- raw[!grepl("AK\\d+$", raw$Symbol),]
raw <- raw[!grepl("AY\\d+$", raw$Symbol),]

#Format the rownames
raw$geneIDs <- paste(raw$Gene, raw$Symbol, sep = " - ")
rownames(raw) <- raw$geneIDs
raw$Geneid<- NULL
raw$Symbol <- NULL
raw$geneIDs <- NULL
raw$Gene <- NULL

raw$X <- NULL
#Create a design table for DESeq2
design <- data.frame(Group = rep(c("A/D:KO", "A/D:WT", "B6D:KO", "B6D:WT"),
                                 c(12,12,11,12)))
rownames(design) <- colnames(raw)

#Check that everything in the design table is in the counts matrix
all(colnames(raw) %in% rownames(design))
all(colnames(raw) == rownames(design))

suppressWarnings({
  dds <- DESeqDataSetFromMatrix(countData = raw,
                                colData = design,
                                design = ~ Group)
})

#Note: whatever you set as the reference group is the denominator according to DESeq2 documentation
dds$Group <- relevel(dds$Group, ref = "A/D:KO")

dds <- DESeq(dds)
saveRDS(dds, file = "allGroups.rds")

#Create a summarized experiment
vsd <- vst(dds)

#Plot PCA
PCA <- plotPCA(vsd, intgroup = "Group") +
  geom_text_repel(aes(label = rownames(design)), size = 0.1, max.overlaps = Inf) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  ggtitle("") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 26),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 26),
        legend.text = element_text(size = 24),
        title = element_text(size = 26)) +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  #scale_fill_manual(values = c("A/D:WT" = "green", "B6D:WT" = "purple"))
PCA
ggsave("ADKO_vs_B6DKO_PCA.tiff", PCA, width = 12, height = 8, dpi = 800)


PCA <- plotPCA(vsd, intgroup = "Group", returnData = TRUE) 
percentVar <- round(100* attr(PCA, "percentVar"))
PCA$group <- factor(PCA$group, levels = c("A/D:WT", "B6D:WT", "A/D:KO", "B6D:KO"))

plot <- ggplot(PCA, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = c("A/D:KO" = "red", "A/D:WT" = "forestgreen", "B6D:WT" = "purple", "B6D:KO" = "blue")) +
  scale_fill_manual(values = c("A/D:KO" = "red", "A/D:WT" = "forestgreen", "B6D:WT" = "purple", "B6D:KO" = "blue")) +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 26),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 26),
        legend.text = element_text(size = 24),
        title = element_text(size = 26))
ggsave("allGroups.tiff", plot, width = 12, height = 8, dpi = 800)








setwd("/users/mikemartinez/Desktop/AJB6_Ulceration/AJB6_RNASeq_Section/DESeq2_Comparison_With_IG_genes_RDS_Files/")
#Read in all the RDS files
AJKOWT <- readRDS("AJKO_vs_AJWT.rds") #A
B6KOWT <- readRDS("B6KO_vs_B6WT.rds") #B
AJB6WT <- readRDS("AJWT_vs_B6WT.rds") #C
AJB6KO <- readRDS("AJKO_vs_B6KO.rds") #D

A <- as.data.frame(results(AJKOWT))
A <- na.omit(A)
B <- as.data.frame(results(B6KOWT))
B <- na.omit(B)
C <- as.data.frame(results(AJB6WT))
C <- na.omit(C)
D <- as.data.frame(results(AJB6KO))
D <- na.omit(D)








