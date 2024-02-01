# Generate DESeq2 RDS files for downstream analysis

library(AnnotationDbi)
library(org.Mm.eg.db)
library(DESeq2)

raw <- read.csv("AJ_B6_KO_counts.csv", header = TRUE, sep = ",")

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
design <- data.frame(Group = rep(c("A/D:KO", "B6D:KO"),
                                 c(12,11)))
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
saveRDS(dds, file = "ADKO_vs_B6DKO.rds")

# Clear environment

raw <- read.csv("AJWT_AJKO_counts.tsv", header = TRUE, sep = "\t")

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
design <- data.frame(Group = rep(c("A/D:KO", "A/D:WT"),
                                 c(12,12)))
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
saveRDS(dds, file = "ADKO_vs_ADWT.rds")


# Clear environment

raw <- read.csv("AJ_B6_WT_counts.csv", header = TRUE, sep = ",")

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
design <- data.frame(Group = rep(c("A/D:WT", "B6D:WT"),
                                 c(12,12)))
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
dds$Group <- relevel(dds$Group, ref = "A/D:WT")

dds <- DESeq(dds)
saveRDS(dds, file = "ADWT_vs_B6DWT.rds")


# Clear environment

raw <- read.csv("B6WT_B6KO_counts.tsv", header = TRUE, sep = "\t")

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
design <- data.frame(Group = rep(c("B6D:WT", "B6D:KO"),
                                 c(12,11)))
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
dds$Group <- relevel(dds$Group, ref = "B6D:WT")

dds <- DESeq(dds)
saveRDS(dds, file = "B6DKO_vs_B6DWT.rds")