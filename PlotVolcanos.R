# The followig script is for making volcano plots for Masako's A/D B6D paper
setwd("/users/michaelmartinez/Desktop/Working_Directory/")

# Specify path to Rds files
dds.files <- list.files("/users/michaelmartinez/Desktop/AJB6_Ulceration/Renamed_Samples_Figures/Renamed_RDS_files/")
path <- "/users/michaelmartinez/Desktop/AJB6_Ulceration/Renamed_Samples_Figures/Renamed_RDS_files/"

# Iterate through the Rds objects and read them in
for (i in dds.files) {
  dds <- readRDS(paste(path, i, sep = ""))
  res <- as.data.frame(results(dds))
  res <- na.omit(res)
  res.ordered <- res[order(res$log2FoldChange, decreasing = TRUE),]
  labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(res.ordered))
  res.ordered$Symbols <- labels
  
  # Get counts and append to results
  counts <- as.data.frame(counts(dds, normalized = TRUE))
  counts$Symbols <- gsub("^[^-]+-(.*)$", "\\1", rownames(counts))
  res.ordered <- merge(res.ordered, counts, by = "Symbols", all = TRUE)
  res.ordered <- na.omit(res.ordered)
  
  volcano <- EnhancedVolcano(res.ordered,
                  lab = res.ordered$Symbols,
                  x = "log2FoldChange",
                  y = "padj",
                  pCutoff = 0.5,
                  labSize = 7,
                  drawConnectors = FALSE,
                  widthConnectors = 0.75,
                  boxedLabels = FALSE,
                  max.overlaps = 59,
                  title = "", subtitle = "",
                  legendPosition = "bottom")
  volcano <- volcano +
    theme(axis.text.x = element_text(size = 24),
          axis.text.y = element_text(size = 24),
          axis.title.x = element_text(size = 26),
          axis.title.y = element_text(size = 26),
          legend.text = element_text(size = 24)) +
    theme(aspect.ratio = 1) 
  
  ggsave(paste(i, "Labelled_Volcano.tiff", sep = "_"), volcano, width = 12, height = 10, dpi = 300)
}



