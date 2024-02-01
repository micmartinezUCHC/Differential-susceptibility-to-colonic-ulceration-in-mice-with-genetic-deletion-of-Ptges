#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Microbiome analysis pipeline: A/J B6 mPGES-1 KO analysis
#  Mike Martinez
#  Rosenberg Lab, University of Connecticut Health Center
#  September 28th, 2023
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----------Load libraries
library(dplyr)
library(tidyverse)
library(broom)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(roxygen2)

#----------Set working directory
setwd("/Users/michaelmartinez/Desktop/test/")
parentDir <- file.path("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section")

#-----------Read in the microbiome data
raw <- read.csv("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Raw_Data/all_shoreline_data_copy.csv",
                header = TRUE, sep = ',')
raw$X <- NULL

#----------Read in the metadata
meta <- read.csv("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Raw_Data/Metadata.csv",
                 header = TRUE, sep = ",")

#-----------Get a list of all the unique tax levels (exclusing the first one because we don't want the root)
tax_levels <- unique(raw$taxlevel)[-1]
names(tax_levels) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
tax_levels <- as.data.frame(tax_levels)
tax_levels$taxonomy <- rownames(tax_levels)
taxonomic_levels <- tax_levels$taxonomy

#----------Initialize empty lists
taxonomy_dfs <- list() #Holds the raw counts for each taxonomic level
long_dfs <- list() #Holds the counts for each taxonomic level pivoted in long format
relative_abundance <- list() #Holds the relative abundances for each taxonomic level in long format

#---------------------------------------------#
###############################################
#####----------Declare functions----------#####
###############################################
#---------------------------------------------#
#'Function to calculate relative abundances
#'@counts a dataframe where column 1 is taxa names at a given taxonomic level followed by raw counts for all samples
#'@relabund return value is a data frame containing the relative abundances of each taxa provided in the input df
#Function to calculate relative abundance
relAbund <- function(counts) {
  sums <- rowSums(counts[,2:ncol(counts)]) + 0.01
  relabund <- counts[,2:ncol(counts)]/sums
  return(relabund)
}

#'Function to calculate top 12 most frequent taxa
#'@relabund a relative abundance data frame that has been pivoted to long format
#'@freq a data frame where column 1 is the taxa name and column 2 is the mean frequency of that taxa
meanFreq <- function(relabund) {
  freq <- relabund %>%
    group_by(Sample_ID, taxon) %>%
    summarise(Count = sum(RelAbund),
              .groups = 'drop') %>%
    group_by(Sample_ID) %>%
    summarise(Freq = Count / sum(Count),
              Taxa = taxon,
              .groups = 'drop') %>%
    group_by(Taxa) %>%
    summarise(mean = mean(Freq),
              .groups = 'drop') %>%
    arrange(desc(mean))
  
  return(as.data.frame(freq))
}

#'Function for relative abundance barplots
#'@x a data frame of relative abundances only for the most frequent taxa
#'@tax_level the taxonomic level you are plotting
#'@taxa_order the desired order of the taxa you are plotting (decreasing mean frequency, as a factor)
barplot <- function(x, tax_level, taxa_order) {
  
  Age_order <- c("8 Weeks", "20 Weeks")
  phenotype_order <- c("WT", "KO")
  
  x <- x %>%
    mutate(Age = factor(Age, levels = Age_order))
  x <- x %>%
    mutate(Phenotype = factor(Phenotype, levels = phenotype_order))
  x <- x %>%
    mutate(taxon = factor(taxon, levels = taxa_order))
  
  barplot <- ggplot(x, (aes(x = Sample_ID, y = RelAbund))) +
    geom_bar(aes(fill = taxon), stat = "identity", position = "fill", width = 1) +
    scale_fill_brewer(palette = "Paired") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          strip.text = element_text(face = "bold", size = 12)) +
    facet_nested_wrap(~Strain + Age + Phenotype, nrow = 1, scale = "free_x", 
                      strip.position = "top") +
    scale_y_continuous(name = "Relative Abundance",
                       labels = scales::percent) +
    theme(strip.background = element_rect(color = "black", fill = "lightgray"),
          panel.spacing = unit(0.2, "lines")) +
    theme(legend.text = element_text(size = 12)) +
    theme(plot.title = element_text(size = 16)) +
    theme(text = element_text(family = "Helvetica")) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
    labs(x = NULL,
         title = paste(tax_levels$taxonomy[long_df], "Level Relative Abundances", sep = " "))
  ggsave(paste(tax_level, "RelAbund_Barplot.pdf", sep = "_"), barplot, width = 12, height = 8)
}

#' #'Function to test for significance at the taxa level
#' #'@x a relative abundance data frame
#' #'@significance return value of significant taxa at a BH-corrected pvalue of 0.05
#' significance <- function(x){
#'   significant <- x %>%
#'     nest(data = -taxon) %>%
#'     mutate(test = map(.x = data, ~kruskal.test(RelAbund ~ Strain + Phenotype, data = .x) %>%
#'                         tidy)) %>%
#'     unnest(test) %>%
#'     mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
#'     filter(p.adj < 0.05)
#'   return(significant)
#' }



#'Function to calculate alpha diversity metrics
#'@x is a dataframe of raw counts in the long format
richness <- function(x){
  sum(x > 0)
}

#'Function to calculate Shannon alpha diversity metric
#'@x is a dataframe of raw counts in the long format
shannon <- function(x){
  rabund <- x[x>0]/sum(x)
  -sum(rabund * log(rabund))
}

#'Function to calculate Simpson alpha diversity metric
#'@x is a dataframe of raw counts in the long format
simpson <- function(x){
  n <- sum(x)
  sum(x * (x-1) / (n * (n-1)))
}

plotAlpha <- function(x, tax_level) {
  
  alphaDiv <- x %>%
    group_by(Phenotype)
  
  Age_order <- c("8 Weeks", "20 Weeks")
  phenotype_order <- c("WT", "KO")
  
  x <- x %>%
    mutate(Age = factor(Age, levels = Age_order))
  x <- x %>%
    mutate(Phenotype = factor(Phenotype, levels = phenotype_order))
  
  Shannons_boxplot <- ggplot(alphaDiv, aes(x = Strain, y = Shannon, fill = Strain)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 2) +
    geom_point(size = 0.3, position = "jitter") +
    facet_nested_wrap(~Age + Phenotype, nrow = 1, scale = "free_x", 
                      strip.position = "top") +
    stat_compare_means(paired = FALSE, label = "p.format") +
    labs(title = paste(tax_levels$taxonomy[df], "Shannon Diversity", sep = " ")) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(paste(tax_level, "Shannon.pdf", sep = "_"), Shannons_boxplot, width = 12, height = 8)
  
  
  Simpsons_boxplot <- ggplot(alphaDiv, aes(x = Strain, y = Simpson, fill = Strain)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 2) +
    geom_point(size = 0.3, position = "jitter") +
    facet_nested_wrap(~Age + Phenotype, nrow = 1, scale = "free_x", 
                      strip.position = "top") +
    stat_compare_means(paired = FALSE, label = "p.format") +
    labs(title = paste(tax_levels$taxonomy[df], "Simpson Diversity", sep = " ")) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(paste(tax_level, "Simpson.pdf", sep = "_"), Simpsons_boxplot, width = 12, height = 8)
  
}



#---------------------------------------------#
###############################################
#----------       Code Blocks       ----------#
###############################################
#---------------------------------------------#

#----------Create a vector of output directories for each taxonomic level
outDirs <- list()
cat("#----- GENERATING OUTPUT DIRECTORIES -----#")
for (tax_level in seq_along(taxonomic_levels)){
  output_dir <- file.path("/Users/michaelmartinez/Desktop/test", taxonomic_levels[tax_level])
  outDirs[[tax_level]] <- output_dir
}
cat("Completed!")

#----------For each taxonomic level, subset the raw counts to only include the i-th taxonomic level
for (i in 1:nrow(tax_levels)) {
  #Create directory for each taxonomic level
  dir.create(outDirs[[i]])
  cat(paste(outDirs[[i]], "created successfully", sep = " "))
  cat("\n")
  temp <- raw[raw$taxlevel == i,]
  taxonomy_dfs[[i]] <- temp
  write.csv(temp, file = paste(outDirs[[i]], paste(tax_levels$taxonomy[i], "Counts.csv", sep = "_"), sep = "/"))
  
  #Reset working directory to parent directory
  setwd(parentDir)
}

#----------Iterate over the taxonomy_dfs list and pivot longer, merge metadata, and store to list
cat("#----- ISOLATING COUNTS AT EACH TAXONOMIC LEVEL, CONVERTING COUNTS TO RELATIVE ABUNDANCES, ALPHA DIVERSITY CALCS -----#")
for(df in 1:length(taxonomy_dfs)) {
  
  #Set working directory
  setwd(outDirs[[df]])
  
  #Get just the taxa names and counts columns
  counts <- taxonomy_dfs[[df]][,5:96] #Only the counts columns
  names <- counts$taxon #Taxon names to append as a new column in the abundances df
  
  #Run the relative abundance function and pivot longer
  cat(paste("Calculating relative abundance at the ", tax_levels$taxonomy[df], "level", sep = " "))
  cat("\n")
  abundances <- relAbund(counts) #Run function
  abundances$taxon <- names #Append names to new column
  cat("-----------Done! \n")
  
  #Pivot longer and merge metadata
  abund.long <- abundances %>%
    pivot_longer(-taxon, names_to = "Sample_ID", values_to = "RelAbund")
  abund.long.merged <- merge(abund.long, meta, by = "Sample_ID", all.y = TRUE)
  relative_abundance[[df]] <- abund.long.merged
  
  #Pivot longer the data frame, only keep taxonomy column and non-control sample columns
  long <- taxonomy_dfs[[df]][,5:96] %>%
    pivot_longer(-taxon, names_to = "Sample_ID", values_to = "count") 
  long.merged <- merge(long, meta, by = "Sample_ID", all.y = TRUE)
  long.merged$group <- paste(long.merged$Strain, long.merged$Phenotype, sep = ":")
  long_dfs[[df]] <- long.merged
  write.csv(long.merged, file = paste(outDirs[[df]], paste(tax_levels$taxonomy[df], "long.meta.csv", sep = "_"), sep = "/"))
  
  
  #Calculate alpha diversity metrics using the alpha diversity functions
  cat(paste("Calculating alpha diversity at the", tax_levels$taxonomy[df], "level", sep = " "))
  cat("\n")
  alpha <- long.merged %>%
    group_by(Sample_ID) %>%
    summarize(sobs = richness(count),
              Shannon = shannon(count),
              Simpson = simpson(count))
  
  #Merge alpha diversity metrics to counts and metadata and write to a csv file
  merged.alpha <- merge(meta, alpha, by = "Sample_ID", all.y = TRUE)
  write.csv(merged.alpha, file = paste(outDirs[[df]], paste(tax_levels$taxonomy[df], "AlphaDiversity.csv", sep = "_"), sep = "/"))
  
  #Plot alpha diversity plots using function
  cat(paste("Plotting Shannon and Simpson diversity at the", tax_levels$taxonomy[df], "level", sep = " "))
  cat("\n")
  plotAlpha(merged.alpha, tax_levels$taxonomy[df])
  cat("-----------Done! \n")
  
  #Reset parent directory
  setwd(parentDir)
  
}
cat("Completed!")

cat("#----- APPENDING METADATA AND WRITING DATA TO OUTPUT DIRECTORIES -----#")
#----------Iterate over the list including metadata-including long dataframes and write as a csv
for (long_df in 1:length(long_dfs)) {
  
  #Set working directory for specified taxonomic level and write the Long.Meta data to a csv file
  setwd(outDirs[[long_df]])
  write.csv(long_dfs[[long_df]], file = paste(outDirs[[long_df]], paste(tax_levels$taxonomy[long_df], "Long.Meta.csv", sep = "_"), sep = "/"))
  
  #Reset parent directory
  setwd(parentDir)
}
cat("Completed!")








































