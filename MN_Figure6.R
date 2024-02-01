# The purpose of this script is to generate finalized timepoint plots for Masako Nakanishi

setwd("/users/michaelmartinez/Desktop/MN_Figure_6/")
library(ggplot2)
library(ggh4x)
library(ggpubr)
library(tidyverse)
library(dplyr)

# -------------------------------------------------------------------------------- #
#                               Figure 6D-J                                        #
# -------------------------------------------------------------------------------- #

# Read in the following taxa files
# Akkermansia muciniphila
amucin <- read.csv("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Species_analysis/Species_timepoint_Plots/Akkermansia_muciniphila.csv", header = TRUE, sep = ",")

# Faecalibaculum unclassified
funclass <- read.csv("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Species_analysis/Species_timepoint_Plots/Faecalibaculum_unclassified.csv", header = TRUE, sep = ",")

# Turibacter unclassified
tunclass <- read.csv("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Species_analysis/Species_timepoint_Plots/Turibacter_unclassified.csv", header = TRUE, sep = ",")

# Bacteroides vulgatus
bvulgatus <- read.csv("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Species_analysis/Species_timepoint_Plots/Bacteroides_vulgatus_counts.csv", header = TRUE, sep = ",")

# Bacteroides paurosaccharolyticus
bpauro <- read.csv("/users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Updated_uBiome_Section/uBiome_Section_Figures/Jan_Figures_Updates/Additional_timepoint_species_counts/Bacteroides_paurosaccharolyticus_counts.csv", header = TRUE, sep = ",")

# Escherichiaunclassified
Eunclass <- read.csv("/users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Updated_uBiome_Section/uBiome_Section_Figures/Jan_Figures_Updates/Additional_timepoint_species_counts/Escherichia_unclassified_counts.csv", header = TRUE, sep = ",")

# Butyrivibrio unclassified
bunclass <- read.csv("/users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Updated_uBiome_Section/uBiome_Section_Figures/Jan_Figures_Updates/Additional_timepoint_species_counts/Butyrivibrio_unclassified_counts.csv", header = TRUE, sep = ",")

# Faecalibacterium unclassified
Faecalibacterium <- read.csv("/users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Updated_uBiome_Section/uBiome_Section_Figures/Jan_Figures_Updates/Additional_timepoint_species_counts/Faecalibacterium_unclassified_counts.csv", header = TRUE, sep = ",")

# Lachnoclostridium unclassified
Lachno <- read.csv("/users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Updated_uBiome_Section/uBiome_Section_Figures/Jan_Figures_Updates/Additional_timepoint_species_counts/Lachnoclostridium_unclassified_counts.csv", header = TRUE, sep = ",")

# Prevotella unclassified
Prevotella <- read.csv("/users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Updated_uBiome_Section/uBiome_Section_Figures/Jan_Figures_Updates/Additional_timepoint_species_counts/Prevotella_unclassified_counts.csv", header = TRUE, sep = ",")

# Rufibacter unclassified
Rufi <- read.csv("/users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Updated_uBiome_Section/uBiome_Section_Figures/Jan_Figures_Updates/Additional_timepoint_species_counts/Rufibacter_unclassified_counts.csv", header = TRUE, sep = ",")



# Write a function to rename and set levels
formatData <- function(x, name) {
  Age_order <- c("8 Weeks",
                 "20 Weeks")
  geno_order <- c("WT", "KO")
  data <- x
  data$rename <- ifelse(data$Strain == "AJ", "A/D", "B6D")
  data <- data %>%
    mutate(Weeks = factor(Age, levels = Age_order))
  data <- data %>%
    mutate(Weeks = factor(Weeks, levels = Age_order))
  data <- data %>%
    mutate(Genotype, factor(Genotype, levels = geno_order))
  data$Group <- paste(data$rename, data$Genotype, sep = ":")
  data$Group <- factor(data$Group, levels = c("A/D:WT", "A/D:KO", "B6D:WT","B6D:KO"))
  data$Species <- name
  return(data)
}

# Run the function on each dataframe to return the formatted dataframe
# Akkermansia muciniphila
amucin <- formatData(amucin, "Akkermansia muciniphila")

# Faecalibaculum unclassified
funclass <- formatData(funclass, "Faecalibaculum unclassified")

# Turibacter unclassified
tunclass <- formatData(tunclass, "Turibacter unclassified")

# Bacteroides vulgatus
bvulgatus <- formatData(bvulgatus, "Bacteroides vulgatus")

# Bacteroides paurosaccharolyticus
bpauro <- formatData(bpauro, "Bacteroides paurosaccharolyticus")

# Escherichiaunclassified
Eunclass <- formatData(Eunclass, "Escherichiaunclassified")

# Butyrivibrio unclassified
bunclass <- formatData(bunclass, "Butyrivibrio unclassified")

# Faecalibacterium unclassified
Faecalibacterium <- formatData(Faecalibacterium, "Faecalibacterium unclassified")
  
# Lachnoclostridium unclassified
Lachno <- formatData(Lachno, "Lachnoclostridium unclassified")

# Prevotella unclassified
Prevotella <- formatData(Prevotella, "Prevotella unclassified")
  
# Rufibacter unclassified
Rufi <- formatData(Rufi, "Rufibacter unclassified")



# # Let's assess the normality of each 
# normality <- function(x, name) {
#   model <- aov(count ~ Group, data = x)
#   residuals <- residuals(model)
#   tiff(paste(name, "Normal Q-Q Plot.tiff", sep = "_"))
#   qqnorm(residuals, main = name)
#   qqline(residuals)
#   dev.off()
# }
# 
# # Assess each dataframe
# normality(amucin, "Akkermansia muciniphila")
# normality(funclass, "Faecalibaculum unclassified")
# normality(tunclass, "Turibacter unclassified")
# normality(bvulgatus, "Bacteroides vulgatus")
# normality(bpauro, "Bacteroides paurosaccharolyticus")
# normality(Eunclass, "Escherichia unclassified")
# normality(bunclass, "Butyrivibrio unclassified")
# normality(Faecalibacterium, "Faecalibacterium unclassified")
# normality(Lachno, "Lachnoclostridium unclassified")
# normality(Prevotella, "FPrevotella unclassified")
# normality(Rufi, "Rufibacter unclassified")
# 
# 
# # Function to run Kruskal-Wallis test
# set.seed(03061999)
# KWtest <- function(x, name) {
#   KWresult <- kruskal.test(count ~ Group, data = x)
#   print(name)
#   print(KWresult)
# }
# 
# KWtest(amucin, "Akkermansia muciniphila")
# KWtest(funclass, "Faecalibaculum unclassified")
# KWtest(tunclass, "Turibacter unclassified")
# KWtest(bvulgatus, "Bacteroides vulgatus")
# KWtest(bpauro, "Bacteroides paurosaccharolyticus")
# KWtest(Eunclass, "Escherichia unclassified")
# KWtest(bunclass, "Butyrivibrio unclassified")
# KWtest(Faecalibacterium, "Faecalibacterium unclassified")
# KWtest(Lachno, "Lachnoclostridium unclassified")
# KWtest(Prevotella, "Prevotella unclassified")
# KWtest(Rufi, "Rufibacter unclassified")


# Function to plot data
plotTimePoint <- function(x, name) {
  x$Group <- factor(x$Group, levels = c("A/D:WT", "A/D:KO", "B6D:WT", "B6D:KO"))
timepoint <- ggplot(x, aes(x = Group, y = count, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 3, position = position_jitter(width = 0.2)) +
  stat_pwc(method = "dunn_test", p.adjust.method = "BH", label = "p.adj.signif", label.size = 8, hide.ns = TRUE) +
  theme_bw() +
  scale_fill_manual(values = c("A/D:WT" = "#ffffff", "A/D:KO" = "#e0dcdc", "B6D:WT" = "#807d7d", "B6D:KO" = "#242300")) +
  facet_nested_wrap(~Species + Weeks, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  theme(axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 26),
        strip.text = element_text(size = 30, face = "bold")) +
  labs(y = "Count",
       x = "") +
  guides(fill = FALSE) 
ggsave(paste(name, "Timepoint.tiff", sep = "_"), timepoint, dpi = 300, width = 10, height = 12)
}


# Run the function on each dataframe to return the plot
# Akkermansia muciniphila
amucin.plot <- plotTimePoint(amucin, "Akkermansia muciniphila")

# Faecalibaculum unclassified
funclass.plot <- plotTimePoint(funclass, "Faecalibaculum unclassified")

# Turibacter unclassified
tunclass.plot <- plotTimePoint(tunclass, "Turibacter unclassified")

# Bacteroides vulgatus
colnames(bvulgatus)[colnames(bvulgatus) == "Count"] <- "count"
bvulgatus.plot <- plotTimePoint(bvulgatus, "Bacteroides vulgatus")

# Bacteroides paurosaccharolyticus
bpauro.plot <- plotTimePoint(bpauro, "Bacteroides paurosaccharolyticus")

# Escherichiaunclassified
Eunclass.plot <- plotTimePoint(Eunclass, "Escherichia unclassified")

# Butyrivibrio unclassified
bunclass.plot <- plotTimePoint(bunclass, "Butyrivibrio unclassified")

# Faecalibacterium unclassified
Faecalibacterium.plot <- plotTimePoint(Faecalibacterium, "Faecalibacterium unclassified")

# Lachnoclostridium unclassified
Lachno.plot <- plotTimePoint(Lachno, "Lachnoclostridium unclassified")

# Prevotella unclassified
Prevotella.plot <- plotTimePoint(Prevotella, "Prevotella unclassified")

# Rufibacter unclassified
Rufi.plot <- plotTimePoint(Rufi, "Rufibacter unclassified")


# -------------------------------------------------------------------------------- #
#                               Figure 6C                                          #
# -------------------------------------------------------------------------------- #

# Significant Family heatmap
# Requires the output from the DAA from MaAslin2 (See `DAA_Analysis.R`) and the microbiome counts at the family level
  # NOT in long format

# Read in the significant results
sigs <- read.csv("/users/michaelmartinez/Desktop/Family_DAA_Analysis/significant_results.tsv", header = TRUE, sep = "\t")
famCounts <- read.csv("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/uBiome/Family/Family_Counts.csv", header = TRUE, sep = ",")

# Extract the taxa names and extract these names from the counts file
sigsFams <- sigs$feature
hmData <- famCounts[famCounts$taxon %in% sigsFams,]

# Set the rownames of hmData to be the taxa names and remove the original column
rownames(hmData) <- hmData$taxon
hmData$taxon <- NULL

# Set the age order
Age_order <- c("8 Weeks", "20 Weeks")

# Create a heatmap annotation dataframe
Strain <- metadata$rename
names(Strain) <- metadata$Sample_ID
Strain <- as.data.frame(Strain)
Strain$Age <- metadata$Age
Strain$Genotype <- metadata$Phenotype
Strain <- Strain %>%
  mutate(Age = factor(Age, levels = Age_order))
Strain <- Strain %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "KO")))
Strain$Strain <- factor(Strain$Strain)
Strain$Age <- factor(Strain$Age)
Strain$Genotype <- factor(Strain$Genotype)
Strain$Group <- paste(Strain$Strain, Strain$Genotype, sep = ":")
Strain[,1:3] <- NULL
Strain$Strain <- Strain$Group
Strain$Group <- NULL
Strain$Age <- metadata$Age
Strain$Age <- factor(Strain$Age, levels = c("8 Weeks", "20 Weeks"))
Strain$Strain <- factor(Strain$Strain, levels = c("A/D:WT", "A/D:KO", "B6D:WT", "B6D:KO"))

# Remove unnecessary metadata columns and controls
hmData[,1:5] <- NULL
hmData[,92:95] <- NULL

# Plot the heatmap
HM <- pheatmap(hmData,
               scale = "row",
               annotation_col = Strain,
               annotation_colors = list(
                 Strain = c("A/D:WT" = "forestgreen", "A/D:KO" = "red", "B6D:WT" = "purple", "B6D:KO" = "blue"),
                 Age = c("8 Weeks" = "gray", "20 Weeks" = "black")),
               show_colnames = FALSE,
               fontsize = 20,
               color = mako(10),
               fontsize_row = 20,
               cluster_cols = FALSE,
               main = "Significant Families", size = 40,
               gaps_col = 44)
ggsave("Strain_Ordered_SignificantFamilies_Heatmap_Jan_HighRes_102323.tiff", HM, width = 18, height = 12, dpi = 300)   


# -------------------------------------------------------------------------------- #
#                               Figure 6B                                          #
# -------------------------------------------------------------------------------- #

# BF Ratio
# Requires microbiome relative abundances at the phylum level

# Read in the phylum relative abundances
phylum <- read.csv("/users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Updated_uBiome_Section/uBiome_Data/Phylum/Phylum_RelAbund.csv", header = TRUE, sep = ",")

# Lets split this into 8 weeks and 20 weeks
taxa_subset <- c("Bacteroidetes", "Firmicutes")
BFRatio <- phylum[phylum$taxon %in% taxa_subset,]

# We need to change the names of the strain
BFRatio$rename <- ifelse(BFRatio$Strain == "AD", "A/D", "B6D")

# Make a group column
BFRatio$group <- paste(BFRatio$Strain, BFRatio$Phenotype, sep = ":")

# Get a dataframe of just Firmicutes and a dataframe of just Bacteroidetes
Firm <- BFRatio[BFRatio$taxon == "Firmicutes",]
Bact <- BFRatio[BFRatio$taxon == "Bacteroidetes",]

# Calculate the F/B ratio
Firm$FBRatio <- Firm$RelAbund / Bact$RelAbund

# Write as a csv so we can manually edit the group names
write.csv(Firm, file = "FB_Ratio.csv")

# Reorder the group
group_order <- c("B6D:KO", "B6D:WT", "A/D:KO", "A/D:WT")

# Re-order
Firm <- Firm %>%
  mutate(group = factor(group, levels = group_order))

# Factor age
#Re-order
Firm <- Firm %>%
  mutate(Age = factor(Age, levels = Age_order))

# Set colors for groups to be consistent with the PCA
custom_colors <- c("B6D:KO" = "blue",
                   "B6D:WT" = "purple",
                   "A/D:KO" = "red",
                   "A/D:WT" = "forestgreen")

# Let's plot the F/B ratio as boxplots
FB_Ratio_plot <- ggplot(Firm, aes(x = FBRatio, y = group, fill = group)) +
  geom_boxplot() +
  facet_grid(~Age) +
  theme_bw() +
  labs(y = "",
       x = "F/B Ratio",
       title = "F/B Ratio") +
  theme(strip.text = element_text(size = 50, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 40),
        title = element_text(size = 30)) +
  scale_fill_manual(values = custom_colors) +
  theme(legend.position = "none")
ggsave("F_B_Ratio_boxplot.tiff", FB_Ratio_plot, width = 10, height = 10, dpi = 300)


# -------------------------------------------------------------------------------- #
#                               Figure 6A                                          #
# -------------------------------------------------------------------------------- #

#-----Phylum level relative abundance
# Requires phylum level relative abundance counts
# Read in the phylum level abundance counts
phyRel <- read.csv("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Updated_uBiome_Section/uBiome_Data/Phylum/Phylum_RelAbund.csv", header = TRUE, sep = ",")

# Function to calculate the mean frequency of phyla
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

# Calculate mean frequency using `meanFreq` function
freq <- meanFreq(phyRel)

# Convert phyla to factors with a specific order based on the mean frequency
Phyla_order <- c("Bacteroidetes",
                 "Proteobacteria",
                 "Firmicutes",
                 "Actinobacteria",
                 "Bacteria_unclassified",
                 "Verrucomicrobia",
                 "Deferribacteres",
                 "Tenericutes",
                 "Cyanobacteria",
                 "Acidobacteria",
                 "Fusobacteria",
                 "Chlorobi")
#Subset
phyRel.subset <- phyRel[phyRel$taxon %in% Phyla_order,]
phyRel.subset$Phyla <- phyRel.subset$taxon

#Re-order
phyRel <- phyRel %>%
  mutate(taxon = factor(taxon, levels = Phyla_order))
phyRel <- phyRel %>%
  mutate(Phyla = factor(Phyla, levels = Phyla_order))

#Convert Age to a factor so 8 weeks comes before 20
Age_order <- c("8 Weeks",
               "20 Weeks")
Age_order <- c("8 Wks", "20 Wks")

#Strain order
pheno_order <- c("WT", "KO")


#Re-order
phyRel.subset <- phyRel.subset %>%
  mutate(Age = factor(Age, levels = Age_order))
#Re-order
phyRel.subset <- phyRel.subset %>%
  mutate(taxon = factor(taxon, levels = Phyla_order))
#Re-order
phyRel.subset <- phyRel.subset %>%
  mutate(Phenotype = factor(Phenotype, levels = pheno_order))

# Add a slash into A/D phenotype
phyRel.subset$rename <- ifelse(phyRel.subset$Strain == "AD", "A/D", "B6D")

# Plot
phylumRelAbund <- ggplot(phyRel.subset, (aes(x = Sample_ID, y = RelAbund))) +
  geom_bar(aes(fill = Phyla), stat = "identity", position = "fill", width = 1) +
  scale_fill_brewer(palette = "Paired") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 20),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 26),
        title = element_text(size = 26),
        legend.text = element_text(size = 24)) +
  scale_y_continuous(name = "Relative Abundance",
                     labels = scales::percent) +
  scale_x_discrete(labels = NULL) +
  facet_nested_wrap(~rename + Age + Phenotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  theme(strip.background = element_rect(color = "black", fill = "lightgray"),
        panel.spacing = unit(0.2, "lines")) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
  labs(title = "Phylum-level Relative Abundance") +
  xlab(NULL)
ggsave("Renamed_PhylumLevel_RelativeAbundance_HiRes_102323.tiff", phylumRelAbund, width = 18, height = 12, dpi = 800)

