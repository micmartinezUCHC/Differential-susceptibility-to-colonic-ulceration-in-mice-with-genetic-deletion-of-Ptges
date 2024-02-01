# The purpose of this script is to regenerate the finalized versions of alpha diversity figures
setwd("/users/michaelmartinez/Desktop/MN_Figure_5/")
library(ggplot2)
library(ggh4x)
library(ggpubr)

# Function to plot the alpha diversity 
plotAlpha <- function(x, taxLevel) {
  
  alphaDiv <- x %>%
    group_by(Phenotype)
  
  # Re name
  alphaDiv$rename <- ifelse(alphaDiv$Strain == "AJ", "A/D", "B6D")
  alphaDiv$Group <- paste(alphaDiv$rename, phylum$Phenotype, sep = ":")
  alphaDiv$Group <- factor(alphaDiv$Group, levels = c("A/D:WT", "A/D:KO", "B6D:WT","B6D:KO"))
  alphaDiv$Level <- taxLevel
  
  Shannons_boxplot <- 
    ggplot(alphaDiv, aes(x = Group, y = Shannon, fill = rename)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = 3, position = position_jitter(width = 0.2)) +
    stat_pwc(method = "tukey_hsd", p.adjust.method = "BH", label = "p.adj.signif", label.size = 7, hide.ns = TRUE) +
    theme_bw() +
    scale_fill_manual(values = c("white", "darkgrey")) +
    facet_nested_wrap(~Level + Age, nrow = 1, scale = "free_x", 
                      strip.position = "top") +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          title = element_text(size = 26),
          strip.text = element_text(size = 24, face = "bold")) +
    labs(y = "Simpson Index",
         x = "") +
    guides(fill = FALSE) 
  ggsave(paste("Renamed", taxLevel, "Shannon_AlphaDiv.tiff", sep = "_"), Shannons_boxplot, dpi = 300, width = 12, height = 8)
  
  
  Simpsons_boxplot <- 
    ggplot(alphaDiv, aes(x = Group, y = Simpson, fill = rename)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = 3, position = position_jitter(width = 0.2)) +
    stat_pwc(method = "tukey_hsd", p.adjust.method = "BH", label = "p.adj.signif", label.size = 7, hide.ns = TRUE) +
    theme_bw() +
    scale_fill_manual(values = c("white", "darkgrey")) +
    facet_nested_wrap(~Level + Age, nrow = 1, scale = "free_x", 
                      strip.position = "top") +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          title = element_text(size = 26),
          strip.text = element_text(size = 24, face = "bold")) +
    labs(y = "Shannon Index",
         x = "") +
    guides(fill = FALSE) 
  ggsave(paste("Renamed", taxLevel, "Simp_AlphaDiv.tiff", sep = "_"), Simpsons_boxplot, dpi = 300, width = 12, height = 8)
}

Age_order <- c("8 Weeks", "20 Weeks")
phenotype_order <- c("WT", "KO")

phylum <- read.csv("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/uBiome/Phylum/Phylum_AlphaDiversity.csv", header = TRUE, sep = ",")
phylum <- phylum %>%
  mutate(Age = factor(Age, levels = Age_order))
phylum <- phylum %>%
  mutate(Phenotype = factor(Phenotype, levels = phenotype_order))
family <- read.csv("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/uBiome/Family/Family_AlphaDiversity.csv", header = TRUE, sep = ",")
family <- family %>%
  mutate(Age = factor(Age, levels = Age_order))
family <- family%>%
  mutate(Phenotype = factor(Phenotype, levels = phenotype_order))
species <- read.csv("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/uBiome/Species/Species_AlphaDiversity.csv", header = TRUE, sep = ",")
species <- species %>%
  mutate(Age = factor(Age, levels = Age_order))
species <- species %>%
  mutate(Phenotype = factor(Phenotype, levels = phenotype_order))

# Run the function
phy <- plotAlpha(phylum, "Phylum")
fam <- plotAlpha(family, "Family")
spec <- plotAlpha(species, "Species")










