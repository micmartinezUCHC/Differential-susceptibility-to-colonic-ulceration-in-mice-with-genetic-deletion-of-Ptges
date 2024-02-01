# The purpose of this script is to do the differential abundance analysis for Masako Nakanishi


speciesCounts <- read.csv("/users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Updated_uBiome_Section/uBiome_Data/Species/Species_Counts.csv",
                          header = TRUE, sep = ",")
speciesCounts[,1:5] <- NULL
speciesCounts[,93:96] <- NULL


duplicated_rows <- duplicated(speciesCounts$taxon)

# Use logical indexing to remove duplicated rows
speciesCounts <- speciesCounts[!duplicated_rows, ]

rownames(speciesCounts) <- speciesCounts$taxon
species <- rownames(speciesCounts)
speciesCounts$taxon <- NULL
sampleIDs <- colnames(speciesCounts)

specMaas <- as.data.frame(t(speciesCounts))

# Now we get some metadata
metadata <- read.csv("/users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Raw_Data/Metadata.csv", header = TRUE, sep = ",")

metadata$rename <- ifelse(metadata$Strain == "AJ", "A/D", "B6D")
metadata$Group <- paste(metadata$rename, metadata$Phenotype, sep = ":")
rownames(metadata) <- metadata$Sample_ID

# Now we can run maaslin
fit_data <- Maaslin2(input_data = specMaas,
                     input_metadata = metadata,
                     min_prevalence = 0,
                     normalization = "NONE",
                     output = "Species_DAA_Analysis",
                     fixed_effects = "rename",
                     reference = "A/D",
                     plot_heatmap = TRUE,
                     plot_scatter = TRUE,
                     heatmap_first_n = 50)

# Now we need to all the species to their famalies
maaslinRes <- read.csv("/users/michaelmartinez/Desktop/Family_DAA_Analysis/significant_results.tsv", header = TRUE, sep = "\t")
shoreline <- read.csv("/users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Raw_Data/all_shoreline_data_copy.csv", header = TRUE, sep = ",")

specNames <- maaslinRes$feature
specLines <- shoreline[shoreline$taxon %in% specNames,]
species <- specLines$taxon
names(species) <- specLines$rankID
species <- as.data.frame(species)
species$RankID <- rownames(species)
species$rankID <- sub("^(\\d+\\.\\d+\\.\\d+\\.\\d+\\.\\d+\\.\\d+).*", "\\1", species$RankID)

familyRankIDs <- species$rankID
species$RankID <- NULL

families <- shoreline[shoreline$taxlevel == 5, ]
families$X <- NULL
family <- families[,c(2,5)]
familyNames <- species$rankID
familySubset <- family[family$rankID %in% familyNames,]

familyNames <- unique(familySubset$rankID)
names(familyNames) <- familySubset$taxon

merged <- merge(family, species, by = "rankID", all.y = FALSE)
uniquemerged <- !duplicated(merged)
merged <- subset(merged, uniquemerged)

colnames(merged) <- c("rankID", "Family", "feature")

finalRes <- merge(maaslinRes, merged, by = "feature")
finalRes <- finalRes[,c(1,11,2:ncol(finalRes))]
finalRes$Family.1 <- NULL
finalRes$DAA <- ifelse(finalRes$coef > 0, "Enriched in B6D", "Enriched in A/D")

write.csv(finalRes, file = "Significant_Species_DAA_analysis.csv")
saveRDS(fit_data, file = "MN_Maaslin_Species.Rds")


# For family
familyCounts <- read.csv("/users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/Updated_uBiome_Section/uBiome_Data/Family/Family_Counts.csv",
                          header = TRUE, sep = ",")
familyCounts[,1:5] <- NULL
familyCounts[,93:96] <- NULL


duplicated_rows <- duplicated(familyCounts$taxon)

# Use logical indexing to remove duplicated rows
familyCounts <- familyCounts[!duplicated_rows, ]

rownames(familyCounts) <- familyCounts$taxon
family <- rownames(familyCounts)
familyCounts$taxon <- NULL
sampleIDs <- colnames(familyCounts)

familyMaas <- as.data.frame(t(familyCounts))

fam_fit_data <- Maaslin2(input_data = familyMaas,
                     input_metadata = metadata,
                     min_prevalence = 0,
                     normalization = "NONE",
                     output = "Family_DAA_Analysis",
                     fixed_effects = "rename",
                     reference = "A/D",
                     plot_heatmap = TRUE,
                     plot_scatter = TRUE,
                     heatmap_first_n = 50)

# Read in the significant results
sigs <- read.csv("/users/michaelmartinez/Desktop/Family_DAA_Analysis/significant_results.tsv", header = TRUE, sep = "\t")
famCounts <- read.csv("/Users/michaelmartinez/Desktop/AJB6_Ulceration/AJB6_uBiome_Section/uBiome/Family/Family_Counts.csv", header = TRUE, sep = ",")

sigsFams <- sigs$feature
hmData <- famCounts[famCounts$taxon %in% sigsFams,]

rownames(hmData) <- hmData$taxon
hmData$taxon <- NULL

Age_order <- c("8 Weeks", "20 Weeks")
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

hmData[,1:5] <- NULL
hmData[,92:95] <- NULL

# order the metadata
metadata$ORDER <- paste(metadata$Group, metadata$Age, sep = "_")
metadata$ORDER <- factor(metadata$Group, levels = c("A/D:WT_8 Weeks", "A/D:KO_8_Weeks", 
                                                    "B6D:WT_8 Weeks", "B6D:KO_8 Weeks",
                                                    "A/D:WT_20 Weeks", "A/D:KO_20_Weeks", 
                                                    "B6D:WT_20 Weeks", "B6D:KO_20 Weeks")) 
desired_order <- c("A/D:WT_8 Weeks", "A/D:KO_8_Weeks", 
                   "B6D:WT_8 Weeks", "B6D:KO_8 Weeks",
                   "A/D:WT_20 Weeks", "A/D:KO_20_Weeks", 
                   "B6D:WT_20 Weeks", "B6D:KO_20 Weeks")

sampleOrder <- rownames(metadata)
hmData <- hmData[,sampleOrder]


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






