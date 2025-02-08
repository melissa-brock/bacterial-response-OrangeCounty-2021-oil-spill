####Calculating diversity of Taxa and changes in community structure - Melissa Brock
####Set up environment####
#Set working directory
getwd()
setwd("C:/SCOP_hpc/")
getwd()

#Load libraries
library(tidyverse)
library(ggpubr)
library(microbiome)
library(vegan)
library(phyloseq)

#Read in files
SCOP_meta <- read.csv("NP SCOP Metadata.csv", header = T, sep = ",")
SCOP_Genus_RPKM <- read.csv("SCOP_Genus_RPKM.csv", header = T, row.names = 1, sep = ",")

####Format files####
#Convert date to date format
SCOP_meta$Date_Combined <- paste(SCOP_meta$Correct_Month, SCOP_meta$Correct_Day, SCOP_meta$Correct_Year, sep = "/")
SCOP_meta$Correct_Date <- as.Date(SCOP_meta$Date_Combined, "%m/%d/%Y")

#Pull out bacterial genera
unique(SCOP_Genus_RPKM$Domain_ConsensusID)
SCOP_Genus_Bacteria_RPKM <- SCOP_Genus_RPKM[which(SCOP_Genus_RPKM$Domain_ConsensusID == "Bacteria"), ]

#Convert to wide format
colnames(SCOP_Genus_Bacteria_RPKM)
SCOP_Genus_Bacteria_RPKM2 <- SCOP_Genus_Bacteria_RPKM[, -2]
SCOP_Genus_Bacteria_RPKM_spread <- spread(SCOP_Genus_Bacteria_RPKM2, Genus_ConsensusID, Genus_SUM)

#Make second dataframe
SCOP_Genus_Bacteria_RPKM_spread2 <- SCOP_Genus_Bacteria_RPKM_spread[, -1]

#Set row names
row.names(SCOP_Genus_Bacteria_RPKM_spread2) <- SCOP_Genus_Bacteria_RPKM_spread$SampleID

#Replace NAs with 0s
SCOP_Genus_Bacteria_RPKM_spread2[is.na(SCOP_Genus_Bacteria_RPKM_spread2)] <- 0

#Transpose
SCOP_Genus_Bacteria_RPKM.t <- as.data.frame(t(SCOP_Genus_Bacteria_RPKM_spread2))

#Remove unclassified and unassigned genera
tail(row.names(SCOP_Genus_Bacteria_RPKM.t), n = 750)
dim(SCOP_Genus_Bacteria_RPKM.t)
SCOP_Genus_Bacteria_RPKM_nounclass <- SCOP_Genus_Bacteria_RPKM.t[!grepl("unclassified",row.names(SCOP_Genus_Bacteria_RPKM.t)),]
dim(SCOP_Genus_Bacteria_RPKM_nounclass)

SCOP_Genus_Bacteria_RPKM_nounass <- SCOP_Genus_Bacteria_RPKM_nounclass[!grepl("<NA>",row.names(SCOP_Genus_Bacteria_RPKM_nounclass)),]
dim(SCOP_Genus_Bacteria_RPKM_nounass)

####Calculate taxonomic richness####
SCOP_Genus_Bacteria_richness <- microbiome::richness(SCOP_Genus_Bacteria_RPKM_nounass)
SCOP_Genus_Bacteria_richness.df <- as.data.frame(SCOP_Genus_Bacteria_richness)

####Calculate taxonomic evenness####
SCOP_Genus_Bacteria_RPKM_nounass.t <- as.data.frame(t(SCOP_Genus_Bacteria_RPKM_nounass))

#Pielou’s evenness J = H′/log(S)
SCOP_Genus_Bacteria_evenness <- vegan::diversity(SCOP_Genus_Bacteria_RPKM_nounass.t, index = "shannon")/log(specnumber(SCOP_Genus_Bacteria_RPKM_nounass.t))
SCOP_Genus_Bacteria_evenness.df <- as.data.frame(SCOP_Genus_Bacteria_evenness)
View(SCOP_Genus_Bacteria_evenness.df)

####Calculate taxonomic Shannon Index####
SCOP_Genus_Bacteria_shannon <- vegan::diversity(SCOP_Genus_Bacteria_RPKM_nounass.t, index = "shannon")
SCOP_Genus_Bacteria_shannon.df <- as.data.frame(SCOP_Genus_Bacteria_shannon)

####Merge diversity metrics with metadata####
SCOP_Genus_Bacteria_richness.df$SampleID <- rownames(SCOP_Genus_Bacteria_richness.df)
colnames(SCOP_Genus_Bacteria_richness.df)
SCOP_Genus_Bacteria_richness_metadata <- merge(SCOP_meta, SCOP_Genus_Bacteria_richness.df, by = "SampleID")
range(SCOP_Genus_Bacteria_richness_metadata$observed) #1889 2306


SCOP_Genus_Bacteria_evenness.df$SampleID <- rownames(SCOP_Genus_Bacteria_evenness.df)
colnames(SCOP_Genus_Bacteria_evenness.df) <- c("Pielou_Evenness", "SampleID")
SCOP_Genus_Bacteria_evenness_metadata <- merge(SCOP_meta, SCOP_Genus_Bacteria_evenness.df, by = "SampleID")
range(SCOP_Genus_Bacteria_evenness_metadata$Pielou_Evenness) #0.4557719 0.6574370


SCOP_Genus_Bacteria_shannon.df$SampleID <- rownames(SCOP_Genus_Bacteria_shannon.df)
colnames(SCOP_Genus_Bacteria_shannon.df) <- c("Shannon", "SampleID")
SCOP_Genus_Bacteria_shannon_metadata <- merge(SCOP_meta, SCOP_Genus_Bacteria_shannon.df, by = "SampleID")
range(SCOP_Genus_Bacteria_shannon_metadata$Shannon) #3.523397 4.959575

####Taxonomic Diversity Plots####
###Richness
SCOP_Genus_Bacteria_richness_metadata$Correct_Date
SCOP_Genus_Bacteria_richness_metadata$Timeline <- c(rep("Before", 6), rep("Week 1", 3), rep("Week 2", 3), rep("Week 3", 3), rep("Week 4", 3), rep("Month 2", 4), rep("Month 3", 5), rep("Month 4", 4))

bacteria.richness.lineplot <- ggplot(SCOP_Genus_Bacteria_richness_metadata, aes(x = Correct_Date, y = observed)) + 
  theme_classic() + 
  labs(y = "Taxonomic Richness", x = "Date", title = " ") + 
  geom_smooth(method = "loess", se = F, span = 0.2, color = "black", linewidth = 1, linetype = "dashed") + 
  geom_point(aes(color = Timeline), shape = 16, size = 5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_date(breaks = SCOP_Genus_Bacteria_richness_metadata$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
bacteria.richness.lineplot


###Evenness
colnames(SCOP_Genus_Bacteria_evenness_metadata)
SCOP_Genus_Bacteria_evenness_metadata$Correct_Date
SCOP_Genus_Bacteria_evenness_metadata$Timeline <- c(rep("Before", 6), rep("Week 1", 3), rep("Week 2", 3), rep("Week 3", 3), rep("Week 4", 3), rep("Month 2", 4), rep("Month 3", 5), rep("Month 4", 4))

bacteria.evenness.lineplot <- ggplot(SCOP_Genus_Bacteria_evenness_metadata, aes(x = Correct_Date, y = Pielou_Evenness)) + 
  theme_classic() + 
  labs(y = "Taxonomic Evenness", x = "Date", title = " ") + 
  geom_smooth(method = "loess", se = F, span = 0.2, color = "black", linewidth = 1, linetype = "dashed") + 
  geom_point(aes(color = Timeline), shape = 16, size = 5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  expand_limits(y = c(0.4, 0.7)) + 
  scale_y_continuous(breaks = seq(0.4, 0.7, 0.05)) + 
  scale_x_date(breaks = SCOP_Genus_Bacteria_evenness_metadata$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
bacteria.evenness.lineplot


###Shannon
colnames(SCOP_Genus_Bacteria_shannon_metadata)
SCOP_Genus_Bacteria_shannon_metadata$Correct_Date
SCOP_Genus_Bacteria_shannon_metadata$Timeline <- c(rep("Before", 6), rep("Week 1", 3), rep("Week 2", 3), rep("Week 3", 3), rep("Week 4", 3), rep("Month 2", 4), rep("Month 3", 5), rep("Month 4", 4))

bacteria.shannon.lineplot <- ggplot(SCOP_Genus_Bacteria_shannon_metadata, aes(x = Correct_Date, y = Shannon)) + 
  theme_classic() + 
  labs(y = "Taxonomic Shannon Index", x = "Date", title = " ") + 
  geom_smooth(method = "loess", se = F, span = 0.2, color = "black", linewidth = 1, linetype = "dashed") + 
  geom_point(aes(color = Timeline), shape = 16, size = 5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  expand_limits(y = c(3, 5)) + 
  scale_y_continuous(breaks = seq(3, 5, 0.5)) + 
  scale_x_date(breaks = SCOP_Genus_Bacteria_shannon_metadata$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
bacteria.shannon.lineplot

####Correlations between Taxonomic Diversity Metrics####
SCOP_Genus_Bacteria_diversity <- merge(SCOP_Genus_Bacteria_richness.df, SCOP_Genus_Bacteria_evenness.df, by = "SampleID")
SCOP_Genus_Bacteria_diversity2 <- merge(SCOP_Genus_Bacteria_diversity, SCOP_Genus_Bacteria_shannon.df, by = "SampleID")

cor.test(y = SCOP_Genus_Bacteria_diversity2$Shannon, x = SCOP_Genus_Bacteria_diversity2$observed, method = "pearson", conf.level = 0.95) #p-value = 0.5372
cor.test(y = SCOP_Genus_Bacteria_diversity2$Shannon, x = SCOP_Genus_Bacteria_diversity2$Pielou_Evenness, method = "pearson", conf.level = 0.95) #p-value < 2.2e-16; cor = 0.9976167
plot(y = SCOP_Genus_Bacteria_diversity2$Shannon, x = SCOP_Genus_Bacteria_diversity2$Pielou_Evenness, type = "p", ylab = "Shannon Index", xlab = "Evenness", main = "Bacteria")

####PCA and PERMANOVA####
#Calculate z-score
SCOP_Genus_Bacteria_RPKM_nounass.t <- as.data.frame(t(SCOP_Genus_Bacteria_RPKM_nounass))
SCOP_Genus_Bacteria_RPKM_nounass.zscore <- as.data.frame(scale(SCOP_Genus_Bacteria_RPKM_nounass.t, center = T, scale = T))

#Make phyloseq otu table
genus.bacteria.rpkm.phy <- otu_table(SCOP_Genus_Bacteria_RPKM_nounass.zscore, taxa_are_rows = F)

#Prepare metadata sheet
row.names(SCOP_meta)<- SCOP_meta$SampleID
SCOP_meta$Timeline <- c(rep("Before", 6), rep("Week 1", 3), rep("Week 2", 3), rep("Week 3", 3), rep("Week 4", 3), rep("Month 2", 4), rep("Month 3", 5), rep("Month 4", 4))
View(SCOP_meta)

#Make phyloseq object for sample metadata
sample.metadata.phy <- sample_data(SCOP_meta)

#Make combined phyloseq object
genus.bacteria.rpkm.phyloseq <- phyloseq(genus.bacteria.rpkm.phy, sample.metadata.phy)

#Ordinate
genus.bacteria.rpkm.pca <- ordinate(genus.bacteria.rpkm.phyloseq, method = "RDA", distance = "euclidean")

#Extract eigenvalues
genus.bacteria.rpkm.pca.eigen.vals <- genus.bacteria.rpkm.pca$CA$eig

#Extract PC1 values
bacteria.loadings <- as.data.frame(genus.bacteria.rpkm.pca$CA$u)
bacteria.loadings$SampleID <- row.names(bacteria.loadings)
bacteria.loadings <- merge(SCOP_meta, bacteria.loadings, by = "SampleID")

#Find variance explained by axes
summary(genus.bacteria.rpkm.pca)
#                           PC1       PC2      
#Eigenvalue            839.5947 126.29660 
#Proportion Explained    0.3577   0.05381 
#Cumulative Proportion   0.3577   0.41154 

#Plot ordination
bacteria.pca.plot <- ggplot(data = bacteria.loadings) + 
  geom_point(aes(y = PC2, x = PC1, color = Timeline), size = 5) + 
  labs(title = "", x = "PC1 [35.8%]", y = "PC2 [5.4%]") + 
  theme_classic() + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans"), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  expand_limits(y = c(-0.4, 0.4), x = c(-0.4, 0.4)) + 
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) + 
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
bacteria.pca.plot

#Plot PC1 by time
bacteria.pc1.plot <- ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) + 
  geom_point(data = bacteria.loadings, aes(x = Correct_Date, y = PC1, color = Timeline), shape = 16, size = 5) + 
  geom_smooth(data = bacteria.loadings, aes(x = Correct_Date, y = PC1), method = "loess", se = T, color = "transparent", linewidth = 0.1) + 
  labs(title = " ", x = "Date", y = "PC1") + 
  theme_classic() + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), text = element_text(family = "sans"), ) +
  expand_limits(y = c(-0.4, 0.4)) + 
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) + 
  scale_x_date(breaks = bacteria.loadings$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
bacteria.pc1.plot

#PERMANOVA
SCOP_Genus_Bacteria_RPKM.dist <- vegdist(SCOP_Genus_Bacteria_RPKM_nounass.zscore, method = "euclidean")
adonis2(SCOP_Genus_Bacteria_RPKM.dist ~ SCOP_meta$Timeline, method = "euclidean") #p = 0.065

####Temporal beta-diversity####
#Calculate distance matrix
SCOP_Genus_Bacteria_eucdist <- as.matrix(vegdist(SCOP_Genus_Bacteria_RPKM_nounass.zscore, method = "euclidean"))

#Extract sub-diagonal of the matrix
SCOP_Genus_Bacteria_eucdist.subdiag <- as.data.frame(diag(SCOP_Genus_Bacteria_eucdist[-1,-ncol(SCOP_Genus_Bacteria_eucdist)]))
SCOP_Genus_Bacteria_eucdist.subdiag$SampleID <- row.names(SCOP_Genus_Bacteria_eucdist)[-1]
colnames(SCOP_Genus_Bacteria_eucdist.subdiag) <- c("Euc.Diss", "SampleID")

#Add row for 08/04 so this plot aligns with the other plots' date axes
SCOP_Genus_Bacteria_eucdist.subdiag[nrow(SCOP_Genus_Bacteria_eucdist.subdiag) + 1,] = c(0,"SCOP_210804")
is.numeric(SCOP_Genus_Bacteria_eucdist.subdiag$Euc.Diss)
SCOP_Genus_Bacteria_eucdist.subdiag$Euc.Diss <- as.numeric(as.character(SCOP_Genus_Bacteria_eucdist.subdiag$Euc.Diss))
SCOP_Genus_Bacteria_eucdist.subdiag[SCOP_Genus_Bacteria_eucdist.subdiag == 0] <- NA

#Combine with metadata
SCOP_Genus_Bacteria_eucdist_meta <- merge(SCOP_Genus_Bacteria_eucdist.subdiag, SCOP_meta, by = "SampleID")

#Plot
bacteria.eucdiss.plot <- ggplot() + 
  geom_line(data = SCOP_Genus_Bacteria_eucdist_meta, aes(x = Correct_Date, y = Euc.Diss), color = "black", linewidth = 1, linetype = "dashed") + 
  geom_point(data = SCOP_Genus_Bacteria_eucdist_meta, aes(x = Correct_Date, y = Euc.Diss, color = Timeline), shape = 16, size = 5) + 
  labs(title = " ", x = "Date", y = "Taxonomic Euclidean Dissimilarity") + 
  theme_classic() + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), text = element_text(family = "sans"), ) +
  expand_limits(y = c(50, 90)) + 
  scale_y_continuous(breaks = seq(50, 90, 10)) + 
  scale_x_date(breaks = SCOP_Genus_Bacteria_eucdist_meta$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
bacteria.eucdiss.plot