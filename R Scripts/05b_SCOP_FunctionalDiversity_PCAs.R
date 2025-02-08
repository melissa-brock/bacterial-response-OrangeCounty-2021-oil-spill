###Functional analysis of SCOP short-read metagenomes processed through the NMDC-Edge Platform
####Setup environment####
getwd()
setwd("/annotation-summaries/")
getwd()

library(tidyverse)
library(phyloseq)
library(vegan)

####File import & formatting####
#Read in metadata file
SCOP.meta <- read.csv("NP SCOP Metadata.csv", header = T, sep = ",", row.names = 1)

#Format metadata file
SCOP.meta$Timeline <- c(rep("Before", 6), rep("Week 1", 3), rep("Week 2", 3), rep("Week 3", 3), rep("Week 4", 3), rep("Month 2", 4), rep("Month 3", 5), rep("Month 4", 4))

#Read in rpkm dataframes
bacteria.ko.rpkm <- read.csv("SCOP_Bacteria_KO_RPKM.csv", header = T, sep = ",", row.names = 1)

#Remove multiple assignments
bacteria.ko.sep.mult <- bacteria.ko.rpkm[grep(",", bacteria.ko.rpkm$ko),]
bacteria.ko.sep.single <- anti_join(bacteria.ko.rpkm, bacteria.ko.sep.mult, by = "ko")

#Convert to wide format
bacteria.ko.sep.single.sub <- bacteria.ko.sep.single[, -4]
bacteria.ko.sep.single.spread <- spread(bacteria.ko.sep.single.sub, ko, KO_SUM)

#Set sample id as row names
bacteria.ko.sep.single.spread.sub <- bacteria.ko.sep.single.spread[, -1]
row.names(bacteria.ko.sep.single.spread.sub) <- bacteria.ko.sep.single.spread$SampleID

#Replace NAs with 0s
bacteria.ko.sep.single.spread.sub[is.na(bacteria.ko.sep.single.spread.sub)] <- 0

#Remove unassigned functions
bacteria.ko.sep.single.spread.sub2 <- bacteria.ko.sep.single.spread.sub[, !grepl("<NA>",colnames(bacteria.ko.sep.single.spread.sub))]

#Calculate z-score for each function
bacteria.ko.sep.single.zscore <- as.data.frame(scale(bacteria.ko.sep.single.spread.sub2, center = T, scale = T))

#Distance matrix
bacteria.ko.sep.single.dist <- as.matrix(vegdist(bacteria.ko.sep.single.zscore, method = "euclidean"))

#Export files
write.csv(bacteria.ko.sep.single.dist, "Bacteria_KO_Dist.csv")

#Make phyloseq objects
SCOP.meta.phy <- sample_data(SCOP.meta)
bacteria.ko.rpkm.phy <- otu_table(bacteria.ko.sep.single.zscore, taxa_are_rows = F)

#Combine phyloseq objects
bacteria.ko.rpkm.physeq <- phyloseq(bacteria.ko.rpkm.phy, SCOP.meta.phy)

####Ordination####
bacteria.ko.rpkm.pca <- ordinate(bacteria.ko.rpkm.physeq, method = "RDA", distance = "euclidean")

#Extract eigenvalues
bacteria.ko.rpkm.pca.eigen.vals <- bacteria.ko.rpkm.pca$CA$eig

#Extract PC values
bacteria.ko.rpkm.pca.loadings <- as.data.frame(bacteria.ko.rpkm.pca$CA$u)

#Export loadings
write.csv(bacteria.ko.rpkm.pca.loadings, "Bacteria_KO_RPKM_PCA_Loadings.csv")

#Find variance explained by axes
summary(bacteria.ko.rpkm.pca)

####Close environment####
save.image("SCOP_Functional_Diversity_PCAs.RData")
quit()
