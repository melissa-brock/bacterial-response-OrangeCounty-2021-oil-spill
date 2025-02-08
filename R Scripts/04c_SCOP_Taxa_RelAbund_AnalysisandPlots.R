####Plotting Relative Abundances of Taxa and Residuals - Melissa Brock - 03/18/2024
####Set up environment####
#Set working directory
getwd()
setwd("C:/SCOP_hpc/")
getwd()

#Load libraries
library(tidyverse)
library(ggpubr)
library(ComplexHeatmap)
library(ggtext)

#Read in files
SCOP_Genus <- read.csv("SCOP_Genus_RelAbund.csv", header = T, row.names = 1, sep = ",")
SCOP_meta <- read.csv("NP SCOP Metadata.csv", header = T, sep = ",")

####Format files####
#Convert date to date format
SCOP_meta$Date_Combined <- paste(SCOP_meta$Correct_Month, SCOP_meta$Correct_Day, SCOP_meta$Correct_Year, sep = "/")
SCOP_meta$Correct_Date <- as.Date(SCOP_meta$Date_Combined, "%m/%d/%Y")

#Pull out bacterial genera
unique(SCOP_Genus$Domain_ConsensusID)
SCOP_Genus_Bacteria <- SCOP_Genus[which(SCOP_Genus$Domain_ConsensusID == "Bacteria"), ]

#Merge taxa table with metadata table based on SampleID
SCOP_Genus_Bacteria_meta <- merge(SCOP_Genus_Bacteria, SCOP_meta, by = "SampleID")

#Check number of unique genera
length(unique(SCOP_Genus_Bacteria_meta$Genus_ConsensusID)) #2496

#Remove NA assignments
SCOP_Genus_Bacteria_meta <- SCOP_Genus_Bacteria_meta[!is.na(SCOP_Genus_Bacteria_meta$Genus_ConsensusID),]

#Convert to wide format
colnames(SCOP_Genus_Bacteria_meta)
SCOP_Genus_Bacteria_meta2 <- SCOP_Genus_Bacteria_meta[, c(1, 3, 5)]
SCOP_Genus_Bacteria_spread <- spread(SCOP_Genus_Bacteria_meta2, Genus_ConsensusID, Genus_Percent)

#Make second dataframe
SCOP_Genus_Bacteria_spread2 <- SCOP_Genus_Bacteria_spread[, -1]

#Set row names
row.names(SCOP_Genus_Bacteria_spread2) <- SCOP_Genus_Bacteria_spread$SampleID

#Replace NAs with 0s
SCOP_Genus_Bacteria_spread2[is.na(SCOP_Genus_Bacteria_spread2)] <- 0

#Subset to top 25 most abundant genera and format dataframe for plotting
dim(SCOP_Genus_Bacteria_spread2) #31 2496
sort(colSums(SCOP_Genus_Bacteria_spread2), decreasing = T)
SCOP_Genus_Bacteria_sorted <- as.data.frame(SCOP_Genus_Bacteria_spread2[, (order(colSums(SCOP_Genus_Bacteria_spread2), decreasing = T))])
head(colnames(SCOP_Genus_Bacteria_sorted), n = 75)
SCOP_Genus_Bacteria_top25_nounclass <- SCOP_Genus_Bacteria_sorted[, c(1, 3, 4, 6, 9, 11, 14, 16, 18, 19, 20, 21, 23:32, 35, 36, 37)]
dim(SCOP_Genus_Bacteria_top25_nounclass) #31 25
SCOP_Genus_Bacteria_top25_nounclass$SampleID <- row.names(SCOP_Genus_Bacteria_top25_nounclass)
SCOP_Genus_Bacteria_top25_gather <- gather(SCOP_Genus_Bacteria_top25_nounclass, Genus_ConsensusID, Genus_Percent, -SampleID)
SCOP_Genus_Bacteria_top25_gather_meta <- merge(SCOP_Genus_Bacteria_top25_gather, SCOP_meta, by = "SampleID")

####Barplots####
SCOP_Genus_Bacteria_top25_barplot <- ggplot(SCOP_Genus_Bacteria_top25_gather_meta, aes(x = SampleID, y = Genus_Percent, fill = Genus_ConsensusID)) +
  geom_bar(width = 0.9, stat = "identity") + 
  theme_classic() + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 20), legend.title = element_text(size = 22), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 10), "seagreen4", "grey30", "grey30", "firebrick3", rep("grey30", 10), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  labs(x = "Date", y = "Relative Abundance (%)", title = " ") +   
    expand_limits(y = c(0, 60)) + 
    scale_y_continuous(breaks = seq(0, 60, 10)) +
  scale_x_discrete(labels = unique(SCOP_Genus_Bacteria_top25_gather_meta$Correct_Date)) +
  scale_fill_manual(guide = guide_legend(ncol = 2), name = "Genus:", values = c("black", "grey40", "grey60", "grey80", "grey90", "seashell1", "wheat3", "tan4", "pink4", "plum3", "violet", "plum1", "darkorange3", "pink1", "deeppink1", "deeppink4", "orchid4", "darkorange1", "darkseagreen1", "seashell4", "royalblue1", "royalblue4", "darkseagreen4",  "lightgoldenrod1", "goldenrod1"))
SCOP_Genus_Bacteria_top25_barplot

#Subset to taxa of interest
colnames(SCOP_Genus_Bacteria_top25_gather_meta)
unique(SCOP_Genus_Bacteria_top25_gather_meta$Genus_ConsensusID)
SCOP_Genus_Bacteria_sub <- c("Synechococcus", "Prochlorococcus", "Candidatus Pelagibacter", "Roseobacter", "Candidatus Thioglobus")

SCOP_Genus_Bacteria_top25_sub <- SCOP_Genus_Bacteria_top25_gather_meta[SCOP_Genus_Bacteria_top25_gather_meta$Genus_ConsensusID %in% SCOP_Genus_Bacteria_sub, ]

SCOP_Genus_Bacteria_top25_sub$Correct_Date <- as.Date(SCOP_Genus_Bacteria_top25_sub$Date_Combined, "%m/%d/%Y")
class(SCOP_Genus_Bacteria_top25_sub$Correct_Date)

SCOP_Genus_Bacteria_top25_sub_barplot <- ggplot(data = SCOP_Genus_Bacteria_top25_sub, mapping = aes(x = Correct_Date, y = Genus_Percent, fill = Genus_ConsensusID, group = Genus_ConsensusID)) +
  geom_bar(width = 1.5, stat = "identity", color = "black", linewidth = 0.25) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  labs(x = "Date", y = "Relative Abundance (%)", title = " ") +  
  expand_limits(y = c(0, 50)) + 
  scale_y_continuous(breaks = seq(0, 50, 10)) + 
  #scale_x_discrete(labels = unique(SCOP_Genus_Bacteria_top25_sub$Correct_Date)) +
  scale_x_date(breaks = as.Date(c("2021-08-04", "2021-08-18", "2021-09-01", "2021-09-15", "2021-09-22", "2021-09-29", "2021-10-04", "2021-10-06", "2021-10-08", "2021-10-11", "2021-10-13", "2021-10-15", "2021-10-18", "2021-10-20", "2021-10-22", "2021-10-25", "2021-10-27", "2021-10-29", "2021-11-03", "2021-11-10", "2021-11-17", "2021-11-24", "2021-12-01", "2021-12-08", "2021-12-15", "2021-12-22", "2021-12-29", "2022-01-05", "2022-01-12", "2022-01-19", "2022-01-26")), date_labels = "%Y-%m-%d") +
  scale_fill_manual(guide = guide_legend(ncol = 1), name = "Genus:", values = c("grey80", "seashell1", "darkseagreen1", "royalblue4", "darkseagreen4"))
SCOP_Genus_Bacteria_top25_sub_barplot

####DWH Enriched Taxa####
DWH.taxa <- c("Acinetobacter", "Alcanivorax", "Alteromonas", "Arcobacter", "Bacillus", "Bartonella", "Colwellia", "Cycloclasticus", "Erythrobacter", "Owenweeksia", "Polaribacter", "Halomonas", "Hyphomonas", "Labrenzia", "Marinobacter", "Marinospirillum", "Methylobacterium", "Methylocella", "Methylococcus", "Methylocystis", "Methylomonas", "Methylophaga", "Methylosinus", "Microbacterium", "Microbulbifer", "Muricauda", "Neptuniibacter", "Oceanicaulis", "Marinobacterium", "Marinomonas", "Oleibacter", "Olleya", "Pseudidiomarina", "Pseudoalteromonas", "Pseudomonas", "Rhodococcus", "Rhodovulum", "Shewanella", "Stappia", "Sulfitobacter", "Tenacibaculum", "Thalassomonas", "Thalassospira", "Vibrio")

#Pull out subset of DWH enriched taxa from dataframe
SCOP_Genus_Bacteria_DWH <- SCOP_Genus_Bacteria_spread2[, colnames(SCOP_Genus_Bacteria_spread2) %in% DWH.taxa]
SCOP_Genus_Bacteria_DWH2 <- SCOP_Genus_Bacteria_DWH
SCOP_Genus_Bacteria_DWH2$Total <- rowSums(SCOP_Genus_Bacteria_DWH2)

#Add SampleID as column 
SCOP_Genus_Bacteria_DWH$SampleID <- row.names(SCOP_Genus_Bacteria_DWH)
SCOP_Genus_Bacteria_DWH2$SampleID <- row.names(SCOP_Genus_Bacteria_DWH2)

#Convert to long format
SCOP_Genus_Bacteria_DWH_gather <- gather(SCOP_Genus_Bacteria_DWH, Genus_ConsensusID, Genus_Percent, -SampleID)

#Add Timeline as variable to metadata dataframe
SCOP_meta$Timeline <- c(rep("Before", 6), rep("Week 1", 3), rep("Week 2", 3), rep("Week 3", 3), rep("Week 4", 3), rep("Month 2", 4), rep("Month 3", 5), rep("Month 4", 4))

#Merge dataframes
SCOP_Genus_Bacteria_DWH_gather_meta <- merge(SCOP_Genus_Bacteria_DWH_gather, SCOP_meta, by = "SampleID")
SCOP_Genus_Bacteria_DWH2 <- merge(SCOP_Genus_Bacteria_DWH2, SCOP_meta, by = "SampleID")
range(SCOP_Genus_Bacteria_DWH2$Total) #3.016261 5.370284

#Line plot
SCOP_Genus_Bacteria_DWH_scatter <- ggplot() + 
  geom_line(data = SCOP_Genus_Bacteria_DWH2, aes(x = Correct_Date, y = Total), linetype = "dotted", linewidth = 1.5, color = "black") + 
  geom_point(data = SCOP_Genus_Bacteria_DWH2, aes(x = Correct_Date, y = Total, color = Timeline), shape = 16, size = 5) + 
  labs(title = " ", x = "Date", y = "Total Relative Abundance (%)") + 
  theme_classic() + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  expand_limits(y = c(2.75, 6.25)) + 
  scale_y_continuous(breaks = seq(3, 6, 0.5)) + 
  scale_x_date(breaks = SCOP_Genus_Bacteria_DWH2$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
SCOP_Genus_Bacteria_DWH_scatter

png("SCOP Total Oil-responders.png", width = 7, height = 4, units = "in", res = 600)
SCOP_Genus_Bacteria_DWH_scatter
dev.off()

#DWH scaled profiles
SCOP_Genus_Bacteria_DWH.t <- t(SCOP_Genus_Bacteria_DWH[-44])
SCOP_Genus_Bacteria_DWH.scale <- as.data.frame(t(t(SCOP_Genus_Bacteria_DWH.t)/colSums(SCOP_Genus_Bacteria_DWH.t))*100)
SCOP_Genus_Bacteria_DWH.scale2 <- as.data.frame(t(SCOP_Genus_Bacteria_DWH.scale))
SCOP_Genus_Bacteria_DWH.scale2$SampleID <- row.names(SCOP_Genus_Bacteria_DWH.scale2)
SCOP_Genus_Bacteria_DWH_scale_gather <- gather(SCOP_Genus_Bacteria_DWH.scale2, Genus_ConsensusID, Genus_Percent, -SampleID)
SCOP_Genus_Bacteria_DWH_scale_meta <- merge(SCOP_Genus_Bacteria_DWH_scale_gather, SCOP_meta, by = "SampleID")

SCOP_Genus_Bacteria_DWH_barplot_scaled <- ggplot(SCOP_Genus_Bacteria_DWH_scale_meta, aes(x = SampleID, y = Genus_Percent, fill = Genus_ConsensusID)) +
  geom_bar(width = 0.9, stat = "identity") + 
  theme_classic() + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 20), legend.title = element_text(size = 22), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 10), "seagreen4", "grey30", "grey30", "firebrick3", rep("grey30", 10), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  labs(x = "Date", y = "Scaled Relative Abundance (%)", title = " ") +  
  scale_x_discrete(labels = unique(SCOP_Genus_Bacteria_DWH_scale_meta$Correct_Date)) + 
  scale_fill_manual(guide = guide_legend(ncol = 2), name = "Genus:", values = c("black", "grey40", "grey60", "grey80", "grey90", "seashell1", "goldenrod1", "wheat3", "tan4", "darkorange1", "plum3", "violet", "plum1", "darkorange3", "coral3",  "lightsalmon1", "pink1", "deeppink1", "deeppink4", "orchid4", "mediumpurple2", "mediumpurple4", "pink4", "lightgoldenrod1", "seashell4", "royalblue1", "royalblue4", "wheat1", "lemonchiffon", "lemonchiffon3", "darkseagreen1", "palegreen3", "darkseagreen4", "lightcyan2", "lightcyan4", "aquamarine1", "aquamarine3", "turquoise3", "turquoise4", "darkslategrey", "lightcoral", "tomato", "tomato4", "lightgreen", "darkolivegreen3"))
SCOP_Genus_Bacteria_DWH_barplot_scaled

#Mean profile
View(SCOP_Genus_Bacteria_DWH_scale_mean_gather)
SCOP_Genus_Bacteria_DWH.scale.mean <- SCOP_Genus_Bacteria_DWH.scale
SCOP_Genus_Bacteria_DWH.scale.mean$Mean <- rowMeans(SCOP_Genus_Bacteria_DWH.scale.mean)

colnames(SCOP_Genus_Bacteria_DWH.scale.mean)
SCOP_Genus_Bacteria_DWH.scale.mean2 <- SCOP_Genus_Bacteria_DWH.scale.mean %>% mutate(across(SCOP_210804:SCOP_220126, ~ . - Mean))

SCOP_Genus_Bacteria_DWH.scale.mean2.t <- as.data.frame(t(SCOP_Genus_Bacteria_DWH.scale.mean2))
SCOP_Genus_Bacteria_DWH.scale.mean2.t$SampleID <- row.names(SCOP_Genus_Bacteria_DWH.scale.mean2.t)

SCOP_Genus_Bacteria_DWH_scale_mean_gather <- gather(SCOP_Genus_Bacteria_DWH.scale.mean2.t, Genus_ConsensusID, Genus_Percent, -SampleID)

SCOP_Genus_Bacteria_DWH_scale_mean_meta <- merge(SCOP_Genus_Bacteria_DWH_scale_mean_gather, SCOP_meta, by = "SampleID")

range(SCOP_Genus_Bacteria_DWH_scale_mean_meta$Genus_Percent)
SCOP_Genus_Bacteria_DWH_scaled_diff <- ggplot(SCOP_Genus_Bacteria_DWH_scale_mean_meta, aes(x = SampleID, y = Genus_Percent, color = Genus_ConsensusID)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1, linetype = "solid") + 
  geom_point(shape = 16, size = 3) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 20), legend.title = element_text(size = 22), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 10), "seagreen4", "grey30", "grey30", "firebrick3", rep("grey30", 10), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  labs(x = "Date", y = "Scaled Change in Relative Abundance from Mean (%)", title = " ")  +
  expand_limits(y = c(-6, 14)) + 
  scale_y_continuous(breaks = seq(-6, 14, 2)) +
  scale_x_discrete(labels = unique(SCOP_Genus_Bacteria_DWH_scale_mean_meta$Correct_Date)) + 
  scale_color_manual(guide = guide_legend(ncol = 2), name = "Genus:", values = c("black", "grey40", "grey60", "grey80", "grey90", "seashell1", "goldenrod1", "wheat3", "tan4", "darkorange1", "plum3", "violet", "plum1", "darkorange3", "coral3",  "lightsalmon1", "pink1", "deeppink1", "deeppink4", "orchid4", "mediumpurple2", "mediumpurple4", "pink4", "lightgoldenrod1", "seashell4", "royalblue1", "royalblue4", "wheat1", "lemonchiffon", "lemonchiffon3", "darkseagreen1", "palegreen3", "darkseagreen4", "lightcyan2", "lightcyan4", "aquamarine1", "aquamarine3", "turquoise3", "turquoise4", "darkslategrey", "lightcoral", "tomato", "tomato4", "lightgreen", "darkolivegreen3"))
SCOP_Genus_Bacteria_DWH_scaled_diff

png("SCOP Change in Oil-responders.png", width = 12, height = 6, units = "in", res = 600)
SCOP_Genus_Bacteria_DWH_scaled_diff
dev.off()

####Ridgeline plots of residuals####
#https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html
library(ggridges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggtext)
library(ggdist)
library(glue)
library(patchwork)

###Formatting 
family.resid <- read.csv("Family_resid.csv", header = T, row.names = 1, sep = ",")
View(family.resid)

#Families of interest:
#Pelagibacteraceae (SAR11)
#Synechococcaceae (Synechococcus)
#Prochloraceae (Prochlorococcus)
#Rhodobacteraceae (Roseobacter)
#Flavobacteriaceae (Polaribacter; Tenacibaculum)
#Halomonadaceae (Halomonas)
#Cryomorphaceae (Owenweeksia)
#Candidatus Thioglobus (unclassified Gammaproteobacteria)

sort(row.names(family.resid))
which(row.names(family.resid) == "Pelagibacteraceae") #2
which(row.names(family.resid) == "Synechococcaceae") #14
which(row.names(family.resid) == "Prochloraceae") #17
which(row.names(family.resid) == "Rhodobacteraceae") #3
which(row.names(family.resid) == "Flavobacteriaceae") #4
which(row.names(family.resid) == "Halomonadaceae") #43
which(row.names(family.resid) == "Cryomorphaceae") #28
which(row.names(family.resid) == "unclassifiedGammaproteobacteria") #5

family.resid.sub <- family.resid[c(2, 3, 4, 14, 17, 28, 43), ]
dim(family.resid.sub)

colnames(family.resid.sub)
family.resid.sub2 <- family.resid.sub[, c(1, 2, 15:26, 38:47, 62:73, 81:90, 103:114, 127:140, 156:165, 178:187, 204:216, 226:267)]

dim(family.resid.sub2)
family.resid.sub2.t <- as.data.frame(t(family.resid.sub2))
row.names(family.resid.sub2.t)
family.resid.sub2.t$SampleID <- row.names(family.resid.sub2.t)
family.resid.sub2.t$Year <- c(rep(2011, 12), rep(2012, 11), rep(2013, 11), rep(2014, 10), rep(2015, 13), rep(2016, 13), rep(2017, 10), rep(2018, 8), rep(2019, 15), rep(2020, 13), rep(2021, 27), rep(2022, 4))
family.resid.sub2.t$Month <- c(1, 1, 8, 8, 9, 9, 10, 10, 11, 11, 11, 12, 1, 1, 8, 8, 9, 10, 10, 10, 11, 11, 12, 1, 8, 8, 9, 10, 10, 10, 11, 11, 12, 12, 1, 1, 8, 8, 9, 10, 11, 11, 12, 12, 1, 1, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 12, 1, 8, 8, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 1, 1, 8, 8, 9, 9, 10, 10, 10, 11, 1, 1, 8, 8, 9, 9, 10, 12, 1, 1, 1, 1, 8, 8, 8, 9, 10, 10, 10, 10, 10, 11, 12, 1, 1, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 12, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 12, 1, 1, 1, 1)
family.resid.sub2.t$Day <- c(12, 26, 17, 24, 7, 21, 12, 19, 3, 16, 30, 14, 11, 25, 8, 22, 12, 3, 10, 24, 14, 28, 5, 9, 7, 14, 11, 2, 16, 30, 13, 27, 18, 24, 8, 22, 6, 20, 17, 1, 12, 26, 10, 17, 6, 14, 5, 19, 2, 16, 14, 28, 11, 25, 9, 13, 23, 20, 3, 7, 17, 31, 14, 28, 12, 26, 9, 23, 7, 21, 18, 25, 9, 23, 6, 20, 11, 18, 25, 15, 12, 17, 8, 22, 5, 19, 2, 5, 9, 16, 23, 30, 6, 7, 21, 4, 2, 9, 16, 23, 30, 13, 18, 8, 22, 12, 26, 9, 23, 7, 21, 18, 25, 2, 9, 16, 4, 18, 1, 15, 22, 29, 4, 6, 8, 11, 13, 15, 18, 20, 22, 25, 27, 29, 3, 10, 17, 24, 1, 8, 15, 22, 29, 5, 12, 19, 26)
family.resid.sub2.t[67, 11] <- 5
family.resid.sub2.t$Date_Combined <- paste(family.resid.sub2.t$Month, family.resid.sub2.t$Day, family.resid.sub2.t$Year, sep = "/")
family.resid.sub2.t$Correct_Date <- as.Date(family.resid.sub2.t$Date_Combined, "%m/%d/%Y")

family.resid.Aug <- family.resid.sub2.t[which(family.resid.sub2.t$Month == 8), ]
family.resid.Sept <- family.resid.sub2.t[which(family.resid.sub2.t$Month == 9), ]
family.resid.Oct <- family.resid.sub2.t[which(family.resid.sub2.t$Month == 10), ]
family.resid.Nov <- family.resid.sub2.t[which(family.resid.sub2.t$Month == 11), ]
family.resid.Dec <- family.resid.sub2.t[which(family.resid.sub2.t$Month == 12), ]
family.resid.Jan <- family.resid.sub2.t[which(family.resid.sub2.t$Month == 1), ]

family.resid.Aug.g <- gather(family.resid.Aug, Family_ConsensusID, Family_Resid, -SampleID, -Correct_Date, -Date_Combined, -Day, -Month, -Year)
family.resid.Sept.g <- gather(family.resid.Sept, Family_ConsensusID, Family_Resid, -SampleID, -Correct_Date, -Date_Combined, -Day, -Month, -Year)
family.resid.Oct.g <- gather(family.resid.Oct, Family_ConsensusID, Family_Resid, -SampleID, -Correct_Date, -Date_Combined, -Day, -Month, -Year)
family.resid.Nov.g <- gather(family.resid.Nov, Family_ConsensusID, Family_Resid, -SampleID, -Correct_Date, -Date_Combined, -Day, -Month, -Year)
family.resid.Dec.g <- gather(family.resid.Dec, Family_ConsensusID, Family_Resid, -SampleID, -Correct_Date, -Date_Combined, -Day, -Month, -Year)
family.resid.Jan.g <- gather(family.resid.Jan, Family_ConsensusID, Family_Resid, -SampleID, -Correct_Date, -Date_Combined, -Day, -Month, -Year)


###Plots
family.resid.Aug.g$Month_Day <- paste(family.resid.Aug.g$Month, family.resid.Aug.g$Day, sep = "/")

Aug.resid.plot <- ggplot(data = subset(family.resid.Aug.g, Year %in% c(2011:2020)), aes(y = Family_ConsensusID, x = as.numeric(Family_Resid), fill = Family_ConsensusID)) +
  geom_density_ridges(alpha = 0.6, quantile_lines = TRUE, quantiles = c(0.05, 0.95), jittered_points = TRUE, point_size = 0.75, point_alpha = 1, point_shape = 21) +
  geom_point(data = subset(family.resid.Aug.g, Year %in% c(2021)), aes(), shape = 8,  position = position_nudge(y = 0.1), size = 2)+
  geom_text(data = subset(family.resid.Aug.g, Year %in% c(2021)), aes(label = Month_Day), size = 2.5, vjust = 0.8) + 
  theme_classic() + 
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  labs(y = "Family", x = "Residuals", title = "August") + 
  scale_y_discrete(labels = c("Cryomorphaceae", "Flavobacteriaceae", "Halomonadaceae", "Pelagibacteraceae", "Prochloraceae", "Rhodobacteraceae", "Synechococcaceae")) +
  scale_fill_manual(values = c("aquamarine3", "darkslategrey", "darkorange1", "grey80", "darkseagreen1", "royalblue3", "darkseagreen4")) + 
  xlim(-3, 4) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1))
Aug.resid.plot



family.resid.Sept.g$Month_Day <- paste(family.resid.Sept.g$Month, family.resid.Sept.g$Day, sep = "/")

Sept.resid.plot <- ggplot(data = subset(family.resid.Sept.g, Year %in% c(2011:2020)), aes(y = Family_ConsensusID, x = as.numeric(Family_Resid), fill = Family_ConsensusID)) +
  geom_density_ridges(alpha = 0.6, quantile_lines = TRUE, quantiles = c(0.05, 0.95), jittered_points = TRUE, point_size = 0.75, point_alpha = 1, point_shape = 21) +
  geom_point(data = subset(family.resid.Sept.g, Year %in% c(2021)), aes(), shape = 8, position = position_nudge(y = 0.1), size = 2)+
  geom_text(data = subset(family.resid.Sept.g, Year %in% c(2021)), aes(label = Month_Day), size = 2.5, vjust = 0.8) + 
  theme_classic() + 
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  labs(y = "Family", x = "Residuals", title = "September") + 
  scale_y_discrete(labels = c("Cryomorphaceae", "Flavobacteriaceae", "Halomonadaceae", "Pelagibacteraceae", "Prochloraceae", "Rhodobacteraceae", "Synechococcaceae")) +
  scale_fill_manual(values = c("aquamarine3", "darkslategrey", "darkorange1", "grey80", "darkseagreen1", "royalblue3", "darkseagreen4"))+
  xlim(-3, 4) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1))
Sept.resid.plot



family.resid.Oct.g$Month_Day <- paste(family.resid.Oct.g$Month, family.resid.Oct.g$Day, sep = "/")
family.resid2.Oct.g$Month_Day <- paste(family.resid2.Oct.g$Month, family.resid2.Oct.g$Day, sep = "/")

Oct.resid.plot <- ggplot(data = subset(family.resid.Oct.g, Year %in% c(2011:2020)), aes(y = Family_ConsensusID, x = as.numeric(Family_Resid), fill = Family_ConsensusID)) +
  geom_density_ridges(alpha = 0.6, quantile_lines = TRUE, quantiles = c(0.05, 0.95), jittered_points = TRUE, point_size = 0.75, point_alpha = 1, point_shape = 21) +
  geom_point(data = subset(family.resid.Oct.g, Year %in% c(2021)), aes(), shape = 8,  position = position_nudge(y = 0.1), size = 2)+
  geom_text(data = subset(family.resid.Oct.g, Year %in% c(2021)), aes(label = Month_Day), size = 2.5, vjust = 0.8) + 
  theme_classic() + 
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  labs(y = "Family", x = "Residuals", title = "October") +
  scale_y_discrete(labels = c("Cryomorphaceae", "Flavobacteriaceae", "Halomonadaceae", "Pelagibacteraceae", "Prochloraceae", "Rhodobacteraceae", "Synechococcaceae")) +
  scale_fill_manual(values = c("aquamarine3", "darkslategrey", "darkorange1", "grey80", "darkseagreen1", "royalblue3", "darkseagreen4", "deeppink1")) + 
  xlim(-3, 4) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1))
Oct.resid.plot



family.resid.Nov.g$Month_Day <- paste(family.resid.Nov.g$Month, family.resid.Nov.g$Day, sep = "/")

Nov.resid.plot <- ggplot(data = subset(family.resid.Nov.g, Year %in% c(2011:2020)), aes(y = Family_ConsensusID, x = as.numeric(Family_Resid), fill = Family_ConsensusID)) +
  geom_density_ridges(alpha = 0.6, quantile_lines = TRUE, quantiles = c(0.05, 0.95), jittered_points = TRUE, point_size = 0.75, point_alpha = 1, point_shape = 21) +
  geom_point(data = subset(family.resid.Nov.g, Year %in% c(2021)), aes(), shape = 8, position = position_nudge(y = 0.1), size = 2)+
  geom_text(data = subset(family.resid.Nov.g, Year %in% c(2021)), aes(label = Month_Day), size = 2.5, vjust = 0.8) + 
  theme_classic() + 
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  labs(y = "Family", x = "Residuals", title = "November") + 
  scale_y_discrete(labels = c("Cryomorphaceae", "Flavobacteriaceae", "Halomonadaceae", "Pelagibacteraceae", "Prochloraceae", "Rhodobacteraceae", "Synechococcaceae")) +
  scale_fill_manual(values = c("aquamarine3", "darkslategrey", "darkorange1", "grey80", "darkseagreen1", "royalblue3", "darkseagreen4")) + 
  xlim(-3, 4) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1))
Nov.resid.plot



family.resid.Dec.g$Month_Day <- paste(family.resid.Dec.g$Month, family.resid.Dec.g$Day, sep = "/")

Dec.resid.plot <- ggplot(data = subset(family.resid.Dec.g, Year %in% c(2011:2020)), aes(y = Family_ConsensusID, x = as.numeric(Family_Resid), fill = Family_ConsensusID)) +
  geom_density_ridges(alpha = 0.6, quantile_lines = TRUE, quantiles = c(0.05, 0.95), jittered_points = TRUE, point_size = 0.75, point_alpha = 1, point_shape = 21) +
  geom_point(data = subset(family.resid.Dec.g, Year %in% c(2021)), aes(), shape = 8, position = position_nudge(y = 0.1), size = 2)+
  geom_text(data = subset(family.resid.Dec.g, Year %in% c(2021)), aes(label = Month_Day), size = 2.5, vjust = 0.8) + 
  theme_classic() + 
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  labs(y = "Family", x = "Residuals", title = "December") + 
  scale_y_discrete(labels = c("Cryomorphaceae", "Flavobacteriaceae", "Halomonadaceae", "Pelagibacteraceae", "Prochloraceae", "Rhodobacteraceae", "Synechococcaceae")) +
  scale_fill_manual(values = c("aquamarine3", "darkslategrey", "darkorange1", "grey80", "darkseagreen1", "royalblue3", "darkseagreen4")) + 
  xlim(-3, 4) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1))
Dec.resid.plot



family.resid.Jan.g$Month_Day <- paste(family.resid.Jan.g$Month, family.resid.Jan.g$Day, sep = "/")
family.resid.Jan.g$Month_Day

Jan.resid.plot <- ggplot(data = subset(family.resid.Jan.g, Year %in% c(2011:2021)), aes(y = Family_ConsensusID, x = as.numeric(Family_Resid), fill = Family_ConsensusID)) +
  geom_density_ridges(alpha = 0.6, quantile_lines = TRUE, quantiles = c(0.05, 0.95), jittered_points = TRUE, point_size = 0.75, point_alpha = 1, point_shape = 21) +
  geom_point(data = subset(family.resid.Jan.g, Year %in% c(2022)), aes(), shape = 8, position = position_nudge(y = 0.1), size = 2)+
  geom_text(data = subset(family.resid.Jan.g, Year %in% c(2022)), aes(label = Month_Day), size = 2.5, vjust = 0.8) + 
  theme_classic() + 
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  labs(y = "Family", x = "Residuals", title = "January") + 
  scale_y_discrete(labels = c("Cryomorphaceae", "Flavobacteriaceae", "Halomonadaceae", "Pelagibacteraceae", "Prochloraceae", "Rhodobacteraceae", "Synechococcaceae")) +
  scale_fill_manual(values = c("aquamarine3", "darkslategrey", "darkorange1", "grey80", "darkseagreen1", "royalblue3", "darkseagreen4")) + 
  xlim(-3, 4) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1))
Jan.resid.plot

####Combine figures####
setwd("~/SCOP/Figures/Taxonomic Patterns and Diversity/")

#residuals
png("SCOP_Genera_Residuals_Density2.png", units = "in", width = 15, height = 8, res = 600)
ggarrange(Aug.resid.plot, Sept.resid.plot, Oct.resid.plot, Nov.resid.plot, Dec.resid.plot, Jan.resid.plot, nrow = 2, ncol = 3, align = "hv", labels = c("a", "b", "c", "d", "e", "f"), font.label = list(size = 16, color = "black", face = "bold", family = "sans"))
dev.off()

svg("SCOP_Genera_Residuals_Density2.svg", width = 15, height = 8)
ggarrange(Aug.resid.plot, Sept.resid.plot, Oct.resid.plot, Nov.resid.plot, Dec.resid.plot, Jan.resid.plot, nrow = 2, ncol = 3, align = "hv", labels = c("a", "b", "c", "d", "e", "f"), font.label = list(size = 16, color = "black", face = "bold", family = "sans"))
dev.off()

#bact div and comp
png("SCOP_BactDiv_Comp2.png", units = "in", width = 21, height = 13, res = 600)
ggarrange(bacteria.evenness.lineplot, bacteria.pc1.plot, SCOP_Genus_Bacteria_top25_sub_barplot, SCOP_Genus_Bacteria_DWH_scatter, ncol = 2, nrow = 2, align = "hv", labels = c("a", "b", "c", "d"), font.label = list(size = 24, color = "black", face = "bold", family = "sans"))
dev.off()

png("SCOP_Genera_Barplots_Supp.png", units = "in", width = 20, height = 24, res = 600)
ggarrange(SCOP_Genus_Bacteria_top25_barplot, SCOP_Genus_Bacteria_DWH_barplot_scaled, SCOP_Genus_Bacteria_DWH_scaled_diff,  ncol = 1, nrow = 3, align = "hv", labels = c("a", "b", "c"), font.label = list(size = 24, color = "black", face = "bold", family = "sans"))
dev.off()