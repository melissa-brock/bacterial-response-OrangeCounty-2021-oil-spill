###Comparison of Genera Abundance from 2020 and 2021 - Melissa Brock
####Set up environment####
#Set working directory
getwd()
setwd("/Users/melissabrock/Desktop/AEM Submission/Figure_Changes/")
getwd()

#Load libraries
library(tidyverse)
library(ggpubr)
library(ggtext)
library(lubridate)

#Read in files
SCOP_Genus <- read.csv("SCOP_Genus_RelAbund.csv", header = T, row.names = 1, sep = ",")
SCOP_meta <- read.csv("NP SCOP Metadata.csv", header = T, sep = ",")

JGI_Genus <- read.csv("JGI_Genus_RelAbund.csv", header = T, row.names = 1, sep = ",")
JGI_meta <- read.csv("7b74ffeaf32d0bae8df863a1.csv", header = T, sep = ",")

####Format files - SCOP####
#Convert date to date format
SCOP_meta$Date_Combined <- paste(SCOP_meta$Correct_Month, SCOP_meta$Correct_Day, SCOP_meta$Correct_Year, sep = "/")
SCOP_meta$Correct_Date <- as.Date(SCOP_meta$Date_Combined, "%m/%d/%Y")
SCOP_meta <- SCOP_meta[c(1:27), ]

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

#Subset to genera of interest
dim(SCOP_Genus_Bacteria_spread2) #27 2495
sort(colSums(SCOP_Genus_Bacteria_spread2), decreasing = T)
SCOP_Genus_Bacteria_sorted <- as.data.frame(SCOP_Genus_Bacteria_spread2[, (order(colSums(SCOP_Genus_Bacteria_spread2), decreasing = T))])
SCOP_Genus_Bacteria_sub <- SCOP_Genus_Bacteria_sorted[, c(1, 3, 41, 27, 18, 16, 6, 4, 20)]

SCOP_Genus_Bacteria_sub$SampleID <- row.names(SCOP_Genus_Bacteria_sub)
SCOP_Genus_Bacteria_sub_gather <- gather(SCOP_Genus_Bacteria_sub, Genus_ConsensusID, Genus_Percent, -SampleID)
SCOP_Genus_Bacteria_sub_gather_meta <- merge(SCOP_Genus_Bacteria_sub_gather, SCOP_meta, by = "SampleID")

####Format files - JGI####
#Convert date to date format
JGI_meta$Date_Combined <- paste(JGI_meta$Month, JGI_meta$Day, JGI_meta$Year, sep = "/")
JGI_meta$Date <- as.Date(JGI_meta$Date_Combined, "%m/%d/%Y")
JGI_meta <- JGI_meta[c(1:267), ]
colnames(JGI_meta)[1] <- "SampleID"

#Pull out bacterial genera
unique(JGI_Genus$Domain_ConsensusID)
JGI_Genus_Bacteria <- JGI_Genus[which(JGI_Genus$Domain_ConsensusID == "Bacteria"), ]

#Merge taxa table with metadata table based on SampleID
JGI_Genus_Bacteria_meta <- merge(JGI_Genus_Bacteria, JGI_meta, by = "SampleID")

#Check number of unique genera
length(unique(JGI_Genus_Bacteria_meta$Genus_ConsensusID)) #2500

#Remove NA assignments
JGI_Genus_Bacteria_meta <- JGI_Genus_Bacteria_meta[!is.na(JGI_Genus_Bacteria_meta$Genus_ConsensusID),]

#Convert to wide format
colnames(JGI_Genus_Bacteria_meta)
JGI_Genus_Bacteria_meta2 <- JGI_Genus_Bacteria_meta[, c(1, 3, 5)]
JGI_Genus_Bacteria_spread <- spread(JGI_Genus_Bacteria_meta2, Genus_ConsensusID, Genus_Percent)

#Make second dataframe
JGI_Genus_Bacteria_spread2 <- JGI_Genus_Bacteria_spread[, -1]

#Set row names
row.names(JGI_Genus_Bacteria_spread2) <- JGI_Genus_Bacteria_spread$SampleID

#Replace NAs with 0s
JGI_Genus_Bacteria_spread2[is.na(JGI_Genus_Bacteria_spread2)] <- 0

#Subset to genera of interest
dim(JGI_Genus_Bacteria_spread2) #11 2499
sort(colSums(JGI_Genus_Bacteria_spread2), decreasing = T)
JGI_Genus_Bacteria_sorted <- as.data.frame(JGI_Genus_Bacteria_spread2[, (order(colSums(JGI_Genus_Bacteria_spread2), decreasing = T))])
JGI_Genus_Bacteria_sub <- JGI_Genus_Bacteria_sorted[, c(1, 3, 7, 8, 10, 28, 30, 44, 60)]

JGI_Genus_Bacteria_sub$SampleID <- row.names(JGI_Genus_Bacteria_sub)
JGI_Genus_Bacteria_sub_gather <- gather(JGI_Genus_Bacteria_sub, Genus_ConsensusID, Genus_Percent, -SampleID)
JGI_Genus_Bacteria_sub_gather_meta <- merge(JGI_Genus_Bacteria_sub_gather, JGI_meta, by = "SampleID")

####Scatter Plot####
colnames(JGI_Genus_Bacteria_sub_gather_meta)
colnames(SCOP_Genus_Bacteria_sub_gather_meta)

JGI_Genus_Bacteria_meta_sub <- JGI_Genus_Bacteria_sub_gather_meta[, c(1:7)]
SCOP_Genus_Bacteria_meta_sub <- SCOP_Genus_Bacteria_sub_gather_meta[, c(1:3, 6:9)]
colnames(SCOP_Genus_Bacteria_meta_sub)[4] <- "Date"
colnames(SCOP_Genus_Bacteria_meta_sub)[5] <- "Year"
colnames(SCOP_Genus_Bacteria_meta_sub)[6] <- "Month"
colnames(SCOP_Genus_Bacteria_meta_sub)[7] <- "Day"

Genus_Bacteria_meta_combined <- rbind(JGI_Genus_Bacteria_meta_sub, SCOP_Genus_Bacteria_meta_sub)
Genus_Bacteria_meta_combined$julian <- yday(Genus_Bacteria_meta_combined$Date)  

genus.compare.plot <- ggplot(Genus_Bacteria_meta_combined, aes(x = julian, y = Genus_Percent, color = as.factor(Year)))+ 
  geom_point(shape = 16, size = 3) + 
  geom_line() + 
  facet_wrap(~Genus_ConsensusID, ncol = 3, scales = "free") + 
  theme_classic() + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 20), legend.title = element_text(size = 22), panel.background = element_rect(color = "black", fill = "transparent")) +
  labs(x = "Julian Day", y = "Relative Abundance (%)", title = " ") +
  expand_limits(x = c(200, 360)) + 
  scale_x_continuous(breaks = seq(200, 360, 40)) +
  scale_color_manual( name = "Year:", values = c("seagreen4", "grey40"))

png("Genus_Comparison.png", width = 12, height = 9, units = "in", res = 600)
pdf("Genus_Comparison.pdf", width = 12, height = 9)
genus.compare.plot
dev.off()
