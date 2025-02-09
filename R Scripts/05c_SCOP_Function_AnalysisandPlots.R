###SCOP Functional Diversity
####Set up environment####
#Set working directory
getwd()
setwd("C:/SCOP_hpc/")
getwd()

#Load libraries
library(tidyverse)
library(ggpubr)
library(vegan)
library(phyloseq)
library(RVAideMemoire)
library(ComplexHeatmap)
library(ggtext)

#Read in files
bact.ko.shannon <- read.csv("SCOP_Bacteria_KO_Shannon.csv", header = T, sep = ",")
bact.ko.richness <- read.csv("SCOP_Bacteria_KO_Richness.csv", header = T, sep = ",")
bact.ko.evenness <- read.csv("SCOP_Bacteria_KO_Evenness.csv", header = T, sep = ",")

SCOP.meta <- read.csv("NP SCOP Metadata.csv", header = T, sep = ",")

####Format files####
#Rename columns
colnames(bact.ko.shannon) <- c("SampleID", "Bacterial_Shannon_KO")
colnames(bact.ko.richness) <- c("SampleID", "Bacterial_Richness_KO")
colnames(bact.ko.evenness) <- c("SampleID", "Bacterial_Evenness_KO")

#Convert date to date format
SCOP.meta$Date_Combined <- paste(SCOP.meta$Correct_Month, SCOP.meta$Correct_Day, SCOP.meta$Correct_Year, sep = "/")
SCOP.meta$Correct_Date <- as.Date(SCOP.meta$Date_Combined, "%m/%d/%Y")

#Merge files
bact.ko.m <- merge(bact.ko.shannon, bact.ko.richness, by = "SampleID")
bact.ko.m2 <- merge(bact.ko.m, bact.ko.evenness, by = "SampleID")
bact.ko.m3 <- merge(bact.ko.m2, SCOP.meta, by = "SampleID")

####Diversity Plots - Bacteria####
range(bact.ko.m3$Bacterial_Shannon_KO) #7.282141 7.427068
range(bact.ko.m3$Bacterial_Richness_KO) #3154 4278
range(bact.ko.m3$Bacterial_Evenness_KO) #0.8813561 0.9089281

colnames(bact.ko.m3)
bact.ko.m3$Correct_Date
bact.ko.m3$Timeline <- c(rep("Before", 6), rep("Week 1", 3), rep("Week 2", 3), rep("Week 3", 3), rep("Week 4", 3), rep("Month 2", 4), rep("Month 3", 5), rep("Month 4", 4))
bact.ko.m3$Timeline2 <- c(rep("Before", 6), rep("M1", 12), rep("M2", 4), rep("M3", 5), rep("M4", 4))

range(bact.ko.m3$Bacterial_Shannon_KO)
bacteria.ko.shannon.plot <- ggplot(bact.ko.m3, aes(x = Correct_Date, y = Bacterial_Shannon_KO)) + 
  theme_classic() + 
  labs(y = "Functional Shannon Index", x = "Date", title = " ") + 
  geom_smooth(method = "loess", se = F, span = 0.3, color = "black", linewidth = 1, linetype = "dashed") + 
  geom_point(aes(color = Timeline), shape = 16, size = 5) + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  scale_x_date(breaks = bact.ko.m3$Correct_Date, date_labels = "%Y-%m-%d") + 
  expand_limits(y = c(7.25, 7.45)) + 
  scale_y_continuous(breaks = seq(7.25, 7.45, 0.05)) + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
bacteria.ko.shannon.plot

range(bact.ko.m3$Bacterial_Richness_KO)
bacteria.ko.richness.plot <- ggplot(bact.ko.m3, aes(x = Correct_Date, y = Bacterial_Richness_KO)) + 
  theme_classic() + 
  labs(y = "Functional Richness", x = "Date", title = " ") + 
  geom_smooth(method = "loess", se = F, span = 0.2, color = "black", linewidth = 1, linetype = "dashed") + 
  geom_point(aes(color = Timeline), shape = 16, size = 5) + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  scale_x_date(breaks = bact.ko.m3$Correct_Date, date_labels = "%Y-%m-%d") + 
  expand_limits(y = c(3000, 4500)) + 
  scale_y_continuous(breaks = seq(3000, 4500, 250)) + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
bacteria.ko.richness.plot

range(bact.ko.m3$Bacterial_Evenness_KO)
bacteria.ko.evenness.plot <- ggplot(bact.ko.m3, aes(x = Correct_Date, y = Bacterial_Evenness_KO)) + 
  theme_classic() + 
  labs(y = "Functional Evenness", x = "Date", title = " ") + 
  geom_smooth(method = "loess", se = F, span = 0.2, color = "black", linewidth = 1, linetype = "dashed") + 
  geom_point(aes(color = Timeline), shape = 16, size = 5) + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  scale_x_date(breaks = bact.ko.m3$Correct_Date, date_labels = "%Y-%m-%d") + 
  expand_limits(y = c(0.87, 0.92)) + 
  scale_y_continuous(breaks = seq(0.87, 0.92, 0.01)) + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
bacteria.ko.evenness.plot

####Correlation of Alpha-diversity metrics####
cor.test(y = bact.ko.m3$Bacterial_Shannon_KO, x = bact.ko.m3$Bacterial_Richness_KO, method = "pearson", conf.level = 0.95) #p-value = 0.0004425; cor = 0.5926958
cor.test(y = bact.ko.m3$Bacterial_Shannon_KO, x = bact.ko.m3$Bacterial_Evenness_KO, method = "pearson", conf.level = 0.95) #p-value = 0.666; cor = -0.08070512

####PCAs####
SCOP.meta$Correct_Date
SCOP.meta$Timeline <- c(rep("Before", 6), rep("Week 1", 3), rep("Week 2", 3), rep("Week 3", 3), rep("Week 4", 3), rep("Month 2", 4), rep("Month 3", 5), rep("Month 4", 4))

bacteria.ko.loadings <- read.csv(file = "Bacteria_KO_RPKM_PCA_Loadings.csv", header = T, row.names = 1, sep= ",")
bacteria.ko.loadings$SampleID <- row.names(bacteria.ko.loadings)
bacteria.ko.loadings <- merge(SCOP.meta, bacteria.ko.loadings, by = "SampleID")

bacteria.ko.pca <- ggplot(data = bacteria.ko.loadings) + 
  geom_point(aes(y = PC2, x = PC1, color = Timeline), size = 5) + 
  labs(title = " ", x = "PC1 [33.4%]", y = "PC2 [11.7%]") + 
  theme_classic() + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24)) +
  expand_limits(y = c(-0.4, 0.4), x = c(-0.4, 0.4)) + 
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) + 
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
bacteria.ko.pca

bacteria.ko.pc1 <- ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  geom_smooth(data = bacteria.ko.loadings, aes(x = Correct_Date, y = PC1), method = "loess", se = T, color = "transparent", size = 0.1) + 
  geom_point(data = bacteria.ko.loadings, aes(x = Correct_Date, y = PC1, color = Timeline), shape = 16, size = 5) + 
  labs(title = " ", x = "Date", y = "PC1") + 
  theme_classic() + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  expand_limits(y = c(-0.4, 0.4)) + 
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) + 
  scale_x_date(breaks = bacteria.ko.loadings$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
bacteria.ko.pc1


###PERMANOVAs
bacteria.ko.dist <- read.csv(file = "Bacteria_KO_Dist.csv", header = T, row.names = 1, sep = ",")
bacteria.ko.dist2 <- as.dist(bacteria.ko.dist)
adonis2(bacteria.ko.dist2 ~ SCOP.meta$Timeline, method = "euclidean") #p = 0.013
bacteria.ko.perm.manova <- pairwise.perm.manova(bacteria.ko.dist2, SCOP.meta$Timeline, nperm = 999, p.method = "hochberg")
bacteria.ko.perm.manova

####Temporal beta-diversity####
#Extract sub-diagonal of the matrix
bacteria.ko.dist.m <- as.matrix(bacteria.ko.dist)
bacteria.ko.dist.subdiag <- as.data.frame(diag(bacteria.ko.dist.m[-1,-ncol(bacteria.ko.dist.m)]))
bacteria.ko.dist.subdiag$SampleID <- row.names(bacteria.ko.dist)[-1]
colnames(bacteria.ko.dist.subdiag) <- c("Euc.Diss", "SampleID")

#Add row for 08/04 so the date axis aligns with the other plots
View(bacteria.ko.dist.subdiag)
bacteria.ko.dist.subdiag[nrow(bacteria.ko.dist.subdiag) + 1,] = c(0,"SCOP_210804")
is.numeric(bacteria.ko.dist.subdiag$Euc.Diss)
bacteria.ko.dist.subdiag$Euc.Diss <- as.numeric(as.character(bacteria.ko.dist.subdiag$Euc.Diss))
bacteria.ko.dist.subdiag[bacteria.ko.dist.subdiag == 0] <- NA

#Combine with metadata
bacteria.ko.dist.meta <- merge(bacteria.ko.dist.subdiag, SCOP.meta, by = "SampleID")

#Plot
bacteria.ko.eucdiss.plot <- ggplot() + 
  geom_line(data = bacteria.ko.dist.meta, aes(x = Correct_Date, y = Euc.Diss), color = "black", linewidth = 1, linetype = "dashed") + 
  geom_point(data = bacteria.ko.dist.meta, aes(x = Correct_Date, y = Euc.Diss, color = Timeline), shape = 16, size = 5) + 
  labs(title = " ", x = "Date", y = "Functional Euclidean Dissimilarity") + 
  theme_classic() + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  expand_limits(y = c(50, 130)) + 
  scale_y_continuous(breaks = seq(50, 130, 10)) + 
  scale_x_date(breaks = bacteria.ko.dist.meta$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
bacteria.ko.eucdiss.plot

####Indicator species analysis####
library(indicspecies)

row.names(bact.ko.rpkm2.spread2.t)
SCOP.meta$Timeline

SCOP_KO_Bacteria_indval <- multipatt(bact.ko.rpkm2.spread2.t, SCOP.meta$Timeline, control = how(nperm=999)) 
summary(SCOP_KO_Bacteria_indval)

#Multilevel pattern analysis
  
#  Association function: IndVal.g
#Significance level (alpha): 0.05

#Total number of species: 5564
#Selected number of species: 389 
#Number of species associated to 1 group: 88 
#Number of species associated to 2 groups: 50 
#Number of species associated to 3 groups: 38 
#Number of species associated to 4 groups: 46 
#Number of species associated to 5 groups: 46 
#Number of species associated to 6 groups: 40 
#Number of species associated to 7 groups: 81 

#List of species associated to each combination: 
  
#  Group Month 2  #sps.  56 
#stat p.value   
#KO:K02828 0.880   0.003 **
#  KO:K12686 0.866   0.002 **
#  KO:K02760 0.866   0.008 **
#  KO:K08969 0.856   0.007 **
#  KO:K08152 0.853   0.002 **
#  KO:K03488 0.830   0.010 **
#  KO:K01222 0.815   0.007 **
#  KO:K02118 0.815   0.010 **
#  KO:K16963 0.813   0.011 * 
#  KO:K18967 0.808   0.021 * 
#  KO:K18209 0.807   0.012 * 
#  KO:K11931 0.805   0.007 **
#  KO:K05823 0.804   0.010 **
#  KO:K19336 0.802   0.006 **
#  KO:K01236 0.798   0.003 **
#  KO:K07290 0.797   0.008 **
#  KO:K09388 0.792   0.012 * 
#  KO:K02123 0.787   0.012 * 
#  KO:K10022 0.782   0.014 * 
#  KO:K07692 0.778   0.035 * 
#  KO:K02826 0.776   0.015 * 
#  KO:K02466 0.774   0.010 **
#  KO:K02402 0.772   0.008 **
#  KO:K03893 0.771   0.032 * 
#  KO:K08168 0.771   0.019 * 
#  KO:K03805 0.770   0.009 **
#  KO:K11936 0.770   0.020 * 
#  KO:K16511 0.768   0.031 * 
#  KO:K05593 0.768   0.013 * 
#  KO:K02829 0.763   0.015 * 
#  KO:K02859 0.760   0.023 * 
#  KO:K06714 0.760   0.030 * 
#  KO:K11614 0.760   0.020 * 
#  KO:K02757 0.759   0.023 * 
#  KO:K13963 0.759   0.011 * 
#  KO:K02759 0.758   0.012 * 
#  KO:K09990 0.755   0.019 * 
#  KO:K02761 0.753   0.046 * 
#  KO:K03897 0.748   0.047 * 
#  KO:K07699 0.746   0.018 * 
#  KO:K03491 0.743   0.037 * 
#  KO:K05518 0.743   0.019 * 
#  KO:K16838 0.742   0.020 * 
#  KO:K11144 0.740   0.033 * 
#  KO:K02351 0.733   0.042 * 
#  KO:K04767 0.727   0.047 * 
#  KO:K06375 0.727   0.022 * 
#  KO:K05311 0.726   0.039 * 
#  KO:K06285 0.726   0.021 * 
#  KO:K11440 0.726   0.036 * 
#  KO:K16044 0.724   0.032 * 
#  KO:K18770 0.723   0.036 * 
#  KO:K03708 0.723   0.045 * 
#  KO:K04775 0.719   0.032 * 
#  KO:K00078 0.707   0.043 * 
#  KO:K18118 0.704   0.048 * 
  
#  Group Month 3  #sps.  4 
#stat p.value   
#KO:K19345 0.775   0.002 **
#KO:K03382 0.733   0.045 * 
#KO:K17105 0.725   0.042 * 
#KO:K10979 0.692   0.045 * 
  
#  Group Month 4  #sps.  6 
#stat p.value   
#KO:K08722 0.866   0.003 **
#KO:K03224 0.707   0.041 * 
#KO:K06632 0.707   0.049 * 
#KO:K07459 0.707   0.041 * 
#KO:K12240 0.707   0.043 * 
#KO:K18434 0.707   0.041 * 
  
#  Group Week 1  #sps.  2 
#stat p.value  
#KO:K08222 0.816   0.033 *
#KO:K01621 0.702   0.036 *
  
#  Group Week 2  #sps.  6 
#stat p.value   
#KO:K10643 0.816   0.026 * 
#KO:K11165 0.816   0.024 * 
#KO:K13573 0.812   0.006 **
#KO:K18892 0.777   0.009 **
#KO:K06869 0.757   0.038 * 
#KO:K17235 0.698   0.048 * 
  
#  Group Week 3  #sps.  2 
#stat p.value  
#KO:K05792 0.816   0.031 *
#KO:K19425 0.759   0.036 *
  
#  Group Week 4  #sps.  12 
#stat p.value   
#KO:K01067 0.868   0.009 **
#KO:K07085 0.849   0.006 **
#KO:K11535 0.822   0.019 * 
#KO:K01446 0.816   0.032 * 
#KO:K01727 0.816   0.032 * 
#KO:K13932 0.816   0.028 * 
#KO:K15533 0.816   0.028 * 
#KO:K16210 0.816   0.030 * 
#KO:K18552 0.816   0.030 * 
#KO:K16444 0.759   0.024 * 
#KO:K03293 0.753   0.030 * 
#KO:K16248 0.734   0.036 * 
  
#  Group Before+Month 3  #sps.  1 
#stat p.value  
#KO:K12976 0.746    0.03 *
  
#  Group Before+Week 4  #sps.  1 
#stat p.value  
#KO:K10857 0.726    0.04 *
  
#  Group Month 2+Month 3  #sps.  14 
#stat p.value   
#KO:K06044 0.848   0.003 **
#KO:K07028 0.838   0.014 * 
#KO:K09684 0.825   0.022 * 
#KO:K03670 0.814   0.010 **
#KO:K08224 0.814   0.028 * 
#KO:K13530 0.804   0.014 * 
#KO:K07639 0.798   0.011 * 
#KO:K11177 0.790   0.025 * 
#KO:K07650 0.788   0.019 * 
#KO:K10038 0.773   0.010 **
#KO:K05875 0.748   0.038 * 
#KO:K11741 0.745   0.031 * 
#KO:K18369 0.745   0.023 * 
#KO:K19290 0.744   0.034 * 
  
#  Group Month 2+Month 4  #sps.  3 
#stat p.value   
#KO:K07777 0.859   0.005 **
#KO:K05916 0.834   0.009 **
#KO:K03201 0.707   0.032 * 
  
#  Group Month 2+Week 1  #sps.  1 
#stat p.value  
#KO:K17762 0.778   0.043 *
  
#  Group Month 2+Week 2  #sps.  1 
#stat p.value  
#KO:K14633 0.769   0.018 *
  
#  Group Month 2+Week 3  #sps.  1 
#stat p.value   
#KO:K02236 0.828   0.006 **
  
#  Group Month 2+Week 4  #sps.  11 
#stat p.value   
#KO:K11907 0.960   0.003 **
#KO:K15342 0.880   0.002 **
#KO:K07248 0.873   0.003 **
#KO:K18697 0.832   0.009 **
#KO:K13955 0.794   0.018 * 
#KO:K16961 0.790   0.018 * 
#KO:K00929 0.777   0.019 * 
#KO:K09155 0.763   0.024 * 
#KO:K01777 0.758   0.035 * 
#KO:K09685 0.756   0.046 * 
#KO:K06439 0.708   0.048 * 
  
#  Group Month 3+Month 4  #sps.  8 
#stat p.value   
#KO:K07054 0.882   0.002 **
#KO:K13652 0.855   0.008 **
#KO:K12511 0.807   0.049 * 
#KO:K02282 0.805   0.016 * 
#KO:K02279 0.787   0.031 * 
#KO:K03304 0.745   0.032 * 
#KO:K09981 0.745   0.028 * 
#KO:K15066 0.725   0.047 * 
  
#  Group Month 4+Week 2  #sps.  4 
#stat p.value  
#KO:K13408 0.778   0.013 *
#KO:K07543 0.756   0.012 *
#KO:K07077 0.748   0.033 *
#KO:K09980 0.696   0.050 *
  
#Group Week 1+Week 2  #sps.  2 
#stat p.value   
#KO:K11537 0.880   0.002 **
#KO:K12263 0.843   0.003 **
  
#  Group Week 1+Week 3  #sps.  1 
#stat p.value  
#KO:K06758 0.809   0.046 *
  
#  Group Week 1+Week 4  #sps.  1 
#stat p.value  
#KO:K15554 0.782   0.028 *
  
#  Group Week 2+Week 3  #sps.  1 
#stat p.value   
#KO:K14647 0.831   0.002 **
  
#  Group Before+Month 3+Week 1  #sps.  3 
#stat p.value   
#KO:K06444 0.869   0.017 * 
#  KO:K10256 0.860   0.009 **
#  KO:K17239 0.756   0.035 * 
  
#  Group Before+Month 4+Week 1  #sps.  1 
#stat p.value  
#KO:K16384 0.784   0.019 *
  
#  Group Month 2+Month 3+Month 4  #sps.  12 
#stat p.value   
#KO:K02827 0.861   0.004 **
#KO:K05559 0.859   0.012 * 
#KO:K05560 0.832   0.009 **
#KO:K03449 0.823   0.013 * 
#KO:K07645 0.798   0.021 * 
#KO:K07644 0.789   0.027 * 
#KO:K05596 0.784   0.025 * 
#KO:K07499 0.784   0.011 * 
#KO:K06957 0.771   0.019 * 
#KO:K05563 0.766   0.049 * 
#KO:K00638 0.752   0.048 * 
#KO:K09891 0.734   0.048 * 
  
#  Group Month 2+Month 3+Week 1  #sps.  1 
#stat p.value  
#KO:K18069 0.742    0.04 *
  
#  Group Month 2+Month 3+Week 2  #sps.  1 
#stat p.value  
#KO:K03319 0.799   0.041 *
  
#  Group Month 2+Month 3+Week 4  #sps.  2 
#stat p.value   
#KO:K02770 0.847   0.017 * 
#KO:K03837 0.838   0.002 **
  
#  Group Month 2+Month 4+Week 2  #sps.  1 
#stat p.value  
#KO:K03166 0.794   0.028 *
  
#  Group Month 2+Month 4+Week 4  #sps.  2 
#stat p.value  
#KO:K00466 0.772   0.049 *
#  KO:K03260 0.748   0.050 *
  
#  Group Month 2+Week 1+Week 4  #sps.  4 
#stat p.value   
#KO:K02073 0.840   0.008 **
#  KO:K04041 0.820   0.015 * 
#  KO:K07816 0.758   0.030 * 
#  KO:K11693 0.751   0.030 * 
  
#  Group Month 2+Week 2+Week 3  #sps.  1 
#stat p.value   
#KO:K00153 0.856   0.003 **
  
#  Group Month 2+Week 2+Week 4  #sps.  2 
#stat p.value  
#KO:K07651 0.854   0.022 *
#  KO:K04065 0.741   0.044 *
  
#  Group Month 3+Month 4+Week 2  #sps.  1 
#stat p.value  
#KO:K01146 0.764   0.019 *
  
#  Group Month 3+Month 4+Week 4  #sps.  1 
#stat p.value  
#KO:K07218 0.764   0.034 *
  
#  Group Month 3+Week 1+Week 2  #sps.  1 
#stat p.value  
#KO:K18615 0.782   0.028 *
  
#  Group Month 3+Week 2+Week 3  #sps.  1 
#stat p.value  
#KO:K12278 0.794   0.024 *
  
#  Group Month 3+Week 2+Week 4  #sps.  1 
#stat p.value  
#KO:K06330 0.739   0.039 *
  
#  Group Month 4+Week 1+Week 3  #sps.  1 
#stat p.value  
#KO:K15255 0.728   0.042 *
  
#  Group Week 1+Week 2+Week 3  #sps.  1 
#stat p.value   
#KO:K15303 0.859   0.002 **
  
#  Group Week 2+Week 3+Week 4  #sps.  1 
#stat p.value  
#KO:K06986 0.771   0.026 *
  
#  Group Before+Month 2+Month 4+Week 4  #sps.  1 
#stat p.value  
#KO:K02736 0.896   0.024 *
  
#  Group Before+Month 3+Week 1+Week 2  #sps.  1 
#stat p.value  
#KO:K08918 0.868    0.02 *
  
#  Group Before+Week 1+Week 2+Week 3  #sps.  1 
#stat p.value   
#KO:K16135 0.855   0.003 **
  
#  Group Month 2+Month 3+Month 4+Week 1  #sps.  1 
#stat p.value  
#KO:K03721 0.75   0.037 *
  
#  Group Month 2+Month 3+Month 4+Week 2  #sps.  7 
#stat p.value   
#KO:K07014 0.903   0.003 **
#  KO:K00697 0.884   0.003 **
#  KO:K07278 0.853   0.015 * 
#  KO:K06968 0.840   0.006 **
#  KO:K11622 0.832   0.012 * 
#  KO:K03669 0.823   0.013 * 
#  KO:K01058 0.787   0.050 * 
  
#  Group Month 2+Month 3+Month 4+Week 3  #sps.  1 
#stat p.value  
#KO:K11618 0.791   0.019 *
  
#  Group Month 2+Month 3+Month 4+Week 4  #sps.  5 
#stat p.value   
#KO:K02855 0.931   0.034 * 
#  KO:K02283 0.924   0.006 **
#  KO:K07347 0.903   0.007 **
#  KO:K05561 0.848   0.039 * 
#  KO:K16012 0.814   0.020 * 
  
#  Group Month 2+Month 3+Week 1+Week 2  #sps.  4 
#stat p.value  
#KO:K12369 0.854   0.022 *
#  KO:K12372 0.850   0.023 *
#  KO:K07720 0.821   0.048 *
#  KO:K00086 0.775   0.023 *
  
#  Group Month 2+Month 3+Week 1+Week 4  #sps.  1 
#stat p.value  
#KO:K03697 0.83   0.041 *
  
#  Group Month 2+Month 3+Week 2+Week 3  #sps.  1 
#stat p.value    
#KO:K02446 0.961   0.001 ***
  
#  Group Month 2+Month 4+Week 1+Week 4  #sps.  1 
#stat p.value  
#KO:K17752 0.782   0.044 *
  
#  Group Month 2+Month 4+Week 2+Week 3  #sps.  1 
#stat p.value  
#KO:K13488 0.735   0.031 *
  
#  Group Month 2+Month 4+Week 2+Week 4  #sps.  5 
#stat p.value    
#KO:K02298 0.903   0.001 ***
#  KO:K19147 0.800   0.043 *  
#  KO:K16648 0.790   0.019 *  
#  KO:K02318 0.763   0.033 *  
#  KO:K10712 0.735   0.044 *  
  
#  Group Month 2+Month 4+Week 3+Week 4  #sps.  1 
#stat p.value  
#KO:K06412 0.852   0.011 *
  
#  Group Month 2+Week 1+Week 2+Week 3  #sps.  3 
#stat p.value    
#KO:K07740 0.924   0.001 ***
#  KO:K16516 0.841   0.023 *  
#  KO:K18136 0.734   0.045 *  
  
#  Group Month 2+Week 1+Week 2+Week 4  #sps.  3 
#stat p.value    
#KO:K04098 0.890   0.001 ***
#  KO:K10549 0.872   0.004 ** 
#  KO:K18843 0.779   0.022 *  
  
#  Group Month 2+Week 1+Week 3+Week 4  #sps.  2 
#stat p.value  
#KO:K07652 0.777   0.042 *
#  KO:K15640 0.734   0.048 *
  
#  Group Month 2+Week 2+Week 3+Week 4  #sps.  2 
#stat p.value    
#KO:K01387 0.921   0.001 ***
#  KO:K16556 0.750   0.034 *  
  
#  Group Month 3+Month 4+Week 2+Week 4  #sps.  1 
#stat p.value    
#KO:K18200 0.941   0.001 ***
  
# Group Month 3+Week 1+Week 2+Week 3  #sps.  1 
#stat p.value  
#KO:K09726 0.843   0.039 *
  
#  Group Month 4+Week 2+Week 3+Week 4  #sps.  2 
#stat p.value  
#KO:K06141 0.877   0.021 *
#  KO:K07966 0.832   0.030 *
  
#  Group Week 1+Week 2+Week 3+Week 4  #sps.  1 
#stat p.value  
#KO:K13283 0.862   0.021 *
  
#  Group Before+Month 2+Month 3+Week 2+Week 3  #sps.  1 
#stat p.value   
#KO:K15059 0.915   0.005 **
  
#  Group Before+Month 3+Month 4+Week 1+Week 2  #sps.  1 
#stat p.value  
#KO:K05795 0.86   0.042 *
  
#  Group Before+Week 1+Week 2+Week 3+Week 4  #sps.  1 
#stat p.value  
#KO:K02100 0.908   0.024 *
  
#  Group Month 2+Month 3+Month 4+Week 1+Week 2  #sps.  1 
#stat p.value  
#KO:K05713 0.881   0.039 *
  
#  Group Month 2+Month 3+Month 4+Week 1+Week 4  #sps.  1 
#stat p.value  
#KO:K16013 0.853   0.017 *
  
#  Group Month 2+Month 3+Month 4+Week 2+Week 3  #sps.  7 
#stat p.value   
#KO:K16698 0.950   0.005 **
#  KO:K16203 0.912   0.018 * 
#  KO:K02481 0.881   0.015 * 
#  KO:K03976 0.877   0.002 **
#  KO:K17883 0.844   0.027 * 
#  KO:K12952 0.836   0.036 * 
#  KO:K13571 0.829   0.044 * 
  
#  Group Month 2+Month 3+Month 4+Week 2+Week 4  #sps.  8 
#stat p.value   
#KO:K09940 0.916   0.006 **
#  KO:K06591 0.913   0.004 **
#  KO:K07283 0.883   0.022 * 
#  KO:K07346 0.861   0.027 * 
#  KO:K05794 0.853   0.017 * 
#  KO:K03737 0.848   0.028 * 
#  KO:K18299 0.837   0.044 * 
#  KO:K00801 0.806   0.026 * 
  
#  Group Month 2+Month 3+Month 4+Week 3+Week 4  #sps.  2 
#stat p.value  
#KO:K16291 0.883   0.015 *
#  KO:K02660 0.877   0.012 *
  
#  Group Month 2+Month 3+Week 1+Week 2+Week 3  #sps.  1 
#stat p.value  
#KO:K00260 0.838   0.042 *
  
#  Group Month 2+Month 3+Week 1+Week 2+Week 4  #sps.  5 
#stat p.value  
#KO:K07491 0.916   0.014 *
#  KO:K10210 0.895   0.046 *
#  KO:K03830 0.891   0.021 *
#  KO:K19430 0.857   0.021 *
#  KO:K19428 0.835   0.019 *
  
#  Group Month 2+Month 3+Week 1+Week 3+Week 4  #sps.  1 
#stat p.value  
#KO:K03781 0.885   0.036 *
  
#  Group Month 2+Month 3+Week 2+Week 3+Week 4  #sps.  6 
#stat p.value   
#KO:K07750 0.878   0.014 * 
#  KO:K13713 0.842   0.047 * 
#  KO:K06937 0.839   0.007 **
#  KO:K16937 0.826   0.041 * 
#  KO:K10004 0.823   0.024 * 
#  KO:K07653 0.782   0.040 * 
  
#  Group Month 2+Month 4+Week 2+Week 3+Week 4  #sps.  3 
#stat p.value   
#KO:K00621 0.916   0.002 **
#  KO:K16317 0.881   0.012 * 
#  KO:K15125 0.862   0.019 * 
  
#  Group Month 2+Week 1+Week 2+Week 3+Week 4  #sps.  6 
#stat p.value   
#KO:K18455 0.949   0.002 **
#  KO:K00371 0.930   0.014 * 
#  KO:K04074 0.927   0.003 **
#  KO:K02299 0.884   0.033 * 
#  KO:K00055 0.882   0.010 **
#  KO:K05304 0.849   0.027 * 
  
#  Group Month 3+Month 4+Week 1+Week 2+Week 3  #sps.  1 
#stat p.value  
#KO:K09932 0.806   0.038 *
  
#  Group Month 3+Month 4+Week 2+Week 3+Week 4  #sps.  1 
#stat p.value  
#KO:K02336 0.901   0.033 *
  
#  Group Before+Month 2+Month 3+Month 4+Week 1+Week 2  #sps.  1 
#stat p.value  
#KO:K07492 0.867   0.049 *
  
#  Group Before+Month 2+Month 3+Week 1+Week 2+Week 3  #sps.  1 
#stat p.value  
#KO:K02713 0.951   0.033 *
  
#  Group Before+Month 2+Month 3+Week 1+Week 2+Week 4  #sps.  1 
#stat p.value   
#KO:K06163 0.927   0.006 **
  
#  Group Before+Month 2+Month 4+Week 1+Week 2+Week 3  #sps.  1 
#stat p.value   
#KO:K18649 0.914   0.008 **
  
# Group Before+Month 2+Week 1+Week 2+Week 3+Week 4  #sps.  7 
#stat p.value    
#KO:K19058 0.970   0.011 *  
#  KO:K09952 0.964   0.001 ***
#  KO:K14213 0.959   0.002 ** 
#  KO:K10216 0.901   0.031 *  
#  KO:K18824 0.879   0.005 ** 
#  KO:K08195 0.877   0.026 *  
#  KO:K00313 0.849   0.031 *  
  
#  Group Before+Month 3+Month 4+Week 1+Week 2+Week 3  #sps.  1 
#stat p.value  
#KO:K18441 0.887   0.019 *
  
#  Group Month 2+Month 3+Month 4+Week 1+Week 2+Week 3  #sps.  4 
#stat p.value  
#KO:K03547 0.929   0.011 *
#  KO:K14664 0.897   0.043 *
#  KO:K00362 0.892   0.038 *
#  KO:K07678 0.847   0.037 *
  
#  Group Month 2+Month 3+Month 4+Week 1+Week 2+Week 4  #sps.  5 
#stat p.value   
#KO:K03311 0.944   0.007 **
#  KO:K13771 0.897   0.017 * 
#  KO:K11475 0.883   0.034 * 
#  KO:K17754 0.875   0.025 * 
#  KO:K06962 0.798   0.047 * 
  
#  Group Month 2+Month 3+Month 4+Week 2+Week 3+Week 4  #sps.  16 
#stat p.value    
#KO:K11782 0.987   0.001 ***
#  KO:K11783 0.954   0.007 ** 
#  KO:K07544 0.939   0.004 ** 
#  KO:K00299 0.932   0.012 *  
#  KO:K07029 0.926   0.022 *  
#  KO:K15538 0.919   0.006 ** 
#  KO:K12977 0.918   0.005 ** 
#  KO:K02413 0.896   0.006 ** 
#  KO:K01728 0.893   0.004 ** 
#  KO:K07181 0.892   0.007 ** 
#  KO:K15723 0.879   0.006 ** 
#  KO:K09607 0.870   0.014 *  
#  KO:K05802 0.862   0.028 *  
#  KO:K00484 0.860   0.033 *  
#  KO:K13009 0.847   0.021 *  
#  KO:K16046 0.826   0.032 *  
  
#  Group Month 2+Month 3+Week 1+Week 2+Week 3+Week 4  #sps.  3 
#stat p.value   
#KO:K19200 0.933   0.006 **
#  KO:K00884 0.908   0.007 **
#  KO:K07321 0.849   0.023 * 
  
# Group Before+Month 2+Month 3+Month 4+Week 1+Week 2+Week 3  #sps.  20 
#stat p.value   
#KO:K02290 0.989   0.018 * 
#  KO:K05574 0.989   0.037 * 
#  KO:K05579 0.988   0.038 * 
#  KO:K15226 0.988   0.032 * 
#  KO:K09835 0.987   0.044 * 
#  KO:K07769 0.987   0.029 * 
#  KO:K02697 0.986   0.024 * 
#  KO:K07079 0.986   0.038 * 
#  KO:K00555 0.986   0.028 * 
#  KO:K05572 0.986   0.047 * 
#  KO:K01664 0.986   0.008 **
#  KO:K02700 0.985   0.035 * 
#  KO:K02096 0.984   0.029 * 
#  KO:K09942 0.984   0.039 * 
#  KO:K02097 0.984   0.023 * 
#  KO:K11145 0.976   0.042 * 
#  KO:K02094 0.972   0.038 * 
#  KO:K02724 0.972   0.035 * 
#  KO:K05585 0.971   0.042 * 
#  KO:K02289 0.968   0.042 * 
  
#  Group Before+Month 2+Month 3+Week 1+Week 2+Week 3+Week 4  #sps.  6 
#stat p.value  
#KO:K19353 0.973   0.036 *
#  KO:K01478 0.970   0.011 *
#  KO:K10118 0.963   0.046 *
#  KO:K03812 0.923   0.028 *
#  KO:K05784 0.923   0.011 *
#  KO:K06166 0.903   0.030 *
  
#  Group Before+Month 2+Month 4+Week 1+Week 2+Week 3+Week 4  #sps.  5 
#stat p.value   
#KO:K03818 0.987   0.003 **
#  KO:K18365 0.983   0.021 * 
#  KO:K02304 0.974   0.007 **
#  KO:K12567 0.949   0.032 * 
#  KO:K02668 0.909   0.045 * 
  
#  Group Month 2+Month 3+Month 4+Week 1+Week 2+Week 3+Week 4  #sps.  50 
#stat p.value    
#KO:K11785 0.999   0.001 ***
#  KO:K01551 0.989   0.001 ***
#  KO:K15735 0.989   0.001 ***
#  KO:K11930 0.986   0.008 ** 
#  KO:K00167 0.985   0.001 ***
#  KO:K07484 0.985   0.028 *  
#  KO:K08365 0.983   0.006 ** 
# KO:K03298 0.980   0.001 ***
#  KO:K19572 0.980   0.001 ***
#  KO:K07172 0.980   0.017 *  
#  KO:K07662 0.979   0.007 ** 
#  KO:K17725 0.978   0.023 *  
#  KO:K03435 0.977   0.023 *  
#  KO:K09946 0.977   0.020 *  
#  KO:K07979 0.977   0.024 *  
#  KO:K18285 0.974   0.012 *  
#  KO:K11784 0.973   0.023 *  
#  KO:K17947 0.970   0.008 ** 
#  KO:K16216 0.969   0.013 *  
#  KO:K08234 0.968   0.045 *  
#  KO:K01408 0.967   0.007 ** 
#  KO:K05847 0.966   0.037 *  
#  KO:K12984 0.962   0.012 *  
#  KO:K09898 0.959   0.019 *  
#  KO:K03762 0.959   0.041 *  
#  KO:K19267 0.955   0.004 ** 
#  KO:K07448 0.953   0.036 *  
#  KO:K06994 0.950   0.007 ** 
#  KO:K14941 0.950   0.045 *  
#  KO:K00370 0.949   0.027 *  
#  KO:K08306 0.949   0.008 ** 
#  KO:K15762 0.946   0.009 ** 
#  KO:K02848 0.943   0.022 *  
#  KO:K00372 0.941   0.028 *  
#  KO:K09792 0.939   0.039 *  
#  KO:K09954 0.937   0.037 *  
#  KO:K07030 0.935   0.009 ** 
#  KO:K02824 0.934   0.009 ** 
#  KO:K03284 0.931   0.010 ** 
#  KO:K12980 0.927   0.018 *  
#  KO:K18941 0.927   0.009 ** 
#  KO:K08989 0.918   0.017 *  
#  KO:K05020 0.915   0.041 *  
#  KO:K08982 0.912   0.018 *  
#  KO:K10622 0.901   0.036 *  
#  KO:K13622 0.885   0.033 *  
#  KO:K10856 0.884   0.036 *  
#  KO:K02584 0.872   0.034 *  
#  KO:K03817 0.872   0.044 *  
#  KO:K15981 0.849   0.042 *  
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

####Heatmap and z-score line plots of indicator KOs####
#Read in file
ind.funct <- read.csv("IndFunc.csv", header = T, sep = ",")

#Subset to indicator KOs
ind.funct.KOs <- ind.funct$KO.ID
SCOP_ind.funct.KOs <- bact.ko.rpkm2.spread2.t[, colnames(bact.ko.rpkm2.spread2.t) %in% ind.funct.KOs]

#Remove NA values
SCOP_ind.funct.KOs[is.na(SCOP_ind.funct.KOs)] <- 0

#Calculate z-scores
SCOP_ind.funct.KOs.zscore <- as.data.frame(scale(SCOP_ind.funct.KOs, center = T, scale = T))

#Format z-score dataframe
row.names(SCOP_ind.funct.KOs.zscore) <- c("August 4", "August 18", "September 1", "September 15", "September 22", "September 29", "October 4", "October 6", "October 8", "October 11", "October 13", "October 15", "October 18", "October 20", "October 22", "October 25", "October 27", "October 29", "November 3", "November 10", "November 17", "November 24", "December 1", "December 8", "December 15", "December 22", "December 29", "January 5", "January 12", "January 19", "January 26")

colnames(SCOP_ind.funct.KOs.zscore) <- paste(colnames(SCOP_ind.funct.KOs.zscore),"(",ind.funct$KO.symbol,")", sep = "") 

#Plot as heatmap
Heatmap(t(as.matrix(SCOP_ind.funct.KOs.zscore)), name = "z-score", column_title = "Date", column_title_side = "bottom", column_title_gp = gpar(fontsize = 12), column_names_gp = gpar(fontsize = 8), cluster_columns = F, show_column_names = T, clustering_distance_rows = "euclidean", clustering_method_rows = "complete", row_dend_width = unit(20, "mm"), row_title = "KEGG Orthologs (gene symbol)", row_title_side = "left", row_title_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 8))

#Plot as line plots
SCOP_ind.funct.KOs.zscore.sub <- SCOP_ind.funct.KOs.zscore
colnames(SCOP_ind.funct.KOs.zscore.sub) <- gsub(r"{\s*\([^\)]+\)}","",as.character(colnames(SCOP_ind.funct.KOs.zscore.sub)))

ind.funct.KOs.sub <- c("KO:K01551", "KO:K02402", "KO:K02759", "KO:K02760", "KO:K02761", "KO:K03805", "KO:K03893", "KO:K06375", "KO:K07692", "KO:K07699", "KO:K08365", "KO:K08969", "KO:K10622", "KO:K11614", "KO:K11931", "KO:K11936", "KO:K15735", "KO:K15762", "KO:K16216", "KO:K17725")

SCOP_ind.funct.KOs.zscore.sub2 <- SCOP_ind.funct.KOs.zscore.sub[, colnames(SCOP_ind.funct.KOs.zscore.sub) %in% ind.funct.KOs.sub, ]
SCOP_ind.funct.KOs.zscore.sub2$Date <- row.names(SCOP_ind.funct.KOs.zscore.sub2)
SCOP_ind.funct.KOs.zscore.sub2$Month <- c(rep(8, 2), rep(9, 4), rep(10, 12), rep(11, 4), rep(12, 5), rep(1, 4))
SCOP_ind.funct.KOs.zscore.sub2$Day <- c(4, 18, 1, 15, 22, 29, 4, 6, 8, 11, 13, 15, 18, 20, 22, 25, 27, 29, 3, 10, 17, 24, 1, 8, 15, 22, 29, 5, 12, 19, 26)
SCOP_ind.funct.KOs.zscore.sub2$Year <- c(rep(2021, 27), rep(2022, 4))
SCOP_ind.funct.KOs.zscore.sub2$Date_Combined <- paste(SCOP_ind.funct.KOs.zscore.sub2$Month, SCOP_ind.funct.KOs.zscore.sub2$Day, SCOP_ind.funct.KOs.zscore.sub2$Year, sep = "/")
SCOP_ind.funct.KOs.zscore.sub2$Correct_Date <- as.Date(SCOP_ind.funct.KOs.zscore.sub2$Date_Combined, "%m/%d/%Y")

SCOP_ind.funct.KOs.zscore.sub2.m1 <- SCOP_ind.funct.KOs.zscore.sub2[, c(13, 18, 19, 26)]
SCOP_ind.funct.KOs.zscore.sub2.m2 <- SCOP_ind.funct.KOs.zscore.sub2[, c(1, 7, 11, 26)]
SCOP_ind.funct.KOs.zscore.sub2.m3 <-SCOP_ind.funct.KOs.zscore.sub2[, c(6, 12, 17, 20, 26)]
SCOP_ind.funct.KOs.zscore.sub2.m4 <-SCOP_ind.funct.KOs.zscore.sub2[, c(2, 3:5, 8:10, 14:16, 26)]

SCOP_ind.funct.KOs.zscore.sub2.m1.g <- gather(SCOP_ind.funct.KOs.zscore.sub2.m1, Kegg.Ind, z.score, -Correct_Date)
SCOP_ind.funct.KOs.zscore.sub2.m2.g <- gather(SCOP_ind.funct.KOs.zscore.sub2.m2, Kegg.Ind, z.score, -Correct_Date)
SCOP_ind.funct.KOs.zscore.sub2.m3.g <- gather(SCOP_ind.funct.KOs.zscore.sub2.m3, Kegg.Ind, z.score, -Correct_Date)
SCOP_ind.funct.KOs.zscore.sub2.m4.g <- gather(SCOP_ind.funct.KOs.zscore.sub2.m4, Kegg.Ind, z.score, -Correct_Date)

#Toluene degradation, aromatic compound degradation, and aldehyde degradation
SCOP_ind.funct.plot1 <- ggplot(data = SCOP_ind.funct.KOs.zscore.sub2.m1.g) + 
  geom_hline(yintercept = 0, color= "black", linetype = "dashed", linewidth = 0.5) + 
  geom_line(aes(x = Correct_Date, y = z.score, color = Kegg.Ind), linetype = "solid", linewidth = 1) + 
  geom_point(aes(x = Correct_Date, y = z.score, color = Kegg.Ind), shape = 16, size = 3) + 
  theme_classic() + 
  labs(y = "z-score", x = "Date", title = "Aromatic Degradation") + 
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 16), legend.title = element_text(size = 18), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  scale_x_date(breaks = SCOP_ind.funct.KOs.zscore.sub2$Correct_Date, date_labels = "%Y-%m-%d") +
  expand_limits(y = c(-2, 4)) +
  scale_y_continuous(breaks = seq(-2, 4, 1)) + 
  scale_color_manual(name = "Kegg Ortholog (gene symbol):", values = c("grey50", "grey80", "black"), labels = c("KO:K10622 (cmtD/dhbB)", "KO:K15762 (tmoC/tbuB/touC)", "KO:K16216 (yueD)"))
SCOP_ind.funct.plot1

#Heavy metal
SCOP_ind.funct.plot2 <- ggplot(data = SCOP_ind.funct.KOs.zscore.sub2.m2.g) + 
  geom_hline(yintercept = 0, color= "black", linetype = "dashed", linewidth = 0.5) + 
  geom_line(aes(x = Correct_Date, y = z.score, color = Kegg.Ind), linetype = "solid", linewidth = 1) + 
  geom_point(aes(x = Correct_Date, y = z.score, color = Kegg.Ind), shape = 16, size = 3) + 
  theme_classic() + 
  labs(y = "z-score", x = "Date", title = "Heavy Metal Tolerance") + 
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 16), legend.title = element_text(size = 18), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  scale_x_date(breaks = SCOP_ind.funct.KOs.zscore.sub2$Correct_Date, date_labels = "%Y-%m-%d") +
  expand_limits(y = c(-2, 4)) +
  scale_y_continuous(breaks = seq(-2, 4, 1)) + 
  scale_color_manual(name = "Kegg Ortholog (gene symbol):", values = c("deeppink2", "deeppink4", "coral"), labels = c("KO:K01551 (arsA/ASNA1/GET3)", "KO:K03893 (arsB)", "KO:K08365 (merR)"))
SCOP_ind.funct.plot2

setwd("G:/My Drive/Dissertation/Presentation/")
png("SCOP Heavy Metal.png", width = 7, height = 4, units = "in", res = 600)
ggplot(data = SCOP_ind.funct.KOs.zscore.sub2.m2.g) + 
  geom_hline(yintercept = 0, color= "black", linetype = "dashed", linewidth = 0.5) + 
  geom_line(aes(x = Correct_Date, y = z.score, color = Kegg.Ind), linetype = "solid", linewidth = 1) + 
  geom_point(aes(x = Correct_Date, y = z.score, color = Kegg.Ind), shape = 16, size = 3) + 
  theme_classic() + 
  labs(y = "z-score", x = "Date", title = "Heavy Metal Tolerance") + 
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 9), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 12), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  scale_x_date(breaks = SCOP_ind.funct.KOs.zscore.sub2$Correct_Date, date_labels = "%Y-%m-%d") +
  expand_limits(y = c(-2, 4)) +
  scale_y_continuous(breaks = seq(-2, 4, 1)) + 
  scale_color_manual(name = "Kegg Ortholog (gene symbol):", values = c("deeppink2", "deeppink4", "coral"), labels = c("KO:K01551 (arsA/ASNA1/GET3)", "KO:K03893 (arsB)", "KO:K08365 (merR)"))
dev.off()

#Sulfur cycling and carbon starvation
SCOP_ind.funct.plot3 <- ggplot(data = SCOP_ind.funct.KOs.zscore.sub2.m3.g) + 
  geom_hline(yintercept = 0, color= "black", linetype = "dashed", linewidth = 0.5) + 
  geom_line(aes(x = Correct_Date, y = z.score, color = Kegg.Ind), linetype = "solid", linewidth = 1) + 
  geom_point(aes(x = Correct_Date, y = z.score, color = Kegg.Ind), shape = 16, size = 3) + 
  theme_classic() + 
  labs(y = "z-score", x = "Date", title = "Sulfur Metabolism and Carbon Starvation") + 
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 16), legend.title = element_text(size = 18), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  scale_x_date(breaks = SCOP_ind.funct.KOs.zscore.sub2$Correct_Date, date_labels = "%Y-%m-%d") +
  expand_limits(y = c(-2, 4)) + 
  scale_y_continuous(breaks = seq(-2, 4, 1)) + 
  scale_color_manual(name = "Kegg Ortholog (gene symbol):", values = c("plum1", "plum3", "darkseagreen4", "plum4"), labels = c("KO:K03805 (dsbG)", "KO:K08969 (mtnE/mtnV)", "KO:K15735 (csiR)", "KO:K17725 (ETHE1)"))
SCOP_ind.funct.plot3

#Biofilm formation, quorum sensing, and spore formation
SCOP_ind.funct.plot4 <- ggplot(data = SCOP_ind.funct.KOs.zscore.sub2.m4.g) + 
  geom_hline(yintercept = 0, color= "black", linetype = "dashed", linewidth = 0.5) + 
  geom_line(aes(x = Correct_Date, y = z.score, color = Kegg.Ind), linetype = "solid", linewidth = 1) + 
  geom_point(aes(x = Correct_Date, y = z.score, color = Kegg.Ind), shape = 16, size = 3) + 
  theme_classic() + 
  labs(y = "z-score", x = "Date", title = "Biofilm Formation, Quorum Sensing, and Spore Formation") + 
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 16), legend.title = element_text(size = 18), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  scale_x_date(breaks = SCOP_ind.funct.KOs.zscore.sub2$Correct_Date, date_labels = "%Y-%m-%d") +
  expand_limits(y = c(-2, 4)) + 
  scale_y_continuous(breaks = seq(-2, 4, 1)) + 
  scale_color_manual(name = "Kegg Ortholog (gene symbol):", values = c("peachpuff1", "paleturquoise2", "paleturquoise3", "paleturquoise4", "palevioletred1", "peachpuff3", "palevioletred3", "mediumorchid3", "coral1", "coral3"), labels = c("KO:K02402 (flhC)", "KO:K02759 (celC/chbA)", "KO:K02760 (celA/chbB)", "KO:K02761 (celB/chbC)", "KO:K06375 (spo0B)", "KO:K07692 (degU)", "KO:K07699 (spo0A)", "KO:K11614 (yufL/malK)", "KO:K11931 (pgaB)", "KO:K11936 (pgaC, icaA)"))
SCOP_ind.funct.plot4

#Export plots
png("SCOP Indicator KOs Z-scores.png", width = 12 , height = 16, units = "in", res = 600)
ggarrange(SCOP_ind.funct.plot1, SCOP_ind.funct.plot2, SCOP_ind.funct.plot3, SCOP_ind.funct.plot4, nrow = 4, ncol = 1, align = "hv", labels = c("a", "b", "c", "d", font.label = list(size = 24, color = "black", face = "bold", family = "sans")))
dev.off()

####Identifying Taxa that Indicator KOs are found in####
ind.funct.taxa <- read.csv("SCOP_Bacteria_KO_Taxa_RPKM.csv", header = T, row.names = 1, sep = ",")

#PAH degradation
#KO:K15762 - toluene
#KO:K10622 - aromatic compounds
#KO:K16216 - aromatic aldehyde
ind.funct.taxa.toluene <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K15762"), ]
sort(unique(ind.funct.taxa.toluene$Genus_ConsensusID))
ind.funct.taxa.toluene %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Candidatus Puniceispirillum

ind.funct.taxa.aro <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K10622"), ]
sort(unique(ind.funct.taxa.aro$Genus_ConsensusID))
ind.funct.taxa.aro %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Beijerinckia; Bosea 

ind.funct.taxa.aro2 <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K16216"), ]
sort(unique(ind.funct.taxa.aro2$Genus_ConsensusID))
taxa.aro2 <- ind.funct.taxa.aro2 %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Pontibacter; Owenweeksia
view(taxa.aro2)

#Heavy metal
#K01551: arsA
#K03893: arsB
#K08365: merR
ind.funct.taxa.arsA <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K01551"), ]
sort(unique(ind.funct.taxa.arsA$Genus_ConsensusID))
taxa.arsA <- ind.funct.taxa.arsA %>% group_by(Genus_ConsensusID) %>% summarize(count=n())
View(taxa.arsA) #Pseudohongiella

ind.funct.taxa.arsB <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K03893"), ]
sort(unique(ind.funct.taxa.arsB$Genus_ConsensusID))
ind.funct.taxa.arsB %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Marinococcus

ind.funct.taxa.merR <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K08365"), ]
sort(unique(ind.funct.taxa.merR$Genus_ConsensusID))
ind.funct.taxa.merR %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Sulfitobacter

#Carbon starvation
#K15735: csiR
ind.funct.taxa.csiR <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K15735"), ]
sort(unique(ind.funct.taxa.csiR$Genus_ConsensusID))
ind.funct.taxa.csiR %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Candidatus Thioglobus

#Sulfur pathways
#KO:K03805 - dsbG
#KO:K17725 - ETHE1
#KO:K08969 - mtnE, mtnV
ind.funct.taxa.dsbG <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K03805"), ]
sort(unique(ind.funct.taxa.dsbG$Genus_ConsensusID))
ind.funct.taxa.dsbG %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Halomonas

ind.funct.taxa.ETHE1 <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K17725"), ]
sort(unique(ind.funct.taxa.ETHE1$Genus_ConsensusID))
ind.funct.taxa.ETHE1 %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Candidatus Thioglobus

ind.funct.taxa.mtnE <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K03805"), ]
sort(unique(ind.funct.taxa.mtnE$Genus_ConsensusID))
ind.funct.taxa.mtnE %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Halomonas

#Spore formation
#KO:K06375 - spo0B
#KO:K07699 - spo0A
ind.funct.taxa.spo0B <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K06375"), ]
sort(unique(ind.funct.taxa.spo0B$Genus_ConsensusID))
ind.funct.taxa.spo0B %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Marinococcus

ind.funct.taxa.spo0A <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K07699"), ]
sort(unique(ind.funct.taxa.spo0A$Genus_ConsensusID))
ind.funct.taxa.spo0A %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Marinococcus

#Quorum sensing and biofilm formation
#K02759 - celC
#K02760 - celA
#K02761 - celB
#K11931 - pgaB
#K11936 - pgaC
#K11614 - yufL
#K02402 - flhC
#K07692 - degU

ind.funct.taxa.celC <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K02759"), ]
sort(unique(ind.funct.taxa.celC$Genus_ConsensusID))
ind.funct.taxa.celC %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Bacillus; Marinococcus

ind.funct.taxa.celA <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K02760"), ]
sort(unique(ind.funct.taxa.celA$Genus_ConsensusID))
ind.funct.taxa.celA %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Bacillus; Marinococcus 

ind.funct.taxa.celB <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K02761"), ]
sort(unique(ind.funct.taxa.celB$Genus_ConsensusID))
ind.funct.taxa.celB %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Marinococcus

ind.funct.taxa.pgaB <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K11931"), ]
sort(unique(ind.funct.taxa.pgaB$Genus_ConsensusID))
ind.funct.taxa.pgaB %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Halomonas

ind.funct.taxa.pgaC <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K11936"), ]
sort(unique(ind.funct.taxa.pgaC$Genus_ConsensusID))
ind.funct.taxa.pgaC %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Halomonas

ind.funct.taxa.yufL <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K11614"), ]
sort(unique(ind.funct.taxa.yufL$Genus_ConsensusID))
ind.funct.taxa.yufL %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Marinococcus

ind.funct.taxa.flhC <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K02402"), ]
sort(unique(ind.funct.taxa.flhC$Genus_ConsensusID))
ind.funct.taxa.flhC %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Halomonas

ind.funct.taxa.degU <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K07692"), ]
sort(unique(ind.funct.taxa.degU$Genus_ConsensusID))
ind.funct.taxa.degU %>% group_by(Genus_ConsensusID) %>% summarize(count=n()) #Marinococcus

#Nitrogen
#KO:K00370
ind.funct.taxa.degU <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K00370"), ]
sort(unique(ind.funct.taxa.degU$Genus_ConsensusID))
ind.funct.taxa.degU %>% group_by(Genus_ConsensusID) %>% summarize(count=n())
                                                          
#KO:K00372
ind.funct.taxa.degU <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K00372"), ]
sort(unique(ind.funct.taxa.degU$Genus_ConsensusID))
ind.funct.taxa.degU %>% group_by(Genus_ConsensusID) %>% summarize(count=n())

#KO:K02584
ind.funct.taxa.degU <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K02584"), ]
sort(unique(ind.funct.taxa.degU$Genus_ConsensusID))
ind.funct.taxa.degU %>% group_by(Genus_ConsensusID) %>% summarize(count=n())

#KO:K19345
ind.funct.taxa.degU <- ind.funct.taxa[which(ind.funct.taxa$ko == "KO:K19345"), ]
sort(unique(ind.funct.taxa.degU$Genus_ConsensusID))
ind.funct.taxa.degU %>% group_by(Genus_ConsensusID) %>% summarize(count=n())

####Changes in Kegg Modules over Time####
#Read in file
KEGG.info <- read.csv("KEGG_Modules_KO_List.csv", header = T, sep = ",")
SCOP.genes <- read.csv("SCOP_Bacteria_KO_Subset_RPKM.csv", header = T, row.names = 1, sep = ",")

#Check unique annotations and subset
SCOP.genes.sub <- SCOP.genes[, c(1, 8, 30)]

#Aggregate based on Sample ID and annotation
SCOP.genes.sub2 <- aggregate(RPKM~SampleID+ko,data=SCOP.genes.sub,FUN=sum)

#Aggregate based on KEGG Module
colnames(KEGG.info) <- c("Module.ID", "Module.Description", "KO", "ko", "KO.Description")

SCOP.kegg.module <- merge(SCOP.genes.sub2, KEGG.info, by = "ko")
SCOP.kegg.module.sub <- SCOP.kegg.module[, c(2:5)]
SCOP.kegg.module.sub2 <- aggregate(RPKM~SampleID+Module.Description, data = SCOP.kegg.module.sub, FUN = sum)

#Convert to wide format
SCOP.kegg.module.spread <- spread(SCOP.kegg.module.sub2, Module.Description, RPKM)

#Replace NAs with 0s 
SCOP.kegg.module.spread[is.na(SCOP.kegg.module.spread)] <- 0
row.names(SCOP.kegg.module.spread) <- SCOP.kegg.module.spread$SampleID

#Calculate z-scores of columns
SCOP.kegg.zscore <- as.data.frame(scale(SCOP.kegg.module.spread[-1], center = T, scale = T))

#Format dataframe
colnames(SCOP.kegg.zscore) <- c("Anoxygenic photosystem II", "Anthranilate degradation", "Assimilatory sulfate reduction", "Benzene degradation", "Benzoate degradation - 1", "Benzoate degradation - 2", "Benzoyl-CoA degradation", "Biphenyl degradation", "Carbazole degradation", "Catechol meta-cleavage", "Catechol ortho-cleavage", "Cobalamin biosynthesis", "Cumate degradation", "Cymene degradation", "Denitrification", "Dissimilatory nitrate reduction", "Dissimilatory sulfate reduction", "Homoprotocatechuate degradation", "Methane oxidation", "Methanogenesis, acetate => methane", "Methanogenesis, CO2 => methane", "Methanogenesis, methanol => methane", "Methanogenesis, methylamine => methane", "Naphthalene degradation", "Nicotinate degradation", "Nicotine degradation", "Nitrate assimilation", "Nitrification", "Photorespiration", "Photosystem I", "Photosystem II", "Phthalate degradation", "Sulfate-sulfur assimilation", "Thiosulfate oxidation", "Toluene degradation, anaerobic", "Toluene degradation") 

row.names(SCOP.kegg.zscore) <- c("August 4", "August 18", "September 1", "September 15", "September 22", "September 29", "October 4", "October 6", "October 8", "October 11", "October 13", "October 15", "October 18", "October 20", "October 22", "October 25", "October 27", "October 29", "November 3", "November 10", "November 17", "November 24", "December 1", "December 8", "December 15", "December 22", "December 29", "January 5", "January 12", "January 19", "January 26")

colnames(SCOP.kegg.zscore)
SCOP.kegg.zscore.nat <- SCOP.kegg.zscore[, c(3, 12, 15:18, 27, 28, 30, 31, 33, 34)]
SCOP.kegg.zscore.nat2 <- SCOP.kegg.zscore.nat
SCOP.kegg.zscore.nat2$Date <- SCOP.meta$Correct_Date

SCOP.kegg.zscore.nat2.m1 <- SCOP.kegg.zscore.nat2[, c(9, 10, 7, 13)]
SCOP.kegg.zscore.nat2.m1.g <- gather(SCOP.kegg.zscore.nat2.m1, Kegg.Module, z.score, -Date)

SCOP.kegg.zscore.nat2.m2 <- SCOP.kegg.zscore.nat2[, c(1, 5, 11:13)]
SCOP.kegg.zscore.nat2.m2.g <- gather(SCOP.kegg.zscore.nat2.m2, Kegg.Module, z.score, -Date)

colnames(SCOP.kegg.zscore)
SCOP.kegg.zscore.oil <- SCOP.kegg.zscore[, c(4, 8:9, 10, 11, 13, 24, 35, 36)]
SCOP.kegg.zscore.oil2 <- SCOP.kegg.zscore.oil
SCOP.kegg.zscore.oil2$Date <- SCOP.meta$Correct_Date
SCOP.kegg.zscore.oil2.m1 <- SCOP.kegg.zscore.oil2[, c(4, 5, 6, 9, 10)]
SCOP.kegg.zscore.oil2.m1.g <- gather(SCOP.kegg.zscore.oil2.m1, Kegg.Module, z.score, -Date)

SCOP.kegg.zscore.nat2.m1.g$Category <- "Module1"
SCOP.kegg.zscore.nat2.m2.g$Category <- "Module2"
SCOP.kegg.zscore.oil2.m1.g$Category <- "Module3"

SCOP.kegg.zscore.rbind <- rbind(SCOP.kegg.zscore.nat2.m1.g, SCOP.kegg.zscore.nat2.m2.g, SCOP.kegg.zscore.oil2.m1.g)

#Plot
SCOP.kegg.zscore.rbind.plot <- ggplot(SCOP.kegg.zscore.rbind, aes(x = Date, y = z.score, color = Kegg.Module))+
  geom_hline(yintercept = 0, color= "black", linetype = "dashed", linewidth = 0.5) + 
  geom_line(linetype = "solid", linewidth = 1) + 
  geom_point(shape = 16, size = 3)+
  facet_wrap(~Category, scales="free_y", ncol=1) + 
  theme_classic() + 
  labs(y = "z-score") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 16), legend.title = element_text(size = 18), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4))), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_x_date(breaks = SCOP.kegg.zscore.rbind$Date, date_labels = "%Y-%m-%d") + 
  expand_limits(y = c(-3, 4)) + 
  scale_y_continuous(breaks = seq(-3, 4, 1)) +
  scale_color_manual(name = "Kegg Module:", values = c("deeppink2", "black", "grey40", "grey60", "deeppink4", "turquoise4", "darkseagreen3", "darkseagreen4", "goldenrod1", "darkorange3", "grey80"))
SCOP.kegg.zscore.rbind.plot

getwd()
setwd("G:/My Drive/SCOP/Figures/Functional Patterns and Diversity/")

png("KeggModules_Time.png", units = "in", width = 12, height = 7, res = 600)
SCOP.kegg.zscore.rbind.plot
dev.off()

####Ridgeline Plots of Kegg Module Residuals####
#Import residuals of KOs
ko.resid <- read.csv("KO_resid.csv", header = T, row.names = 1, sep = ",")

#Import list of KOs in each Kegg Module of interest
kegg.module.kos <- read.csv("KEGG_Modules_KO_List_sub.csv", header = T, sep = ",")

#Subset KO residuals to those that are in Kegg Modules of interest
kegg.module.ko.resid <- ko.resid[which(row.names(ko.resid) %in% kegg.module.kos$KO.ID), ]

#Format dataframe 
kegg.module.ko.resid.t <- as.data.frame(t(kegg.module.ko.resid))
kegg.module.ko.resid.t2 <- kegg.module.ko.resid.t[c(1, 2, 15:26, 38:47, 62:73, 81:90, 103:114, 127:140, 156:165, 178:187, 204:216, 226:267), ]
kegg.module.ko.resid.t2$SampleID <- row.names(kegg.module.ko.resid.t2)
kegg.module.ko.resid.t2$Year <- c(rep(2011, 12), rep(2012, 11), rep(2013, 11), rep(2014, 10), rep(2015, 13), rep(2016, 13), rep(2017, 10), rep(2018, 8), rep(2019, 15), rep(2020, 13), rep(2021, 27), rep(2022, 4))
kegg.module.ko.resid.t2$Month <- c(1, 1, 8, 8, 9, 9, 10, 10, 11, 11, 11, 12, 1, 1, 8, 8, 9, 10, 10, 10, 11, 11, 12, 1, 8, 8, 9, 10, 10, 10, 11, 11, 12, 12, 1, 1, 8, 8, 9, 10, 11, 11, 12, 12, 1, 1, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 12, 1, 8, 8, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 1, 1, 8, 8, 9, 9, 10, 10, 10, 11, 1, 1, 8, 8, 9, 9, 10, 12, 1, 1, 1, 1, 8, 8, 8, 9, 10, 10, 10, 10, 10, 11, 12, 1, 1, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 12, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 12, 1, 1, 1, 1)
kegg.module.ko.resid.t2$Day <- c(12, 26, 17, 24, 7, 21, 12, 19, 3, 16, 30, 14, 11, 25, 8, 22, 12, 3, 10, 24, 14, 28, 5, 9, 7, 14, 11, 2, 16, 30, 13, 27, 18, 24, 8, 22, 6, 20, 17, 1, 12, 26, 10, 17, 6, 14, 5, 19, 2, 16, 14, 28, 11, 25, 9, 13, 23, 20, 3, 7, 17, 31, 14, 28, 12, 26, 9, 23, 7, 21, 18, 25, 9, 23, 6, 20, 11, 18, 25, 15, 12, 17, 8, 22, 5, 19, 2, 5, 9, 16, 23, 30, 6, 7, 21, 4, 2, 9, 16, 23, 30, 13, 18, 8, 22, 12, 26, 9, 23, 7, 21, 18, 25, 2, 9, 16, 4, 18, 1, 15, 22, 29, 4, 6, 8, 11, 13, 15, 18, 20, 22, 25, 27, 29, 3, 10, 17, 24, 1, 8, 15, 22, 29, 5, 12, 19, 26)

#Put date in date format
kegg.module.ko.resid.t2$Date_Combined <- paste(kegg.module.ko.resid.t2$Month, kegg.module.ko.resid.t2$Day, kegg.module.ko.resid.t2$Year, sep = "/")
kegg.module.ko.resid.t2$Correct_Date <- as.Date(kegg.module.ko.resid.t2$Date_Combined, "%m/%d/%Y")

#Subset to each month
kegg.ko.resid.Aug <- kegg.module.ko.resid.t2[which(kegg.module.ko.resid.t2$Month == 8), ]
kegg.ko.resid.Sept <- kegg.module.ko.resid.t2[which(kegg.module.ko.resid.t2$Month == 9), ]
kegg.ko.resid.Oct <- kegg.module.ko.resid.t2[which(kegg.module.ko.resid.t2$Month == 10), ]
kegg.ko.resid.Nov <- kegg.module.ko.resid.t2[which(kegg.module.ko.resid.t2$Month == 11), ]
kegg.ko.resid.Dec <- kegg.module.ko.resid.t2[which(kegg.module.ko.resid.t2$Month == 12), ]
kegg.ko.resid.Jan <- kegg.module.ko.resid.t2[which(kegg.module.ko.resid.t2$Month == 1), ]

#Convert to long format
kegg.ko.resid.Aug.g <- gather(kegg.ko.resid.Aug, KO.ID, KO_Resid, -SampleID, -Correct_Date, -Date_Combined, -Day, -Month, -Year)
kegg.ko.resid.Sept.g <- gather(kegg.ko.resid.Sept, KO.ID, KO_Resid, -SampleID, -Correct_Date, -Date_Combined, -Day, -Month, -Year)
kegg.ko.resid.Oct.g <- gather(kegg.ko.resid.Oct, KO.ID, KO_Resid, -SampleID, -Correct_Date, -Date_Combined, -Day, -Month, -Year)
kegg.ko.resid.Nov.g <- gather(kegg.ko.resid.Nov, KO.ID, KO_Resid, -SampleID, -Correct_Date, -Date_Combined, -Day, -Month, -Year)
kegg.ko.resid.Dec.g <- gather(kegg.ko.resid.Dec, KO.ID, KO_Resid, -SampleID, -Correct_Date, -Date_Combined, -Day, -Month, -Year)
kegg.ko.resid.Jan.g <- gather(kegg.ko.resid.Jan, KO.ID, KO_Resid, -SampleID, -Correct_Date, -Date_Combined, -Day, -Month, -Year)

#Combine with metadata
kegg.ko.resid.Aug.g.meta <- merge(kegg.ko.resid.Aug.g, kegg.module.kos, by = "KO.ID")
kegg.ko.resid.Sept.g.meta <- merge(kegg.ko.resid.Sept.g, kegg.module.kos, by = "KO.ID")
kegg.ko.resid.Oct.g.meta <- merge(kegg.ko.resid.Oct.g, kegg.module.kos, by = "KO.ID")
kegg.ko.resid.Nov.g.meta <- merge(kegg.ko.resid.Nov.g, kegg.module.kos, by = "KO.ID")
kegg.ko.resid.Dec.g.meta <- merge(kegg.ko.resid.Dec.g, kegg.module.kos, by = "KO.ID")
kegg.ko.resid.Jan.g.meta <- merge(kegg.ko.resid.Jan.g, kegg.module.kos, by = "KO.ID")

#Make individual dataframes for each Kegg Module for each month - recombine according to broader categories
#August
kegg.ko.resid.Aug.g.meta.ps1 <- kegg.ko.resid.Aug.g.meta[which(kegg.ko.resid.Aug.g.meta$Module.Description == "Photosystem I"), ]
kegg.ko.resid.Aug.g.meta.ps2 <- kegg.ko.resid.Aug.g.meta[which(kegg.ko.resid.Aug.g.meta$Module.Description == "Photosystem II"), ]
kegg.ko.resid.Aug.g.meta.no3 <- kegg.ko.resid.Aug.g.meta[which(kegg.ko.resid.Aug.g.meta$Module.Description == "Nitrate assimilation"), ]
kegg.ko.resid.Aug.g.meta.tol <- kegg.ko.resid.Aug.g.meta[which(kegg.ko.resid.Aug.g.meta$Module.Description == "Toluene degradation, toluene => benzoate"), ]
kegg.ko.resid.Aug.g.meta.cat1 <- kegg.ko.resid.Aug.g.meta[which(kegg.ko.resid.Aug.g.meta$Module.Description == "Catechol meta-cleavage, catechol => acetyl-CoA / 4-methylcatechol => propanoyl-CoA"), ]
kegg.ko.resid.Aug.g.meta.cumate <- kegg.ko.resid.Aug.g.meta[which(kegg.ko.resid.Aug.g.meta$Module.Description == "Cumate degradation, p-cumate => 2-oxopent-4-enoate + 2-methylpropanoate"), ]
kegg.ko.resid.Aug.g.meta.cat2 <- kegg.ko.resid.Aug.g.meta[which(kegg.ko.resid.Aug.g.meta$Module.Description == "Catechol ortho-cleavage, catechol => 3-oxoadipate"), ]
kegg.ko.resid.Aug.g.meta.sulf1 <- kegg.ko.resid.Aug.g.meta[which(kegg.ko.resid.Aug.g.meta$Module.Description == "Assimilatory sulfate reduction, sulfate => H2S"), ]
kegg.ko.resid.Aug.g.meta.sulf2 <- kegg.ko.resid.Aug.g.meta[which(kegg.ko.resid.Aug.g.meta$Module.Description == "Dissimilatory sulfate reduction, sulfate => H2S"), ]
kegg.ko.resid.Aug.g.meta.sulf3 <- kegg.ko.resid.Aug.g.meta[which(kegg.ko.resid.Aug.g.meta$Module.Description == "Sulfate-sulfur assimilation"), ]
kegg.ko.resid.Aug.g.meta.sulf4 <- kegg.ko.resid.Aug.g.meta[which(kegg.ko.resid.Aug.g.meta$Module.Description == "Thiosulfate oxidation by SOX complex, thiosulfate => sulfate"), ]

kegg.ko.resid.Aug.g.meta.c <- rbind(kegg.ko.resid.Aug.g.meta.ps1, kegg.ko.resid.Aug.g.meta.ps2, kegg.ko.resid.Aug.g.meta.no3)
kegg.ko.resid.Aug.g.meta.c2 <- rbind(kegg.ko.resid.Aug.g.meta.tol, kegg.ko.resid.Aug.g.meta.cat1, kegg.ko.resid.Aug.g.meta.cat2, kegg.ko.resid.Aug.g.meta.cumate)
kegg.ko.resid.Aug.g.meta.c3 <- rbind(kegg.ko.resid.Aug.g.meta.sulf1, kegg.ko.resid.Aug.g.meta.sulf2, kegg.ko.resid.Aug.g.meta.sulf3, kegg.ko.resid.Aug.g.meta.sulf4)

kegg.ko.resid.Aug.g.meta.c4 <- rbind(kegg.ko.resid.Aug.g.meta.c, kegg.ko.resid.Aug.g.meta.c2, kegg.ko.resid.Aug.g.meta.c3)
kegg.ko.resid.Aug.g.meta.median <- kegg.ko.resid.Aug.g.meta.c4 %>% group_by(SampleID, Module.Description) %>% summarize(Median_KO_Resid = median(KO_Resid, na.rm = TRUE))

kegg.ko.resid.Aug.meta <- kegg.ko.resid.Aug[, c(139:144)]
kegg.ko.resid.Aug.g.meta.median2 <- merge(kegg.ko.resid.Aug.g.meta.median, kegg.ko.resid.Aug.meta, by = "SampleID")
kegg.ko.resid.Aug.g.meta.median2$Month_Day <- paste(kegg.ko.resid.Aug.g.meta.median2$Month, kegg.ko.resid.Aug.g.meta.median2$Day, sep = "/")


#September
kegg.ko.resid.Sept.g.meta.ps1 <- kegg.ko.resid.Sept.g.meta[which(kegg.ko.resid.Sept.g.meta$Module.Description == "Photosystem I"), ]
kegg.ko.resid.Sept.g.meta.ps2 <- kegg.ko.resid.Sept.g.meta[which(kegg.ko.resid.Sept.g.meta$Module.Description == "Photosystem II"), ]
kegg.ko.resid.Sept.g.meta.no3 <- kegg.ko.resid.Sept.g.meta[which(kegg.ko.resid.Sept.g.meta$Module.Description == "Nitrate assimilation"), ]
kegg.ko.resid.Sept.g.meta.tol <- kegg.ko.resid.Sept.g.meta[which(kegg.ko.resid.Sept.g.meta$Module.Description == "Toluene degradation, toluene => benzoate"), ]
kegg.ko.resid.Sept.g.meta.cat1 <- kegg.ko.resid.Sept.g.meta[which(kegg.ko.resid.Sept.g.meta$Module.Description == "Catechol meta-cleavage, catechol => acetyl-CoA / 4-methylcatechol => propanoyl-CoA"), ]
kegg.ko.resid.Sept.g.meta.cat2 <- kegg.ko.resid.Sept.g.meta[which(kegg.ko.resid.Sept.g.meta$Module.Description == "Catechol ortho-cleavage, catechol => 3-oxoadipate"), ]
kegg.ko.resid.Sept.g.meta.cumate <- kegg.ko.resid.Sept.g.meta[which(kegg.ko.resid.Sept.g.meta$Module.Description == "Cumate degradation, p-cumate => 2-oxopent-4-enoate + 2-methylpropanoate"), ]
kegg.ko.resid.Sept.g.meta.sulf1 <- kegg.ko.resid.Sept.g.meta[which(kegg.ko.resid.Sept.g.meta$Module.Description == "Assimilatory sulfate reduction, sulfate => H2S"), ]
kegg.ko.resid.Sept.g.meta.sulf2 <- kegg.ko.resid.Sept.g.meta[which(kegg.ko.resid.Sept.g.meta$Module.Description == "Dissimilatory sulfate reduction, sulfate => H2S"), ]
kegg.ko.resid.Sept.g.meta.sulf3 <- kegg.ko.resid.Sept.g.meta[which(kegg.ko.resid.Sept.g.meta$Module.Description == "Sulfate-sulfur assimilation"), ]
kegg.ko.resid.Sept.g.meta.sulf4 <- kegg.ko.resid.Sept.g.meta[which(kegg.ko.resid.Sept.g.meta$Module.Description == "Thiosulfate oxidation by SOX complex, thiosulfate => sulfate"), ]

kegg.ko.resid.Sept.g.meta.c <- rbind(kegg.ko.resid.Sept.g.meta.ps1, kegg.ko.resid.Sept.g.meta.ps2, kegg.ko.resid.Sept.g.meta.no3)
kegg.ko.resid.Sept.g.meta.c2 <- rbind(kegg.ko.resid.Sept.g.meta.tol, kegg.ko.resid.Sept.g.meta.cat1, kegg.ko.resid.Sept.g.meta.cat2, kegg.ko.resid.Sept.g.meta.cumate)
kegg.ko.resid.Sept.g.meta.c3 <- rbind(kegg.ko.resid.Sept.g.meta.sulf1, kegg.ko.resid.Sept.g.meta.sulf2, kegg.ko.resid.Sept.g.meta.sulf3, kegg.ko.resid.Sept.g.meta.sulf4)

kegg.ko.resid.Sept.g.meta.c4 <- rbind(kegg.ko.resid.Sept.g.meta.c, kegg.ko.resid.Sept.g.meta.c2, kegg.ko.resid.Sept.g.meta.c3)
kegg.ko.resid.Sept.g.meta.median <- kegg.ko.resid.Sept.g.meta.c4 %>% group_by(SampleID, Module.Description) %>% summarize(Median_KO_Resid = median(KO_Resid, na.rm = TRUE))

kegg.ko.resid.Sept.meta <- kegg.ko.resid.Sept[, c(139:144)]
kegg.ko.resid.Sept.g.meta.median2 <- merge(kegg.ko.resid.Sept.g.meta.median, kegg.ko.resid.Sept.meta, by = "SampleID")
kegg.ko.resid.Sept.g.meta.median2$Month_Day <- paste(kegg.ko.resid.Sept.g.meta.median2$Month, kegg.ko.resid.Sept.g.meta.median2$Day, sep = "/")


#October
kegg.ko.resid.Oct.g.meta.ps1 <- kegg.ko.resid.Oct.g.meta[which(kegg.ko.resid.Oct.g.meta$Module.Description == "Photosystem I"), ]
kegg.ko.resid.Oct.g.meta.ps2 <- kegg.ko.resid.Oct.g.meta[which(kegg.ko.resid.Oct.g.meta$Module.Description == "Photosystem II"), ]
kegg.ko.resid.Oct.g.meta.no3 <- kegg.ko.resid.Oct.g.meta[which(kegg.ko.resid.Oct.g.meta$Module.Description == "Nitrate assimilation"), ]
kegg.ko.resid.Oct.g.meta.tol <- kegg.ko.resid.Oct.g.meta[which(kegg.ko.resid.Oct.g.meta$Module.Description == "Toluene degradation, toluene => benzoate"), ]
kegg.ko.resid.Oct.g.meta.cat1 <- kegg.ko.resid.Oct.g.meta[which(kegg.ko.resid.Oct.g.meta$Module.Description == "Catechol meta-cleavage, catechol => acetyl-CoA / 4-methylcatechol => propanoyl-CoA"), ]
kegg.ko.resid.Oct.g.meta.cat2 <- kegg.ko.resid.Oct.g.meta[which(kegg.ko.resid.Oct.g.meta$Module.Description == "Catechol ortho-cleavage, catechol => 3-oxoadipate"), ]
kegg.ko.resid.Oct.g.meta.cumate <- kegg.ko.resid.Oct.g.meta[which(kegg.ko.resid.Oct.g.meta$Module.Description == "Cumate degradation, p-cumate => 2-oxopent-4-enoate + 2-methylpropanoate"), ]
kegg.ko.resid.Oct.g.meta.sulf1 <- kegg.ko.resid.Oct.g.meta[which(kegg.ko.resid.Oct.g.meta$Module.Description == "Assimilatory sulfate reduction, sulfate => H2S"), ]
kegg.ko.resid.Oct.g.meta.sulf2 <- kegg.ko.resid.Oct.g.meta[which(kegg.ko.resid.Oct.g.meta$Module.Description == "Dissimilatory sulfate reduction, sulfate => H2S"), ]
kegg.ko.resid.Oct.g.meta.sulf3 <- kegg.ko.resid.Oct.g.meta[which(kegg.ko.resid.Oct.g.meta$Module.Description == "Sulfate-sulfur assimilation"), ]
kegg.ko.resid.Oct.g.meta.sulf4 <- kegg.ko.resid.Oct.g.meta[which(kegg.ko.resid.Oct.g.meta$Module.Description == "Thiosulfate oxidation by SOX complex, thiosulfate => sulfate"), ]

kegg.ko.resid.Oct.g.meta.c <- rbind(kegg.ko.resid.Oct.g.meta.ps1, kegg.ko.resid.Oct.g.meta.ps2, kegg.ko.resid.Oct.g.meta.no3)
kegg.ko.resid.Oct.g.meta.c2 <- rbind(kegg.ko.resid.Oct.g.meta.tol, kegg.ko.resid.Oct.g.meta.cat1, kegg.ko.resid.Oct.g.meta.cat2, kegg.ko.resid.Oct.g.meta.cumate)
kegg.ko.resid.Oct.g.meta.c3 <- rbind(kegg.ko.resid.Oct.g.meta.sulf1, kegg.ko.resid.Oct.g.meta.sulf2, kegg.ko.resid.Oct.g.meta.sulf3, kegg.ko.resid.Oct.g.meta.sulf4)

kegg.ko.resid.Oct.g.meta.c4 <- rbind(kegg.ko.resid.Oct.g.meta.c, kegg.ko.resid.Oct.g.meta.c2, kegg.ko.resid.Oct.g.meta.c3)
kegg.ko.resid.Oct.g.meta.median <- kegg.ko.resid.Oct.g.meta.c4 %>% group_by(SampleID, Module.Description) %>% summarize(Median_KO_Resid = median(KO_Resid, na.rm = TRUE))

kegg.ko.resid.Oct.meta <- kegg.ko.resid.Oct[, c(139:144)]
kegg.ko.resid.Oct.g.meta.median2 <- merge(kegg.ko.resid.Oct.g.meta.median, kegg.ko.resid.Oct.meta, by = "SampleID")
kegg.ko.resid.Oct.g.meta.median2$Month_Day <- paste(kegg.ko.resid.Oct.g.meta.median2$Month, kegg.ko.resid.Oct.g.meta.median2$Day, sep = "/")


#November
kegg.ko.resid.Nov.g.meta.ps1 <- kegg.ko.resid.Nov.g.meta[which(kegg.ko.resid.Nov.g.meta$Module.Description == "Photosystem I"), ]
kegg.ko.resid.Nov.g.meta.ps2 <- kegg.ko.resid.Nov.g.meta[which(kegg.ko.resid.Nov.g.meta$Module.Description == "Photosystem II"), ]
kegg.ko.resid.Nov.g.meta.no3 <- kegg.ko.resid.Nov.g.meta[which(kegg.ko.resid.Nov.g.meta$Module.Description == "Nitrate assimilation"), ]
kegg.ko.resid.Nov.g.meta.tol <- kegg.ko.resid.Nov.g.meta[which(kegg.ko.resid.Nov.g.meta$Module.Description == "Toluene degradation, toluene => benzoate"), ]
kegg.ko.resid.Nov.g.meta.cat1 <- kegg.ko.resid.Nov.g.meta[which(kegg.ko.resid.Nov.g.meta$Module.Description == "Catechol meta-cleavage, catechol => acetyl-CoA / 4-methylcatechol => propanoyl-CoA"), ]
kegg.ko.resid.Nov.g.meta.cat2 <- kegg.ko.resid.Nov.g.meta[which(kegg.ko.resid.Nov.g.meta$Module.Description == "Catechol ortho-cleavage, catechol => 3-oxoadipate"), ]
kegg.ko.resid.Nov.g.meta.cumate <- kegg.ko.resid.Nov.g.meta[which(kegg.ko.resid.Nov.g.meta$Module.Description == "Cumate degradation, p-cumate => 2-oxopent-4-enoate + 2-methylpropanoate"), ]
kegg.ko.resid.Nov.g.meta.sulf1 <- kegg.ko.resid.Nov.g.meta[which(kegg.ko.resid.Nov.g.meta$Module.Description == "Assimilatory sulfate reduction, sulfate => H2S"), ]
kegg.ko.resid.Nov.g.meta.sulf2 <- kegg.ko.resid.Nov.g.meta[which(kegg.ko.resid.Nov.g.meta$Module.Description == "Dissimilatory sulfate reduction, sulfate => H2S"), ]
kegg.ko.resid.Nov.g.meta.sulf3 <- kegg.ko.resid.Nov.g.meta[which(kegg.ko.resid.Nov.g.meta$Module.Description == "Sulfate-sulfur assimilation"), ]
kegg.ko.resid.Nov.g.meta.sulf4 <- kegg.ko.resid.Nov.g.meta[which(kegg.ko.resid.Nov.g.meta$Module.Description == "Thiosulfate oxidation by SOX complex, thiosulfate => sulfate"), ]

kegg.ko.resid.Nov.g.meta.c <- rbind(kegg.ko.resid.Nov.g.meta.ps1, kegg.ko.resid.Nov.g.meta.ps2, kegg.ko.resid.Nov.g.meta.no3)
kegg.ko.resid.Nov.g.meta.c2 <- rbind(kegg.ko.resid.Nov.g.meta.tol, kegg.ko.resid.Nov.g.meta.cat1, kegg.ko.resid.Nov.g.meta.cat2, kegg.ko.resid.Nov.g.meta.cumate)
kegg.ko.resid.Nov.g.meta.c3 <- rbind(kegg.ko.resid.Nov.g.meta.sulf1, kegg.ko.resid.Nov.g.meta.sulf2, kegg.ko.resid.Nov.g.meta.sulf3, kegg.ko.resid.Nov.g.meta.sulf4)

kegg.ko.resid.Nov.g.meta.c4 <- rbind(kegg.ko.resid.Nov.g.meta.c, kegg.ko.resid.Nov.g.meta.c2, kegg.ko.resid.Nov.g.meta.c3)
kegg.ko.resid.Nov.g.meta.median <- kegg.ko.resid.Nov.g.meta.c4 %>% group_by(SampleID, Module.Description) %>% summarize(Median_KO_Resid = median(KO_Resid, na.rm = TRUE))

kegg.ko.resid.Nov.meta <- kegg.ko.resid.Nov[, c(139:144)]
kegg.ko.resid.Nov.g.meta.median2 <- merge(kegg.ko.resid.Nov.g.meta.median, kegg.ko.resid.Nov.meta, by = "SampleID")
kegg.ko.resid.Nov.g.meta.median2$Month_Day <- paste(kegg.ko.resid.Nov.g.meta.median2$Month, kegg.ko.resid.Nov.g.meta.median2$Day, sep = "/")


#December
kegg.ko.resid.Dec.g.meta.ps1 <- kegg.ko.resid.Dec.g.meta[which(kegg.ko.resid.Dec.g.meta$Module.Description == "Photosystem I"), ]
kegg.ko.resid.Dec.g.meta.ps2 <- kegg.ko.resid.Dec.g.meta[which(kegg.ko.resid.Dec.g.meta$Module.Description == "Photosystem II"), ]
kegg.ko.resid.Dec.g.meta.no3 <- kegg.ko.resid.Dec.g.meta[which(kegg.ko.resid.Dec.g.meta$Module.Description == "Nitrate assimilation"), ]
kegg.ko.resid.Dec.g.meta.tol <- kegg.ko.resid.Dec.g.meta[which(kegg.ko.resid.Dec.g.meta$Module.Description == "Toluene degradation, toluene => benzoate"), ]
kegg.ko.resid.Dec.g.meta.cat1 <- kegg.ko.resid.Dec.g.meta[which(kegg.ko.resid.Dec.g.meta$Module.Description == "Catechol meta-cleavage, catechol => acetyl-CoA / 4-methylcatechol => propanoyl-CoA"), ]
kegg.ko.resid.Dec.g.meta.cat2 <- kegg.ko.resid.Dec.g.meta[which(kegg.ko.resid.Dec.g.meta$Module.Description == "Catechol ortho-cleavage, catechol => 3-oxoadipate"), ]
kegg.ko.resid.Dec.g.meta.cumate <- kegg.ko.resid.Dec.g.meta[which(kegg.ko.resid.Dec.g.meta$Module.Description == "Cumate degradation, p-cumate => 2-oxopent-4-enoate + 2-methylpropanoate"), ]
kegg.ko.resid.Dec.g.meta.sulf1 <- kegg.ko.resid.Dec.g.meta[which(kegg.ko.resid.Dec.g.meta$Module.Description == "Assimilatory sulfate reduction, sulfate => H2S"), ]
kegg.ko.resid.Dec.g.meta.sulf2 <- kegg.ko.resid.Dec.g.meta[which(kegg.ko.resid.Dec.g.meta$Module.Description == "Dissimilatory sulfate reduction, sulfate => H2S"), ]
kegg.ko.resid.Dec.g.meta.sulf3 <- kegg.ko.resid.Dec.g.meta[which(kegg.ko.resid.Dec.g.meta$Module.Description == "Sulfate-sulfur assimilation"), ]
kegg.ko.resid.Dec.g.meta.sulf4 <- kegg.ko.resid.Dec.g.meta[which(kegg.ko.resid.Dec.g.meta$Module.Description == "Thiosulfate oxidation by SOX complex, thiosulfate => sulfate"), ]

kegg.ko.resid.Dec.g.meta.c <- rbind(kegg.ko.resid.Dec.g.meta.ps1, kegg.ko.resid.Dec.g.meta.ps2, kegg.ko.resid.Dec.g.meta.no3)
kegg.ko.resid.Dec.g.meta.c2 <- rbind(kegg.ko.resid.Dec.g.meta.tol, kegg.ko.resid.Dec.g.meta.cat1, kegg.ko.resid.Dec.g.meta.cat2, kegg.ko.resid.Dec.g.meta.cumate)
kegg.ko.resid.Dec.g.meta.c3 <- rbind(kegg.ko.resid.Dec.g.meta.sulf1, kegg.ko.resid.Dec.g.meta.sulf2, kegg.ko.resid.Dec.g.meta.sulf3, kegg.ko.resid.Dec.g.meta.sulf4)

kegg.ko.resid.Dec.g.meta.c4 <- rbind(kegg.ko.resid.Dec.g.meta.c, kegg.ko.resid.Dec.g.meta.c2, kegg.ko.resid.Dec.g.meta.c3)
kegg.ko.resid.Dec.g.meta.median <- kegg.ko.resid.Dec.g.meta.c4 %>% group_by(SampleID, Module.Description) %>% summarize(Median_KO_Resid = median(KO_Resid, na.rm = TRUE))

kegg.ko.resid.Dec.meta <- kegg.ko.resid.Dec[, c(139:144)]
kegg.ko.resid.Dec.g.meta.median2 <- merge(kegg.ko.resid.Dec.g.meta.median, kegg.ko.resid.Dec.meta, by = "SampleID")
kegg.ko.resid.Dec.g.meta.median2$Month_Day <- paste(kegg.ko.resid.Dec.g.meta.median2$Month, kegg.ko.resid.Dec.g.meta.median2$Day, sep = "/")


#January
kegg.ko.resid.Jan.g.meta.ps1 <- kegg.ko.resid.Jan.g.meta[which(kegg.ko.resid.Jan.g.meta$Module.Description == "Photosystem I"), ]
kegg.ko.resid.Jan.g.meta.ps2 <- kegg.ko.resid.Jan.g.meta[which(kegg.ko.resid.Jan.g.meta$Module.Description == "Photosystem II"), ]
kegg.ko.resid.Jan.g.meta.no3 <- kegg.ko.resid.Jan.g.meta[which(kegg.ko.resid.Jan.g.meta$Module.Description == "Nitrate assimilation"), ]
kegg.ko.resid.Jan.g.meta.tol <- kegg.ko.resid.Jan.g.meta[which(kegg.ko.resid.Jan.g.meta$Module.Description == "Toluene degradation, toluene => benzoate"), ]
kegg.ko.resid.Jan.g.meta.cat1 <- kegg.ko.resid.Jan.g.meta[which(kegg.ko.resid.Jan.g.meta$Module.Description == "Catechol meta-cleavage, catechol => acetyl-CoA / 4-methylcatechol => propanoyl-CoA"), ]
kegg.ko.resid.Jan.g.meta.cat2 <- kegg.ko.resid.Jan.g.meta[which(kegg.ko.resid.Jan.g.meta$Module.Description == "Catechol ortho-cleavage, catechol => 3-oxoadipate"), ]
kegg.ko.resid.Jan.g.meta.cumate <- kegg.ko.resid.Jan.g.meta[which(kegg.ko.resid.Jan.g.meta$Module.Description == "Cumate degradation, p-cumate => 2-oxopent-4-enoate + 2-methylpropanoate"), ]
kegg.ko.resid.Jan.g.meta.sulf1 <- kegg.ko.resid.Jan.g.meta[which(kegg.ko.resid.Jan.g.meta$Module.Description == "Assimilatory sulfate reduction, sulfate => H2S"), ]
kegg.ko.resid.Jan.g.meta.sulf2 <- kegg.ko.resid.Jan.g.meta[which(kegg.ko.resid.Jan.g.meta$Module.Description == "Dissimilatory sulfate reduction, sulfate => H2S"), ]
kegg.ko.resid.Jan.g.meta.sulf3 <- kegg.ko.resid.Jan.g.meta[which(kegg.ko.resid.Jan.g.meta$Module.Description == "Sulfate-sulfur assimilation"), ]
kegg.ko.resid.Jan.g.meta.sulf4 <- kegg.ko.resid.Jan.g.meta[which(kegg.ko.resid.Jan.g.meta$Module.Description == "Thiosulfate oxidation by SOX complex, thiosulfate => sulfate"), ]

kegg.ko.resid.Jan.g.meta.c <- rbind(kegg.ko.resid.Jan.g.meta.ps1, kegg.ko.resid.Jan.g.meta.ps2, kegg.ko.resid.Jan.g.meta.no3)
kegg.ko.resid.Jan.g.meta.c2 <- rbind(kegg.ko.resid.Jan.g.meta.tol, kegg.ko.resid.Jan.g.meta.cat1, kegg.ko.resid.Jan.g.meta.cat2, kegg.ko.resid.Jan.g.meta.cumate)
kegg.ko.resid.Jan.g.meta.c3 <- rbind(kegg.ko.resid.Jan.g.meta.sulf1, kegg.ko.resid.Jan.g.meta.sulf2, kegg.ko.resid.Jan.g.meta.sulf3, kegg.ko.resid.Jan.g.meta.sulf4)

kegg.ko.resid.Jan.g.meta.c4 <- rbind(kegg.ko.resid.Jan.g.meta.c, kegg.ko.resid.Jan.g.meta.c2, kegg.ko.resid.Jan.g.meta.c3)
kegg.ko.resid.Jan.g.meta.median <- kegg.ko.resid.Jan.g.meta.c4 %>% group_by(SampleID, Module.Description) %>% summarize(Median_KO_Resid = median(KO_Resid, na.rm = TRUE))

kegg.ko.resid.Jan.meta <- kegg.ko.resid.Jan[, c(139:144)]
kegg.ko.resid.Jan.g.meta.median2 <- merge(kegg.ko.resid.Jan.g.meta.median, kegg.ko.resid.Jan.meta, by = "SampleID")
kegg.ko.resid.Jan.g.meta.median2$Month_Day <- paste(kegg.ko.resid.Jan.g.meta.median2$Month, kegg.ko.resid.Jan.g.meta.median2$Day, sep = "/")


##Ridgeline plots
#August
Aug.resid.plot <- ggplot(data = subset(kegg.ko.resid.Aug.g.meta.median2, Year %in% c(2011:2020)), aes(y = Module.Description, x = Median_KO_Resid, fill = Module.Description)) +
  geom_density_ridges(alpha = 0.6, quantile_lines = TRUE, quantiles = c(0.05, 0.95), jittered_points = TRUE, point_size = 0.75, point_alpha = 1, point_shape = 21) + 
  geom_point(data = subset(kegg.ko.resid.Aug.g.meta.median2, Year %in% c(2021)), aes(), shape = 8,  position = position_nudge(y = 0.1), size = 2)+
  geom_text(data = subset(kegg.ko.resid.Aug.g.meta.median2, Year %in% c(2021)), aes(label = Month_Day), size = 2.5, vjust = 0.8) + 
  theme_classic() +
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  labs(y = "Kegg Module", x = "Median of Kegg Ortholog Residuals", title = "August") +
    scale_y_discrete(labels = c("Assimilatory sulfate reduction", "Catechol meta-cleavage", "Catechol ortho-cleavage", "Cumate degradation", "Dissimilatory sulfate reduction", "Nitrate assimilation", "Photosystem I", "Photosystem II", "Sulfate-sulfur assimilation", "Thiosulfate oxidation", "Toluene degradation")) +
  scale_fill_manual(values = c("deeppink2", "black", "grey30", "grey60", "deeppink4", "turquoise4", "darkseagreen3", "darkseagreen4", "goldenrod1", "darkorange3", "grey80")) + 
  xlim(-3, 4) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1))
Aug.resid.plot


#September
Sept.resid.plot <- ggplot(data = subset(kegg.ko.resid.Sept.g.meta.median2, Year %in% c(2011:2020)), aes(y = Module.Description, x = Median_KO_Resid, fill = Module.Description)) +
  geom_density_ridges(alpha = 0.6, quantile_lines = TRUE, quantiles = c(0.05, 0.95), jittered_points = TRUE, point_size = 0.75, point_alpha = 1, point_shape = 21) + 
  geom_point(data = subset(kegg.ko.resid.Sept.g.meta.median2, Year %in% c(2021)), aes(), shape = 8,  position = position_nudge(y = 0.1), size = 2)+
  geom_text(data = subset(kegg.ko.resid.Sept.g.meta.median2, Year %in% c(2021)), aes(label = Month_Day), size = 2.5, vjust = 0.8) + 
  theme_classic() + 
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  labs(y = "Kegg Module", x = "Median of Kegg Ortholog Residuals", title = "September") + 
  scale_y_discrete(labels = c("Assimilatory sulfate reduction", "Catechol meta-cleavage", "Catechol ortho-cleavage", "Cumate degradation", "Dissimilatory sulfate reduction", "Nitrate assimilation", "Photosystem I", "Photosystem II", "Sulfate-sulfur assimilation", "Thiosulfate oxidation", "Toluene degradation")) +
  scale_fill_manual(values = c("deeppink2", "black", "grey30", "grey60", "deeppink4", "turquoise4", "darkseagreen3", "darkseagreen4", "goldenrod1", "darkorange3", "grey80"))+ 
  xlim(-3, 4) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1))
Sept.resid.plot


#October
Oct.resid.plot <- ggplot(data = subset(kegg.ko.resid.Oct.g.meta.median2, Year %in% c(2011:2020)), aes(y = Module.Description, x = Median_KO_Resid, fill = Module.Description)) +
  geom_density_ridges(alpha = 0.6, quantile_lines = TRUE, quantiles = c(0.05, 0.95), jittered_points = TRUE, point_size = 0.75, point_alpha = 1, point_shape = 21) + 
  geom_point(data = subset(kegg.ko.resid.Oct.g.meta.median2, Year %in% c(2021)), aes(), shape = 8,  position = position_nudge(y = 0.1), size = 2)+
  geom_text(data = subset(kegg.ko.resid.Oct.g.meta.median2, Year %in% c(2021)), aes(label = Month_Day), size = 2.5, vjust = 0.8) + 
  theme_classic() + 
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  labs(y = "Kegg Module", x = "Median of Kegg Ortholog Residuals", title = "October") +
  scale_y_discrete(labels = c("Assimilatory sulfate reduction", "Catechol meta-cleavage", "Catechol ortho-cleavage", "Cumate degradation", "Dissimilatory sulfate reduction", "Nitrate assimilation", "Photosystem I", "Photosystem II", "Sulfate-sulfur assimilation", "Thiosulfate oxidation", "Toluene degradation")) +
  scale_fill_manual(values = c("deeppink2", "black", "grey30", "grey60", "deeppink4", "turquoise4", "darkseagreen3", "darkseagreen4", "goldenrod1", "darkorange3", "grey80"))+ 
  xlim(-3, 4) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1))
Oct.resid.plot


#November
Nov.resid.plot <- ggplot(data = subset(kegg.ko.resid.Nov.g.meta.median2, Year %in% c(2011:2020)), aes(y = Module.Description, x = Median_KO_Resid, fill = Module.Description)) +
  geom_density_ridges(alpha = 0.6, quantile_lines = TRUE, quantiles = c(0.05, 0.95), jittered_points = TRUE, point_size = 0.75, point_alpha = 1, point_shape = 21) + 
  geom_point(data = subset(kegg.ko.resid.Nov.g.meta.median2, Year %in% c(2021)), aes(), shape = 8,  position = position_nudge(y = 0.1), size = 2)+
  geom_text(data = subset(kegg.ko.resid.Nov.g.meta.median2, Year %in% c(2021)), aes(label = Month_Day), size = 2.5, vjust = 0.8) + 
  theme_classic() + 
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  labs(y = "Kegg Module", x = "Median of Kegg Ortholog Residuals", title = "November") +
  scale_y_discrete(labels = c("Assimilatory sulfate reduction", "Catechol meta-cleavage", "Catechol ortho-cleavage", "Cumate degradation", "Dissimilatory sulfate reduction", "Nitrate assimilation", "Photosystem I", "Photosystem II", "Sulfate-sulfur assimilation", "Thiosulfate oxidation", "Toluene degradation")) +
  scale_fill_manual(values = c("deeppink2", "black", "grey30", "grey60", "deeppink4", "turquoise4", "darkseagreen3", "darkseagreen4", "goldenrod1", "darkorange3", "grey80"))+ 
  xlim(-3, 4) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1))
Nov.resid.plot


#December
Dec.resid.plot <- ggplot(data = subset(kegg.ko.resid.Dec.g.meta.median2, Year %in% c(2011:2020)), aes(y = Module.Description, x = Median_KO_Resid, fill = Module.Description)) +
  geom_density_ridges(alpha = 0.6, quantile_lines = TRUE, quantiles = c(0.05, 0.95), jittered_points = TRUE, point_size = 0.75, point_alpha = 1, point_shape = 21) + 
  geom_point(data = subset(kegg.ko.resid.Dec.g.meta.median2, Year %in% c(2021)), aes(), shape = 8,  position = position_nudge(y = 0.1), size = 2)+
  geom_text(data = subset(kegg.ko.resid.Dec.g.meta.median2, Year %in% c(2021)), aes(label = Month_Day), size = 2.5, vjust = 0.8) + 
  theme_classic() + 
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  labs(y = "Kegg Module", x = "Median of Kegg Ortholog Residuals", title = "December") +
  scale_y_discrete(labels = c("Assimilatory sulfate reduction", "Catechol meta-cleavage", "Catechol ortho-cleavage", "Cumate degradation", "Dissimilatory sulfate reduction", "Nitrate assimilation", "Photosystem I", "Photosystem II", "Sulfate-sulfur assimilation", "Thiosulfate oxidation", "Toluene degradation")) +
  scale_fill_manual(values = c("deeppink2", "black", "grey30", "grey60", "deeppink4", "turquoise4", "darkseagreen3", "darkseagreen4", "goldenrod1", "darkorange3", "grey80"))+ 
  xlim(-3, 4) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1))
Dec.resid.plot


#January
Jan.resid.plot <- ggplot(data = subset(kegg.ko.resid.Jan.g.meta.median2, Year %in% c(2011:2021)), aes(y = Module.Description, x = Median_KO_Resid, fill = Module.Description)) +
  geom_density_ridges(alpha = 0.6, quantile_lines = TRUE, quantiles = c(0.05, 0.95), jittered_points = TRUE, point_size = 0.75, point_alpha = 1, point_shape = 21) + 
  geom_point(data = subset(kegg.ko.resid.Jan.g.meta.median2, Year %in% c(2022)), aes(), shape = 8,  position = position_nudge(y = 0.1), size = 2)+
  geom_text(data = subset(kegg.ko.resid.Jan.g.meta.median2, Year %in% c(2022)), aes(label = Month_Day), size = 2.5, vjust = 0.8) + 
  theme_classic() + 
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  labs(y = "Kegg Module", x = "Median of Kegg Ortholog Residuals", title = "January") +
  scale_y_discrete(labels = c("Assimilatory sulfate reduction", "Catechol meta-cleavage", "Catechol ortho-cleavage", "Cumate degradation", "Dissimilatory sulfate reduction", "Nitrate assimilation", "Photosystem I", "Photosystem II", "Sulfate-sulfur assimilation", "Thiosulfate oxidation", "Toluene degradation")) +
  scale_fill_manual(values = c("deeppink2", "black", "grey30", "grey60", "deeppink4", "turquoise4", "darkseagreen3", "darkseagreen4", "goldenrod1", "darkorange3", "grey80"))+ 
  xlim(-3, 4) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1))
Jan.resid.plot

####KEGG Module correlations with other Modules####
#Photosystem I
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem I`, x = SCOP.kegg.zscore.cor2$`Toluene degradation`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem I`, x = SCOP.kegg.zscore.cor2$`Cumate degradation`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem I`, x = SCOP.kegg.zscore.cor2$`Catechol meta-cleavage`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem I`, x = SCOP.kegg.zscore.cor2$`Catechol ortho-cleavage`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem I`, x = SCOP.kegg.zscore.cor2$`Sulfate-sulfur assimilation`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem I`, x = SCOP.kegg.zscore.cor2$`Thiosulfate oxidation`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem I`, x = SCOP.kegg.zscore.cor2$`Assimilatory sulfate reduction`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem I`, x = SCOP.kegg.zscore.cor2$`Dissimilatory sulfate reduction`, method = "pearson", conf.level = 0.95) #p-value = 0.03926; cor = 0.372122 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem I`, x = SCOP.kegg.zscore.cor2$`Homoprotocatechuate degradation`, method = "pearson", conf.level = 0.95) 


#Photosystem II
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem II`, x = SCOP.kegg.zscore.cor2$`Toluene degradation`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem II`, x = SCOP.kegg.zscore.cor2$`Cumate degradation`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem II`, x = SCOP.kegg.zscore.cor2$`Catechol meta-cleavage`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem II`, x = SCOP.kegg.zscore.cor2$`Catechol ortho-cleavage`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem II`, x = SCOP.kegg.zscore.cor2$`Sulfate-sulfur assimilation`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem II`, x = SCOP.kegg.zscore.cor2$`Thiosulfate oxidation`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem II`, x = SCOP.kegg.zscore.cor2$`Assimilatory sulfate reduction`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem II`, x = SCOP.kegg.zscore.cor2$`Dissimilatory sulfate reduction`, method = "pearson", conf.level = 0.95) #p-value = 0.01693; cor = 0.4257975 
cor.test(y = SCOP.kegg.zscore.cor2$`Photosystem II`, x = SCOP.kegg.zscore.cor2$`Homoprotocatechuate degradation`, method = "pearson", conf.level = 0.95) 


#Thiosulfate oxidation
cor.test(y = SCOP.kegg.zscore.cor2$`Thiosulfate oxidation`, x = SCOP.kegg.zscore.cor2$`Toluene degradation`, method = "pearson", conf.level = 0.95) #p-value = 0.001751; 0.5391152
cor.test(y = SCOP.kegg.zscore.cor2$`Thiosulfate oxidation`, x = SCOP.kegg.zscore.cor2$`Cumate degradation`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Thiosulfate oxidation`, x = SCOP.kegg.zscore.cor2$`Catechol meta-cleavage`, method = "pearson", conf.level = 0.95) #p-value = 4.892e-07; cor = 0.7667898 
cor.test(y = SCOP.kegg.zscore.cor2$`Thiosulfate oxidation`, x = SCOP.kegg.zscore.cor2$`Catechol ortho-cleavage`, method = "pearson", conf.level = 0.95) #p-value = 3.512e-10; cor = 0.8648735
cor.test(y = SCOP.kegg.zscore.cor2$`Thiosulfate oxidation`, x = SCOP.kegg.zscore.cor2$`Sulfate-sulfur assimilation`, method = "pearson", conf.level = 0.95) #p-value = 4.253e-06; cor = 0.7234689 
cor.test(y = SCOP.kegg.zscore.cor2$`Thiosulfate oxidation`, x = SCOP.kegg.zscore.cor2$`Assimilatory sulfate reduction`, method = "pearson", conf.level = 0.95) #p-value = 1.615e-10; cor = 0.8723635
cor.test(y = SCOP.kegg.zscore.cor2$`Thiosulfate oxidation`, x = SCOP.kegg.zscore.cor2$`Dissimilatory sulfate reduction`, method = "pearson", conf.level = 0.95) #0.009941; cor = 0.4559541
cor.test(y = SCOP.kegg.zscore.cor2$`Thiosulfate oxidation`, x = SCOP.kegg.zscore.cor2$`Homoprotocatechuate degradation`, method = "pearson", conf.level = 0.95) #p-value = 4.314e-07; cor = 0.7690588


#Sulfate-sulfur assimilation
cor.test(y = SCOP.kegg.zscore.cor2$`Sulfate-sulfur assimilation`, x = SCOP.kegg.zscore.cor2$`Toluene degradation`, method = "pearson", conf.level = 0.95) #p-value = 0.02371; 0.4052862 
cor.test(y = SCOP.kegg.zscore.cor2$`Sulfate-sulfur assimilation`, x = SCOP.kegg.zscore.cor2$`Cumate degradation`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Sulfate-sulfur assimilation`, x = SCOP.kegg.zscore.cor2$`Catechol meta-cleavage`, method = "pearson", conf.level = 0.95) #p-value = 4.892e-07; cor = 0.5491481 
cor.test(y = SCOP.kegg.zscore.cor2$`Sulfate-sulfur assimilation`, x = SCOP.kegg.zscore.cor2$`Catechol ortho-cleavage`, method = "pearson", conf.level = 0.95) #p-value = 0.0007547; cor = 0.57297 
cor.test(y = SCOP.kegg.zscore.cor2$`Sulfate-sulfur assimilation`, x = SCOP.kegg.zscore.cor2$`Assimilatory sulfate reduction`, method = "pearson", conf.level = 0.95) #p-value = 1.146e-07; cor = 0.7915076
cor.test(y = SCOP.kegg.zscore.cor2$`Sulfate-sulfur assimilation`, x = SCOP.kegg.zscore.cor2$`Dissimilatory sulfate reduction`, method = "pearson", conf.level = 0.95) #p-value = 1.482e-05; cor = 0.6941731
cor.test(y = SCOP.kegg.zscore.cor2$`Sulfate-sulfur assimilation`, x = SCOP.kegg.zscore.cor2$`Homoprotocatechuate degradation`, method = "pearson", conf.level = 0.95) #p-value = 0.0006592; cor = 0.5780902 


#Assimilatory sulfate reduction
cor.test(y = SCOP.kegg.zscore.cor2$`Assimilatory sulfate reduction`, x = SCOP.kegg.zscore.cor2$`Toluene degradation`, method = "pearson", conf.level = 0.95) #p-value = 0.0005232; cor = 0.5866405 
cor.test(y = SCOP.kegg.zscore.cor2$`Assimilatory sulfate reduction`, x = SCOP.kegg.zscore.cor2$`Cumate degradation`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Assimilatory sulfate reduction`, x = SCOP.kegg.zscore.cor2$`Catechol meta-cleavage`, method = "pearson", conf.level = 0.95) #p-value = 3.147e-10; cor = 0.8659595 
cor.test(y = SCOP.kegg.zscore.cor2$`Assimilatory sulfate reduction`, x = SCOP.kegg.zscore.cor2$`Catechol ortho-cleavage`, method = "pearson", conf.level = 0.95) #p-value = 8.98e-11; cor = 0.8777321
cor.test(y = SCOP.kegg.zscore.cor2$`Assimilatory sulfate reduction`, x = SCOP.kegg.zscore.cor2$`Dissimilatory sulfate reduction`, method = "pearson", conf.level = 0.95) #p-value = 2.832e-06; cor = 0.7322906
cor.test(y = SCOP.kegg.zscore.cor2$`Assimilatory sulfate reduction`, x = SCOP.kegg.zscore.cor2$`Homoprotocatechuate degradation`, method = "pearson", conf.level = 0.95) #p-value = 0.0001508; cor = 0.6289539 


#Dissimilatory sulfate reduction
cor.test(y = SCOP.kegg.zscore.cor2$`Dissimilatory sulfate reduction`, x = SCOP.kegg.zscore.cor2$`Toluene degradation`, method = "pearson", conf.level = 0.95) #p-value = 0.02353; cor = 0.4057432 
cor.test(y = SCOP.kegg.zscore.cor2$`Dissimilatory sulfate reduction`, x = SCOP.kegg.zscore.cor2$`Cumate degradation`, method = "pearson", conf.level = 0.95) 
cor.test(y = SCOP.kegg.zscore.cor2$`Dissimilatory sulfate reduction`, x = SCOP.kegg.zscore.cor2$`Catechol meta-cleavage`, method = "pearson", conf.level = 0.95) #p-value = 1.027e-05; cor = 0.7031396 
cor.test(y = SCOP.kegg.zscore.cor2$`Dissimilatory sulfate reduction`, x = SCOP.kegg.zscore.cor2$`Catechol ortho-cleavage`, method = "pearson", conf.level = 0.95) #p-value = 0.0003463; cor = 0.6013528
cor.test(y = SCOP.kegg.zscore.cor2$`Dissimilatory sulfate reduction`, x = SCOP.kegg.zscore.cor2$`Homoprotocatechuate degradation`, method = "pearson", conf.level = 0.95)


####Plot export####
library(ggpubr)

#functional diversity
png("SCOP_Bacterial_FunctionalDiversity.png", units = "in", width = 20, height = 13, res = 600)
ggarrange(bacteria.ko.shannon.plot, bacteria.ko.richness.plot, bacteria.ko.evenness.plot, bacteria.ko.pca, bacteria.ko.pc1, bacteria.ko.eucdiss.plot, nrow = 2, ncol = 3, common.legend = T, legend = "bottom", labels = c("a", "b", "c", "d", "e", "f"), font.label = list(size = 24, color = "black", face = "bold", family = "sans"), align = "hv")
dev.off()


#ridgeline plots
png("SCOP_KOModules_Residuals_Density.png", units = "in", width = 15, height = 8, res = 600)
ggarrange(Aug.resid.plot, Sept.resid.plot, Oct.resid.plot, Nov.resid.plot, Dec.resid.plot, Jan.resid.plot, nrow = 2, ncol = 3, align = "hv", labels = c("a", "b", "c", "d", "e", "f"), font.label = list(size = 16, color = "black", face = "bold", family = "sans"))
dev.off()

svg("SCOP_KOModules_Residuals_Density.svg", width = 15, height = 8)
ggarrange(Aug.resid.plot, Sept.resid.plot, Oct.resid.plot, Nov.resid.plot, Dec.resid.plot, Jan.resid.plot, nrow = 2, ncol = 3, align = "hv", labels = c("a", "b", "c", "d", "e", "f"), font.label = list(size = 16, color = "black", face = "bold", family = "sans"))
dev.off()
