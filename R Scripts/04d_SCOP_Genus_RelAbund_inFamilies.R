###Script to calculate relative abundances of genera within each family's residuals that are plotted in the ridgeline plots - Melissa Brock
####Setup environment####
getwd()
setwd("/annotation-summaries/")
getwd()

library(tidyverse)
library(phyloseq)
library(plyr)

####File import & formatting####
#Read in filter statistics metadata sheet
NP.meta <- read.csv("NPfilterstats.csv", header = T, sep = ",")

#Define function that adds unique filenames as a column
read_csv_filename <- function(filename){
  ret <- read.csv(filename)
  ret$Filename <- filename #Edit
  ret
}

#Make list of files that contain a specific pattern
filenames <- list.files(path = ".", pattern = "SCOP_2", all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)

#Import list of files 
import.list <- ldply(filenames, read_csv_filename)

#Split into multiple objects based on filename
sample_list <- group_split(import.list, Filename)

#Subset to Bacteria only
sample_list2 <- lapply(sample_list, function(x) x[which(x$Domain_ConsensusID == "Bacteria"),])

#Remove duplicate contig values within each sample
#DO NOT DO THIS WHEN ANALYZING FUNCTION; ONLY USE FOR ANALYSIS OF TAXA
sample_list2 <- lapply(sample_list, function(x) x[!duplicated(x[c("contig")], fromLast = F), ])

#Remove "_annotation.csv" from end of "Filename" to create new column called "SampleID"
sample_list2 <- map(sample_list2, ~.x[]%>% mutate(SampleID = str_replace_all(Filename,"_annotation.csv", "")))

#Merge with filter statistics
merge.func <- function(x,y){merge(x, y, by="SampleID")}
sample_list2 <- lapply(sample_list2, merge.func, NP.meta)

####Normalize read counts using RPKM (reads per kilobase million)####
#RPKM = numReads / (geneLength/1000 * totalNumReads/1,000,000)
#numReads        - number of reads mapped to a gene sequence
#geneLength      - length of the gene sequence
#totalNumReads   - total number of mapped reads of a sample (outputReads)
#AvgFold = numReads/geneLength
#Calculation becomes: RPKM = Avg_fold / (totalNumReads/1000*1,000,000)

#Calculate RPKM
rpkm.func <- function(x){x$RPKM <- x$Avg_fold/(x$outputReads/(1000*1000000));x}
sample_list2 <- lapply(sample_list2, rpkm.func)

#housekeeping point to check that calculation is working properly when performed on a list of dataframes
range(sample_list2[[1]]$RPKM)
range(sample_list2[[2]]$RPKM) 

####Subset to taxa of interest####
taxa_families <- c("Synechococcaceae", "Prochloraceae", "Rhodobacteraceae", "Pelagibacteraceae", "Halomonadaceae", "Flavobacteriaceae", "Cryomorphaceae")

sample_list3 <- lapply(sample_list2, function(x) x[which(x$Family_ConsensusID %in% taxa_families),])

####Calculate relative abundances of genera within each family using RPKM values####
#Unload packages; this MUST be done to use the map and group_by functions correctly
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

#Load only necessary packages
library(vegan)
library(tidyr)
library(purrr)
library(dplyr)

test.genus <- map(sample_list3, ~.x %>% group_by(SampleID, Family_ConsensusID, Genus_ConsensusID) %>% summarize(Genus_SUM=sum(RPKM)) %>% mutate(Genus_Percent=(Genus_SUM/sum(Genus_SUM))*100))

####Combine relative abundance dataframes####
genus.combined <- Reduce(full_join, test.genus)

####Split into dataframes based on Family Consensus ID####
syn.df <- subset(genus.combined, Family_ConsensusID == 'Synechococcaceae')
pro.df <- subset(genus.combined, Family_ConsensusID == 'Prochloraceae')
rhodo.df <- subset(genus.combined, Family_ConsensusID == 'Rhodobacteraceae')
pelagi.df <- subset(genus.combined, Family_ConsensusID == 'Pelagibacteraceae')
halo.df <- subset(genus.combined, Family_ConsensusID == 'Halomonadaceae')
flavo.df <- subset(genus.combined, Family_ConsensusID == 'Flavobacteriaceae')
cryo.df <- subset(genus.combined, Family_ConsensusID == 'Cryomorphaceae')

####Format subsetted dataframes####
syn.spread <- syn.df %>% ungroup() %>% select(-c(Family_ConsensusID, Genus_SUM)) %>% spread(key = Genus_ConsensusID, value = Genus_Percent)

pro.spread <- pro.df %>% ungroup() %>% select(-c(Family_ConsensusID, Genus_SUM)) %>% spread(key = Genus_ConsensusID, value = Genus_Percent)

rhodo.spread <- rhodo.df %>% ungroup() %>% select(-c(Family_ConsensusID, Genus_SUM)) %>% spread(key = Genus_ConsensusID, value = Genus_Percent)

pelagi.spread <- pelagi.df %>% ungroup() %>% select(-c(Family_ConsensusID, Genus_SUM)) %>% spread(key = Genus_ConsensusID, value = Genus_Percent)

halo.spread <- halo.df %>% ungroup() %>% select(-c(Family_ConsensusID, Genus_SUM)) %>% spread(key = Genus_ConsensusID, value = Genus_Percent)

flavo.spread <- flavo.df %>% ungroup() %>% select(-c(Family_ConsensusID, Genus_SUM)) %>% spread(key = Genus_ConsensusID, value = Genus_Percent)

cryo.spread <- cryo.df %>% ungroup() %>% select(-c(Family_ConsensusID, Genus_SUM)) %>% spread(key = Genus_ConsensusID, value = Genus_Percent)

####Export dataframes####
write.csv(syn.spread, "SCOP_Genus_RelAbund_inSyn.csv")
write.csv(pro.spread, "SCOP_Genus_RelAbund_inPro.csv")
write.csv(rhodo.spread, "SCOP_Genus_RelAbund_inRhodo.csv")
write.csv(pelagi.spread, "SCOP_Genus_RelAbund_inPelagi.csv")
write.csv(halo.spread, "SCOP_Genus_RelAbund_inHalo.csv")
write.csv(flavo.spread, "SCOP_Genus_RelAbund_inFlavo.csv")
write.csv(cryo.spread, "SCOP_Genus_RelAbund_inCryo.csv")

####Save environment####
save.image("SCOP_Genus_RelAbund_inFamilies.RData")
quit()