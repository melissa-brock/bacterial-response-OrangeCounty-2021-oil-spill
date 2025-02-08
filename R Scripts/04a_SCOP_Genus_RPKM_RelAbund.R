##
###Taxonomy analysis of short-read metagenomes processed through the NMDC-Edge Platform
##
####Setup environment####
getwd()
setwd("/NMDC_Output/annotation-summaries/")
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

####Calculate relative abundances of taxa at each taxonomic level using RPKM values####
#Unload packages; this MUST be done to use the map and group_by functions correctly
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

#Load only necessary packages
library(vegan)
library(tidyr)
library(purrr)
library(dplyr)

#Group by Sample and unique taxa assignment; calculate sum of RPKM for each unique taxa; divide by total RPKM; multiply by 100 to get relative abundance of each taxa
genus.RPKM <- map(sample_list2, ~.x %>% group_by(SampleID, Domain_ConsensusID, Genus_ConsensusID) %>% summarize(Genus_SUM=sum(RPKM)))

genus.RelAbund <- map(sample_list2, ~.x %>% group_by(SampleID, Domain_ConsensusID, Genus_ConsensusID) %>% summarize(Genus_SUM=sum(RPKM)) %>% mutate(Genus_Percent=(Genus_SUM/sum(Genus_SUM))*100))

####Combine dataframes####
genus.RPKM.combined <- Reduce(full_join, genus.RPKM)
genus.RelAbund.combined <- Reduce(full_join, genus.RelAbund)

####Export dataframes####
write.csv(genus.RPKM.combined, "Genus_RPKM.csv")
write.csv(genus.RelAbund.combined, "Genus_RelAbund.csv")

####Close environment####
save.image("Genus_RPKM_RelAbund.RData")
quit()
