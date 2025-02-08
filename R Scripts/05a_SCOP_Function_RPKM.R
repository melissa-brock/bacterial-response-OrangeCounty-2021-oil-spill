###Functional analysis of SCOP short-read metagenomes processed through the NMDC-Edge Platform
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

####Calculate RPKM of genes for each functional assignment system####
#Unload packages; this MUST be done to use the map and group_by functions correctly
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

#Load only necessary packages
library(vegan)
library(tidyr)
library(purrr)
library(dplyr)

#Group by Sample and unique function assignment; calculate sum of RPKM for each unique function; divide by total RPKM; multiply by 100 to get relative abundance of each function
test.ko <- map(sample_list2, ~.x %>% group_by(SampleID, ko) %>% summarize(KO_SUM=sum(RPKM)) %>% mutate(KO_Percent=(KO_SUM/sum(KO_SUM))*100))

####Combine RPKM dataframes for each functional assignment system####
ko.combined <- Reduce(full_join, test.ko)

####Subset based on number of assignments to each gene####
ko.sep.mult <- ko.combined[grep(",", ko.combined$ko),]
ko.sep.single <- anti_join(ko.combined, ko.sep.mult, by = "ko")

####Calculate relative abundance of multiple assignments####
sum(ko.combined$KO_Percent)
sum(ko.sep.single$KO_Percent)
sum(ko.sep.mult$KO_Percent)

####Export dataframes of single annotation RPKMs####
write.csv(ko.sep.single, "SCOP_Bacteria_KO_RPKM.csv")

####Convert to wide format and format dataframe####
ko.sep.single.sub <- ko.sep.single[, -c(4)]
ko.sep.single.wide <- as.data.frame(spread(ko.sep.single.sub, SampleID, KO_SUM, fill = 0))
ko.sep.single.wide$ko[is.na(ko.sep.single.wide$ko)] <- "Unassigned"
ko.sep.single.wide.sub <- ko.sep.single.wide[, -1]
row.names(ko.sep.single.wide.sub) <- ko.sep.single.wide$ko
ko.sep.single.wide.sub2 <- ko.sep.single.wide.sub[which(row.names(ko.sep.single.wide.sub) != "Unassigned"), ]

####Calculate Shannon Index, richness, and evenness####
ko.shannon <- as.data.frame(diversity(t(ko.sep.single.wide.sub2), index = "shannon"))
ko.evenness <- as.data.frame(diversity(t(ko.sep.single.wide.sub2), index = "shannon")/log(specnumber(t(ko.sep.single.wide.sub2))))
ko.richness <- as.data.frame(specnumber(t(ko.sep.single.wide.sub2)))

####Export dataframes of diversity metrics####
write.csv(ko.shannon, "SCOP_Bacteria_KO_Shannon.csv")
write.csv(ko.evenness, "SCOP_Bacteria_KO_Evenness.csv")
write.csv(ko.richness, "SCOP_Bacteria_KO_Richness.csv")

####Close environment####
save.image("SCOP_BacteriaFunction_RPKM.RData")
quit()