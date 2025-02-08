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

####Subset to list of genes####
gene.list <- read.csv("KEGG_Modules_KO_List.csv", header = T, sep = ",")
sample_list3 <- lapply(sample_list2, function(x) x[x$ko %in% gene.list$KO.ID,])

####Combine RPKM dataframes for subset of functional genes####
sample_list_combined <- Reduce(full_join, sample_list3)

####Export dataframes of functional genes####
write.csv(sample_list_combined, "SCOP_Bacteria_KO_Subset_RPKM.csv")

####Close environment####
save.image("SCOP_BacteriaFunction_Comp.RData")
quit()
