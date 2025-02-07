General Overview:  
1. Run python scripts on output files from NMDC-Edge.
2. Run R scripts on output files from python scripts. 

#Python Scripts  
"01_scop-consensustaxonomy.py" contains python code that calculates consensus taxonomy using a majority rule.  
"02_scop-annotate-combo.py" contains python code combining contig coverage information, functional annotations, and consensus taxonomy into a single file for each sample.  

#R Scripts  
"03_SCOP_Genus_RPKM_RelAbund.R" contains R code that calculates RPKM and relative abundance of each taxa at the genus level.  
