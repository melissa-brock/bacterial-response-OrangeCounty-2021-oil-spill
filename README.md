General Overview:  
1. Run python scripts on output files from NMDC-Edge.
2. Run MatLab scripts on frequency tables to calculate residuals.
3. Run R scripts on output files from python and MatLab scripts. 

#Python Scripts  
"01_scop-consensustaxonomy.py" contains python code that calculates consensus taxonomy using a majority rule.  
"02_scop-annotate-combo.py" contains python code combining contig coverage information, functional annotations, and consensus taxonomy into a single file for each sample.  

#MatLab Scripts  
"03a_EachAnalysisTable_wResidual.m" contains MatLab code that calculates residuals for Kegg Orthologs and for taxa at the family level.   
"03b_create_analysis_tables_wResiduals.m" contains MatLab code that combines and exports the frequency tables with the residuals.    

#R Scripts  
"04a_SCOP_Genus_RPKM_RelAbund.R" contains R code that calculates RPKM and relative abundance of each taxa at the genus level.  
