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
"04b_SCOP_Genus_Diversity.R" contains R code that calculates taxa diversity metrics (richness, evenness, and Shannon), calculates changes in community structure, and visualizes these trends.  
"04c_SCOP_Taxa_RelAbund_AnalysisandPlots.R" contains R code that plots changes in relative abundance of the top 25 most abundant genera, that plots changes in taxa found to be enriched during deepwater horizon, and that plots the residuals of family abundances across the ten-year time-series.  
"04d_SCOP_Genus_RelAbund_inFamilies.R" contains R code that calculates the relative abundance of genera assigned to each family that was plotted in the ridgeline plots.  

"05a_SCOP_Function_RPKM.R" contains R code that calculates RPKM of each Kegg Ortholog.  
"05b_SCOP_FunctionalDiversity_PCAs.R" contains R code that calculates functional diversity metrics (richness, evenness, and Shannon) and performs ordination.   
"05c_SCOP_Function_AnalysisandPlots.R" contains R code for plotting functional diversity metrics, constructing PCA, showing changes in time in Kegg Orthologs and Kegg Moduels, and plotting the residuals of KEGG Modules as ridgeline plots.  
"05d_SCOP_KOs_in_Taxa.R" contains R code that identifies what taxa the indicator functions are found in.    
"05e_SCOP_KOs_in_KEGGModules.R" contains R code that matches KOs with the Kegg Module that they are found in.  

"06_SCOP_EnvironmentalData_and_PAHs.R" contains R code for plotting environmental conditions and polycyclic aromatic hydrocarbon (PAH) concentrations.  
"07_SCOP_FlowCytometry_and_Microscopy.R" contains R code for plotting flow cytometry and microscopy data.  
