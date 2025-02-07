#!/usr/bin/python

import string
import sys
import re
import os
import glob
import pandas as pd

path = r'/NMDC_Output/SAMPLENAME/'

## Import contig coverage
fullcoverage_df = pd.read_table(os.path.join(path, "covstats.tsv"), delimiter='\t', header=1)
fullcoverage_df.columns = ['contig','Avg_fold','Length','Ref_GC','Covered_percent','Covered_bases','Plus_reads','Minus_reads','Read_GC','Median_fold','Std_Dev']

coverage_df = fullcoverage_df[['contig','Avg_fold','Median_fold','Std_Dev','Length','Covered_percent','Read_GC']]


## Import functional annotation 
anno_file = glob.glob(os.path.join(path, "*_functional_annotation.gff"))

df = pd.read_table(anno_file[0], delimiter='\t',header=None)

df.columns = ['contig','source','feature','start','end','score','strand','frame','attribute']

def find_match(string_list, wanted):
    for string in string_list:
        if string.startswith(wanted):
            return string.replace(wanted,'')
    return None
    
matches = ['ID=','product=','product_source=','cath_funfam=','cog=','ko=','ec_number=','pfam=','smart=','superfamily=','tigrfam=']

new_rows = [] 

for index, row in df.iterrows():
	row_str = row['attribute']
	row_list = row_str.split(';')
	search_result = []
	for match in matches:
		result = find_match(row_list,match)
		search_result.append(result)
	new_rows.append(search_result)

colnames = ['gene','annotation','annotation_source','cath_funfam','cog','ko','ec_number','pfam','smart','superfamily','tigrfam']
	
new_df = pd.DataFrame(new_rows,columns=colnames)

annotate_df = pd.concat([df['contig'],new_df],axis=1)

## Import consensus taxonomy 

taxa_file = glob.glob(os.path.join(path, "*_consensus_taxonomy.csv"))

taxa_df = pd.read_csv(taxa_file[0],header=0)
taxa_df.columns = ['contig','Domain_ConsensusID','Phylum_ConsensusID','Class_ConsensusID','Order_ConsensusID','Family_ConsensusID','Genus_ConsensusID','Species_ConsensusID','Strain_ConsensusID']

## Final merged annotation 

mergeanno_df = pd.merge(annotate_df,coverage_df,on='contig')

final_df = pd.merge(mergeanno_df,taxa_df,on='contig')

final_df.to_csv("SAMPLENAME_annotation.csv",na_rep='NA',index=False)


