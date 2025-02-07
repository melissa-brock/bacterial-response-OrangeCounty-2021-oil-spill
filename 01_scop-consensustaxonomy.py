#!/usr/bin/env python

## Import packages
import pandas as pd
import numpy as np
from functools import reduce

import string
import sys
import re
import os
import glob


## Define functions 

def check_taxon(taxon, df):
	#Calculate the percentage of reads from each contig that are assigned to each Strain
	df_taxa = df.loc[:, ["ID", "ContigID", "TaxaID", taxon]]
	df_taxa['gene_count'] = df_taxa.groupby('ContigID')['ContigID'].transform('size').astype(int)
	df_taxa['taxon_count'] = df_taxa.groupby(['ContigID', taxon])[taxon].transform('size').astype(int)
	df_taxa['taxon_perc'] = ((df_taxa['taxon_count']/df_taxa['gene_count'])*100).astype(float)
	df_taxa2 = df_taxa.drop_duplicates(subset=['ContigID', taxon])
	df_taxa3 = df_taxa2.loc[:]
	#Assign taxonomy: >= 50% of reads in a single contig match to a specific taxonomic assignment
	df_taxa3[str(taxon)] = df_taxa3.loc[df_taxa3['taxon_perc'] >= 50, taxon]
	#Pull out assignment for each contig
	df_taxa4 = df_taxa3[['ContigID','TaxaID',str(taxon),'taxon_perc']]
	df_taxa5 = df_taxa4.loc[:]
	#Drop rows with duplicate taxa IDs per contig
	df_taxa5 = df_taxa5.loc[df_taxa5.astype(str).drop_duplicates(subset=['ContigID',str(taxon)]).index]
	#Drop rows with NA taxa IDs
	df_taxa5 = df_taxa5[~df_taxa5['ContigID'].duplicated(keep=False) | df_taxa5[[str(taxon)]].notnull().any(axis=1)]
	df_taxa5.columns = ['ContigID',str(taxon)+'_TaxaID',taxon,str(taxon)+'_perc']
	return(df_taxa5)

def consensus_tax(row):
	if pd.isnull(row['Level']): 
		return row['Strain_TaxaID']
	elif row['Level'] == 'Domain': 
		return row['Domain_TaxaID']
	else:
		return row[str(taxonomies[taxonomies.index(str(row['Level']))+1]+'_TaxaID')]

def level_repl(row):
	if not pd.isnull(row['Level']):
			idx = row.index.get_loc(str(row['Level']))
			row.iloc[idx:len(row)] = None 
			row.fillna(value=np.nan, inplace=True)
			return


## Read in file
path = r'/NMDC_Output/SAMPLENAME/'

phylo_file = glob.glob(os.path.join(path, "*_gene_phylogeny.tsv"))

gene_phylo = pd.read_table(phylo_file[0], delimiter='\t',header=None)

## Name columns in gene phylogeny file
gene_phylo.columns = ["ID", "Unk1", "Unk2", "Percent_Match", "TaxaID"]

## Split TaxaID into separate columns
TaxaID = gene_phylo['TaxaID']
gene_phylo[["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]] = TaxaID.str.split(";", n = 7, expand=True)

## Make new column with contig ID in it
gene_phylo["Contig"].astype(str)
gene_phylo["ContigID"] = gene_phylo["ID"].str.split("_",6).str[0:6].str.join('_')

## Define taxonomic levels 
taxonomies = ['Strain', 'Species', 'Genus','Family','Order','Class','Phylum','Domain']
rev_taxonomies = taxonomies[::-1]

## Create empty taxonomy dataframe 
taxonomy_check = pd.DataFrame(columns = ["ContigID"])

for taxon in taxonomies:
	result = check_taxon(taxon, gene_phylo)
	taxonomy_dfs = [taxonomy_check, result]
	taxonomy_check = reduce(lambda  left,right: pd.merge(left,right,on=['ContigID'], how='outer'), taxonomy_dfs)

## Keep all unique assignments 
taxonomy_unique = taxonomy_check.drop_duplicates(subset=['ContigID'],keep=False)

## Create separate dataframe with duplicated assignments 
taxonomy_dup = taxonomy_check.loc[taxonomy_check.duplicated(subset='ContigID',keep=False)]

## Remove rows with unresolved taxonomy at the Domain level
taxonomy_dup = taxonomy_dup.loc[taxonomy_dup['Domain_perc'] > 50]

## Set duplicated taxonomies <= 50% to NaN 
for taxon in taxonomies:
	taxonomy_dup.loc[taxonomy_dup[str(taxon)+'_perc'] <= 50,str(taxon)] = np.nan 

## Create df of now unduplicated taxonomies 
taxonomy_undup = taxonomy_dup.drop_duplicates(subset=['ContigID']+taxonomies,keep='first')

## Concatenate unique dfs 
taxonomy_concat = pd.concat([taxonomy_unique,taxonomy_undup])

## Find highest level NaN 
fintax_rev = taxonomy_concat[taxonomy_concat.columns[::-1]]
bool_tax = fintax_rev.isnull()
bool_test = bool_tax[bool_tax.any(axis=1)].idxmax(axis=1)

fintax_level = pd.concat([taxonomy_concat,bool_test],axis=1)
fintax_level.columns = [*fintax_level.columns[:-1], 'Level']

## Find taxon assignment one level up from highest NaN 
pd.options.mode.chained_assignment = None

fintax_level['Result'] = fintax_level.apply(consensus_tax,axis=1)

## Reformat data 
fintax_subset = fintax_level[['ContigID','Level','Result']]
taxonomy_select = fintax_subset['Result']
fintax_subset[["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]] = taxonomy_select.str.split(";", n = 7, expand=True)

## Replace all row values (taxa assignments) at highest level NaN and below with NaN 
for index,row in fintax_subset.iterrows():
	if not pd.isnull(row['Level']):
			idx = row.index.get_loc(str(row['Level']))
			row.iloc[idx:len(row)] = None 
			row.fillna(value=np.nan, inplace=True)


## Output taxonomy dataframe
final_tax = fintax_subset[['ContigID']+rev_taxonomies]
final_taxonomy = final_tax.sort_values(by=['ContigID'])

pd.DataFrame.to_csv(final_taxonomy, os.path.join(path,"SAMPLENAME_consensus_taxonomy.csv"), sep=',', na_rep='NA', index=False)


