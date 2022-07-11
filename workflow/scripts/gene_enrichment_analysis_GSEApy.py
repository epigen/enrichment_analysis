#!/bin/env python

# load libraries
import pandas as pd
import json
# import pickle
import os
import numpy as np
import gseapy as gp


# utils for manual odds ratio calculation
def overlap_converter(overlap_str, bg_n, gene_list_n):
    overlap_n, gene_set_n = str(overlap_str).split('/')
    return odds_ratio_calc(bg_n, gene_list_n, int(gene_set_n), int(overlap_n))

def odds_ratio_calc(bg_n, gene_list_n, gene_set_n, overlap_n): 
    import scipy.stats as stats
    table=np.array([[gene_set_n, bg_n-gene_set_n],[overlap_n, gene_list_n-overlap_n]])
    oddsratio, pvalue = stats.fisher_exact(table)
    return (1/oddsratio)


# configs
query_genes_path = snakemake.input['query_genes']
background_genes_path = snakemake.input['background_genes']
database_path = snakemake.input['enrichr_database']
db = snakemake.params["database"]
result_path = snakemake.output['result_file']
dir_results = os.path.dirname(result_path)

if not os.path.exists(dir_results):
    os.mkdir(dir_results)

# check if genes file exists & load or handle exception
if os.path.exists(query_genes_path):
    genes = open(query_genes_path, "r")
    gene_list = genes.read()
    gene_list = gene_list.split('\n')
    genes.close()
else:
    with open(os.path.join(dir_results,"no_genes_found.txt"), 'w') as f:
        f.write('no genes found')
    quit()

# load background genes
bg_file = open(background_genes_path, "r")
background = bg_file.read()
background = background.split('\n')
bg_file.close()

# load database .pkl file
# with open(enrichr_databases, 'rb') as f:
#     db_dict = pickle.load(f)
with open(database_path) as json_file:
    db_dict = json.load(json_file)
    
# convert gene lists to upper case
gene_list=[str(x).upper() for x in list(gene_list)]
background=[str(x).upper() for x in list(background)]
    
# perform enrichment of every database with GSEApy (plots are generated automatically)
bg_n = len(background)

res = dict()

#for db in db_dict.keys():
res = gp.enrichr(gene_list=gene_list,
                 gene_sets=db_dict,#db_dict[db],
                 background=background,
#                      organism='mouse',
                 outdir=os.path.join(dir_results),
                 top_term=25,
                 cutoff=0.05,
                 format='png',
                 verbose=False,
                )

# # move on if result is empty
# if res.results.shape[0]!=0:
#     continue

# annotate used gene set
res.results['Gene_set'] = db

# odds ratio calculation 
gene_list_n=len(gene_list)
res.results['Odds Ratio'] = res.results['Overlap'].apply(overlap_converter, args=(bg_n, gene_list_n))

# make column names language agnostic (ie R compatible)
column_names = list(res.results.columns.values)
column_names = [col.replace(" ","_") for col in column_names]
column_names = [col.replace("-","_") for col in column_names]
res.results.columns = column_names

# separate export
res.results.to_csv(result_path)

