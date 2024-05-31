#!/bin/env python

# load libraries
import pandas as pd
import json
# import pickle
import os
import numpy as np
import gseapy as gp
import sys

# # utils for manual odds ratio calculation -> not used anymore
# def overlap_converter(overlap_str, bg_n, gene_list_n):
#     overlap_n, gene_set_n = str(overlap_str).split('/')
#     return odds_ratio_calc(bg_n, gene_list_n, int(gene_set_n), int(overlap_n))

# def odds_ratio_calc(bg_n, gene_list_n, gene_set_n, overlap_n): 
#     import scipy.stats as stats
#     table=np.array([[gene_set_n, bg_n-gene_set_n],[overlap_n, gene_list_n-overlap_n]])
#     oddsratio, pvalue = stats.fisher_exact(table)
#     return (1/oddsratio)


# configs

# input
query_genes_path = snakemake.input['query_genes']
background_genes_path = snakemake.input['background_genes']
database_path = snakemake.input['database']

# output
result_path = snakemake.output['result_file']

# parameters
db = snakemake.params["database"]

dir_results = os.path.dirname(result_path)

if not os.path.exists(dir_results):
    os.mkdir(dir_results)

# check if genes file exists & load or handle exception
if os.path.exists(query_genes_path):
    genes = open(query_genes_path, "r")
    gene_list = genes.read()
    gene_list = gene_list.split('\n')
    gene_list.remove('')
    genes.close()
else:
    with open(os.path.join(dir_results,"no_genes_found.txt"), 'w') as f:
        f.write('no genes found')
    quit()

# load background genes
bg_file = open(background_genes_path, "r")
background = bg_file.read()
background = background.split('\n')
background.remove('')
bg_file.close()

# move on if query-genes are empty
if len(gene_list)==0:
    open(result_path, mode='a').close()
    sys.exit(0)

# load database .pkl file
# with open(enrichr_databases, 'rb') as f:
#     db_dict = pickle.load(f)
# with open(database_path) as json_file:
#     db_dict = json.load(json_file)

# load database GMT file
db_dict = gp.parser.read_gmt(database_path)
    
# convert gene lists and database to upper case
gene_list=[str(x).upper() for x in list(gene_list)]
background=[str(x).upper() for x in list(background)]
db_dict = {key: [ele.upper() for ele in db_dict[key] ] for key in db_dict}
    
# count number of background genes for odds-ratio calculation
# bg_n = len(background)

# if background-genes are empty provide number of genes as 20,000 as heuristic
# this excpetion only occurs if background region set exceeds 500,000 regions
if len(background)==0:
    background = 20000

# perform ORA (hypergeometric test) in database using GSEApy (barplots are generated automatically)
try:
    res = gp.enrich(gene_list=gene_list,
                     gene_sets=db_dict,
                     background=background,
                     outdir=None,#os.path.join(dir_results),
                     top_term=25,
                     cutoff=0.05,
                     format='png',
                     verbose=True,
                    ).res2d
except ValueError:
    print("Result is empty")
    res = pd.DataFrame()

# move on if result is empty
if res.shape[0]==0:
    open(result_path, mode='a').close()
    sys.exit(0)

# annotate used gene set
res['Gene_set'] = db

# odds ratio calculation -> results still differ, but correlation is >0.99 and with large OR values the difference shows up at third decimal place
# gene_list_n=len(gene_list)
# res['Odds Ratio'] = res['Overlap'].apply(overlap_converter, args=(bg_n, gene_list_n))

# make column names language agnostic (i.e., R compatible)
column_names = list(res.columns.values)
column_names = [col.replace(" ","_") for col in column_names]
column_names = [col.replace("-","_") for col in column_names]
res.columns = column_names

# separate export
res.to_csv(result_path)

