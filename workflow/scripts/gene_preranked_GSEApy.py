#!/bin/env python

# load libraries
import pandas as pd
import json
# import pickle
import os
import numpy as np
import gseapy as gp
import sys


# configs

# input
query_genes_path = snakemake.input['query_genes']
database_path = snakemake.input['database']

# output
result_path = snakemake.output['result_file']

# parameters
db = snakemake.params["database"]

dir_results = os.path.dirname(result_path)

if not os.path.exists(dir_results):
    os.mkdir(dir_results)
    

# load gene-score file
genes = pd.read_csv(query_genes_path, index_col=0)

# load database JSON file
with open(database_path) as json_file:
    db_dict = json.load(json_file)
    
# convert all genes to upper case
genes.index = [str(x).upper() for x in list(genes.index)]
db_dict = {key: [ele.upper() for ele in db_dict[key] ] for key in db_dict}

# remove duplicates: only retain largest absolute value
# sort by absolute value of "score" column (i.e., first column)
genes = genes.iloc[abs(genes.iloc[:, 0]).argsort()[::-1]]
# drop duplicates based on index
genes = genes[~genes.index.duplicated(keep='first')]

# run prerank GSEA of database with GSEApy
res = gp.prerank(rnk=genes, 
                 gene_sets=db_dict,
                 #threads=4,
                 min_size=1, # Minimum allowed number of genes from gene set also the data set. Default: 15.
                 max_size=100000, # Maximum allowed number of genes from gene set also the data set. Defaults: 500.
                 permutation_num=1000, # Number of permutations. Reduce number to speed up testing;  Default: 1000. Minimial possible nominal p-value is about 1/nperm.
                 outdir=os.path.join(dir_results),
                 graph_num = 25, # Plot graphs for top sets of each phenotype.
                 format='png',
                 seed=42,
                 verbose=True,
                ).res2d

# move on if result is empty
if res.shape[0]==0:
    open(result_path, mode='a').close()
    sys.exit(0)

# annotate used gene set
res['Gene_set'] = db

# make column names language agnostic (i.e., R compatible)
column_names = list(res.columns.values)
column_names = [col.replace(" %","") for col in column_names]
column_names = [col.replace(" ","_") for col in column_names]
column_names = [col.replace("-","_") for col in column_names]
res.columns = column_names

# separate export
res.to_csv(result_path)
