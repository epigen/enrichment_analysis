#!/bin/env python

# load libraries
import pandas as pd
import os
import numpy as np

# configs

# input
result_paths = snakemake.input['enrichment_results']

# output
results_all_path = snakemake.output['results_all']

# parameters
group = snakemake.wildcards["group"]
tool = snakemake.wildcards["tool"]
db = snakemake.wildcards["db"]

term_col = snakemake.config["column_names"][tool]["term"]
adjp_col = snakemake.config["column_names"][tool]["adj_pvalue"]

adjp_th = snakemake.config["adjp_th"][tool]

dir_results = os.path.dirname(results_all_path)
if not os.path.exists(dir_results):
    os.mkdir(dir_results)

# load results
results_list = list()
for result_path in result_paths:
    if os.stat(result_path).st_size != 0:
        tmp_name = os.path.basename(result_path).replace("_{}.csv".format(db),"")
        tmp_res = pd.read_csv(result_path, index_col=0)
        tmp_res['name'] = tmp_name
        results_list.append(tmp_res)
        
        
# move on if results are empty
if len(results_list)==0:
    open(results_all_path, mode='a').close()
    sys.exit(0)
        
# concatenate all results into one results dataframe
result_df = pd.concat(results_list, axis=0)

# save all enirchment results
result_df.to_csv(results_all_path)