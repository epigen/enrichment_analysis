#!/bin/env python

# load libraries
import pandas as pd
import os
import numpy as np

# configs
result_paths = snakemake.input['enrichment_results']

results_all_path = snakemake.output['results_all']
results_sig_path = snakemake.output['results_sig']

group = snakemake.wildcards["group"] #"testgroup"
tool = snakemake.wildcards["tool"] #"GSEApy"
db = snakemake.wildcards["db"] #"WikiPathways_2019_Human"

term_col = snakemake.config["column_names"][tool]["term"] #'Term'
adjp_col = snakemake.config["column_names"][tool]["adj_pvalue"] #'Adjusted_P_value'
or_col = snakemake.config["column_names"][tool]["odds_ratio"] #'Odds_Ratio'

adjp_th = snakemake.config["adjp_th"][tool] #0.05

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
    open(results_sig_path, mode='a').close()
    sys.exit(0)
        
# concatenate all results into one results dataframe
result_df = pd.concat(results_list, axis=0)

# save all enirchment results
result_df.to_csv(results_all_path)

# find union of statistically significant terms
sig_terms = result_df.loc[result_df[adjp_col]<adjp_th, term_col].unique()

# filter by significant terms
result_sig_df = result_df.loc[result_df[term_col].isin(sig_terms), :]

# save filtered enirchment results by significance
result_sig_df.to_csv(results_sig_path)
