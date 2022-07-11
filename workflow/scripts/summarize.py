#!/bin/env python

# load libraries
import pandas as pd
import os
import numpy as np


# configs
result_paths = snakemake.input['enrichment_results'] #[os.path.join("/research/home/sreichl/projects/genomic_region_enrichment/test/enrichment_analysis/IRF8_down/GSEApy/WikiPathways_2019_Human/IRF8_down_WikiPathways_2019_Human.csv")]

adjp_path = snakemake.output['adjpvalues']
or_path = snakemake.output['oddsratios']

group = snakemake.wildcards["group"] #"testgroup"
tool = snakemake.wildcards["tool"] #"GSEApy"
db = snakemake.wildcards["db"] #"WikiPathways_2019_Human"

term_col = snakemake.config["column_names"][tool]["term"] #'Term'
adjp_col = snakemake.config["column_names"][tool]["adj_pvalue"] #'Adjusted_P_value'
or_col = snakemake.config["column_names"][tool]["odds_ratio"] #'Odds_Ratio'

adjp_th = snakemake.config["adjp_th"][tool] #0.05

dir_results = os.path.dirname(adjp_path)
if not os.path.exists(dir_results):
    os.mkdir(dir_results)

# load results
results_dict = dict()
for result_path in result_paths:
    tmp_name = os.path.basename(result_path).replace("_{}.csv".format(db),"")
    
    if os.stat(result_path).st_size != 0:
        results_dict[tmp_name] = pd.read_csv(result_path)
    
# find union of stat. sign. terms
sign_terms = list()
for result in results_dict.keys():
#     tmp_terms = results_dict[result][term_col][results_dict[result][adjp_col]<adjp_th].to_list()
    tmp_terms = results_dict[result].loc[results_dict[result][adjp_col]<adjp_th, term_col].to_list()
    sign_terms = sign_terms + tmp_terms
sign_terms = list(set(sign_terms))

# aggregate data of stat. sign. terms in summary data frames (adjp and or)
adjp_df=pd.DataFrame(index=sign_terms, columns=results_dict.keys())
or_df=pd.DataFrame(index=sign_terms, columns=results_dict.keys())
for result in results_dict.keys():
    # determine sign. terms within the result
    idx_intersect = list(set(results_dict[result].set_index(term_col).index).intersection(set(sign_terms)))
    
    # fill data frames
    adjp_df.loc[idx_intersect, result]=results_dict[result].set_index(term_col).loc[idx_intersect,adjp_col]
    or_df.loc[idx_intersect, result]=results_dict[result].set_index(term_col).loc[idx_intersect,or_col]

# save results
adjp_df.to_csv(adjp_path)
or_df.to_csv(or_path)
