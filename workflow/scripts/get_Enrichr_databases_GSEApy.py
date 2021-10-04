#!/bin/env python
import gseapy as gp
import pickle
import os
    
# parameters from snakemake
results_path = snakemake.output['result_file']

# get enrichr database names
databases = gp.get_library_name()

# download all databases to local dict
db_dict = dict()
for db in databases:
    db_dict[db]=gp.parser.gsea_gmt_parser(db, min_size=0, max_size=100000,)

# save as pickle 
with open(results_path, 'wb') as f:
        pickle.dump(db_dict, f, pickle.HIGHEST_PROTOCOL)