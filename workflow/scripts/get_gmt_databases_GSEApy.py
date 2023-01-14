#!/bin/env python       
import json
import gseapy as gp
import os
    
# config

# input
gmt_path = snakemake.input[0]

# output
results_path = snakemake.output['result_file']

# parameters
db = snakemake.params["database"]

# download the GMT database to local dict
db_dict = gp.parser.read_gmt(gmt_path)

# save as json
with open(results_path, 'w') as fp:
    json.dump(db_dict, fp,  indent=4)