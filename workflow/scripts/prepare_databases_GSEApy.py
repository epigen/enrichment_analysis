#!/bin/env python       
import json
import gseapy as gp
import os
import shutil
    
# config

# input
db_path = snakemake.input[0]

# output
results_path = snakemake.output['db_file']

# parameters
db = snakemake.wildcards["database"]


# if GMT, just copy
if db_path.endswith('.gmt'):
    shutil.copy(db_path, results_path)
elif db_path.endswith('.json'):
    # JSON load and save as GMT
    with open(db_path, 'r') as f:
        data = json.load(f)
    
    with open(results_path, 'w') as f:
        for key, values in data.items():
            f.write(f"{key}\t\t" + "\t".join(values) + "\n")
else:
    print("Error: Please provide a .gmt or .json database file.")
