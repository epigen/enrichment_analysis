
# libraries
# from pathlib import Path
# import sys
import pyranges as pr
# import os
# import matplotlib.pyplot as plt
import pickle

import pycistarget
from pycistarget.motif_enrichment_cistarget import *

# configs

# input
regions_path = snakemake.input["regions"]
ctx_db_path =  snakemake.input["ctx_db"]
motif2tf_path = snakemake.input["motif2tf"]

# output
motif_csv_path = snakemake.output['motif_csv']
motif_html_path = snakemake.output['motif_html']
motif_pickle_path = snakemake.output['motif_pickle']

# parameters
genome = snakemake.config["genome"]
region_set_name = snakemake.wildcards["region_set"]
pycistarget_params = snakemake.config["pycistarget_parameters"]
threads = snakemake.threads

if genome=="hg19" or genome=="hg38":
    specie = 'homo_sapiens'
elif genome=="mm9" or genome=="mm11":
    specie = 'mus_musculus'
else:
    print("Provide supported genome.")

# load region set
region_set = {region_set_name: pr.read_bed(regions_path)}

# run pycisTarget
cistarget_dict = run_cistarget(ctx_db = ctx_db_path,
                               region_sets = region_set,
                               specie = specie,
                               fraction_overlap_w_cistarget_database = pycistarget_params["fraction_overlap_w_cistarget_database"],
                               auc_threshold = pycistarget_params["auc_threshold"],
                               nes_threshold = pycistarget_params["nes_threshold"],
                               rank_threshold = pycistarget_params["rank_threshold"],
                               motif_similarity_fdr = pycistarget_params["motif_similarity_fdr"],
                               path_to_motif_annotations = motif2tf_path,
                               annotations_to_use = pycistarget_params["annotations_to_use"],
                               annotation_version = pycistarget_params["annotation_version"],
                               n_cpu = threads, 
                               _temp_dir = pycistarget_params["temp_dir"]
                              )

# save motif enrichments
# as CSV
cistarget_dict[region_set_name].motif_enrichment.to_csv(motif_csv_path)
# as HTML
cistarget_dict[region_set_name].motif_enrichment.to_html(open(motif_html_path, 'w'), escape=False, col_space=80)
# as pickle
with open(motif_pickle_path, 'wb') as f:
    pickle.dump(cistarget_dict, f)    
