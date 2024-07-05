
# libraries
import os
import pycistarget
from pycistarget.input_output import read_hdf5

# configs

# input
motif_hdf5_path = snakemake.input['motif_hdf5']

# output
motif_csv_path = snakemake.output['motif_csv']

# parameters
region_set_name = snakemake.wildcards["region_set"]
term_col = snakemake.config["pycistarget_parameters"]["annotations_to_use"][0]

# quit early if file is empty
if os.path.getsize(motif_hdf5_path) == 0:
    open(motif_csv_path, 'w').close()
    quit()

# load pycisTarget results from hdf5
results = read_hdf5(motif_hdf5_path)

# extract results
results_df = results[region_set_name].motif_enrichment

# reformat
results_df.index.name = "motif"
results_df.reset_index(inplace=True)
results_df["description"] = results_df["motif"] + "(" + results_df[term_col] + ")"

# save motif enrichments as CSV for downstream processing and plotting
results_df.to_csv(motif_csv_path)
