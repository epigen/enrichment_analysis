# libraries
import os
import pandas as pd
import pycistarget
from pycistarget.input_output import read_hdf5

# configs

# input
motif_hdf5_path = snakemake.input['motif_hdf5']

# output
motif_csv_path = snakemake.output['motif_csv']
hits_csv_path = snakemake.output['hits_csv']
cistrome_csv_path = snakemake.output['cistrome_csv']

# parameters
region_set_name = snakemake.wildcards["region_set"]
term_col = snakemake.config["pycistarget_parameters"]["annotations_to_use"][0]

# quit early if file is empty
if os.path.getsize(motif_hdf5_path) == 0:
    open(motif_csv_path, 'w').close()
    open(hits_csv_path, 'w').close()
    open(cistrome_csv_path, 'w').close()
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


# Save motif hits
motifs_df = pd.DataFrame(results[region_set_name].motif_hits)
motifs_df["database"] = motifs_df["database"].apply(
    lambda x: ";".join(x) if type(x) == list else ""
)
motifs_df["region_set"] = motifs_df["region_set"].apply(
    lambda x: ";".join(x) if type(x) == list else ""
)
motifs_df.to_csv(hits_csv_path)

# Save cistrome hits
cistrome_df = pd.DataFrame(results[region_set_name].cistromes)
cistrome_df["database"] = cistrome_df["database"].apply(
    lambda x: ";".join(x) if type(x) == list else ""
)
cistrome_df["region_set"] = cistrome_df["region_set"].apply(
    lambda x: ";".join(x) if type(x) == list else ""
)
cistrome_df.to_csv(cistrome_csv_path)
