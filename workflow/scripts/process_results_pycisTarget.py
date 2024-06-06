
# libraries
import pycistarget
from pycistarget.input_output import read_hdf5

# configs

# input
motif_hdf5_path = snakemake.input['motif_hdf5']

# output
motif_csv_path = snakemake.output['motif_csv']

# parameters
region_set_name = snakemake.wildcards["region_set"]

# load pycisTarget results from hdf5
results = read_hdf5(motif_hdf5_path)

# extract results
results_df = results[region_set_name].motif_enrichment

# reformat

# save motif enrichments as CSV for downstream processing and plotting
results_df.to_csv(motif_csv_path)
