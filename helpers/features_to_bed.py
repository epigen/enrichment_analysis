# convert feature lists to BED files for genomic region enrichment analysis

#### libraries ####
import os
import pandas as pd

#### configs ####
# inputs
consensus_annotation_path = os.path.join("path/to/consensus_annotation.csv")
features_txt_path = os.path.join("path/to/features.txt")

# output
features_bed_path = os.path.join("path/to/features.bed")

# parameters
# columns in consensus_annotation
region_col = "peak_id"
chr_col = "gencode_chr"
start_col = "gencode_start"
end_col = "gencode_end"


#### load data ####
consensus_df = pd.read_csv(consensus_annotation_path)
features_df = pd.read_csv(features_txt_path, header=None, names=[region_col])

#### merge and select relevant columns ####
merged_df = pd.merge(features_df, consensus_df, on=region_col, how="inner")
# Select and reorder columns to match the BED format (chr, start, end, name)
bed_df = merged_df[[chr_col, start_col, end_col, region_col]]

#### save result ####
bed_df.to_csv(features_bed_path, sep="\t", header=False, index=False)
