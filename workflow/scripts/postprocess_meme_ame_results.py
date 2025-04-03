import pandas as pd
import numpy as np
import os
import sys


########### FUNCTIONS ###########
def check_empty_tsv_minus_comments(file_path, comment_char='#'):
    with open(file_path, 'r') as file:
        for line in file:
            # Skip lines that start with the comment character
            if not line.strip().startswith(comment_char) and line.strip():
                return False  # Found non-comment, non-empty line
    return True  # All lines were either comments or empty


#### Process MEME-AME results to aggregate and clean the data.
input_file = snakemake.input['ame_results']
output_file = snakemake.output['aggregated_results']


# Ensure output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# quit early if file is empty
if check_empty_tsv_minus_comments(input_file):
    open(output_file, 'w').close()
    quit()
else:
    # Process the MEME-AME results
    try:
        # Read the file, skipping blank lines and comments
        df = pd.read_csv(input_file, sep='\t', comment='#', skip_blank_lines=True)

        # Drop any completely empty rows
        df.dropna(how='all', inplace=True)
        
        # Calculate effect size (i.e. log2 fold enrich)
        df["FoldEnrich"] = df["%TP"] / df["%FP"]
        
        # Clean up colum names so that it's R friendly
        df.columns = df.columns.str.replace("-", "_").str.replace("%", "perc")

        # Save the cleaned and aggregated results as a CSV
        df.to_csv(output_file, index=False)

        print(f"Processed MEME-AME results saved to {output_file}")

    except Exception as e:
        print(f"Error processing MEME-AME results: {e}")
        sys.exit(1)

