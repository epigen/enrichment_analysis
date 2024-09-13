You need one configuration file and one annotation file to run the complete workflow. You can use the provided example as starting point. If in doubt read the comments in the config and/or try the default values.

- project configuration (config/config.yaml): different for every project/dataset and configures the analyses to be performed and databases to be used. The fields are described within the file.
- annotation: CSV file consisting of **5 mandatory columns**
    - name: unique(!) name of the query gene/region set
    - features_path: one of the following
        - path to a query region set as **.bed** file -> will be analyzed using LOLA, GREAT, pycisTarget and ORA_GSEApy, RcisTarget (using the region-gene association provided by GREAT)
        - path to a query gene set as **.txt** file with one gene per line -> will be analyzed using ORA_GSEApy and RcisTarget
        - path to a query (preranked) gene-score table with a header as **.csv** file, where the first column consists of gene-symbols/names and the second of corresponding gene-scores (e.g., from differential expression analysis results) -> will be analyzed using preranked_GSEApy
    - background_name: name of the background gene/region set (only required for region- and gene-sets, leave empty for gene-score tables)
    - background_path: path to the background/universe gene/region-set as .txt/.bed file (only required for region- and gene-sets, leave empty for gene-score tables)
    - group: enrichment results are aggregated and visualized per analysis and database based on this group variable (e.g., gene/region-sets resulting from the same analysis)

Set workflow-specific `resources` or command line arguments (CLI) in the workflow profile `workflow/profiles/default.config.yaml`, which supersedes global Snakemake profiles.