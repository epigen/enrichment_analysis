You need one configuration file and one annotation file to run the complete workflow. You can use the provided example as starting point. Always use absolute paths. If in doubt read the comments in the config and/or try the default values.

- project configuration (config/config.yaml): different for every project/dataset and configures the analyses to be performed and databases to be used. The fields are described within the file.
- annotation: CSV file consisting of 5 mandatory columns
    - name: name of the query gene/region set, has to be unique
    - features_path: path to a query gene set as .txt with one gene per line OR query region set as .bed file
    - background_name: name of the background gene/region set
    - background_path: path to the background/universe gene/region set as .txt/.bed file
    - group: enrichment results are aggregated and visualized per database based on this group variable (eg gene/region sets from the same analysis)
