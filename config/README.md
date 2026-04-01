You need one configuration file and one annotation file to run the complete workflow. You can use the provided example as starting point. If in doubt read the comments in the config and/or try the default values.

- project configuration (config/config.yaml): different for every project/dataset and configures the analyses to be performed and databases to be used. The fields are described within the file.
- annotation: CSV file consisting of **5 mandatory columns**
    - name: unique(!) name of the query gene/region set
    - features_path: one of the following
        - path to a query genomic region set as [**.bed**](https://en.wikipedia.org/wiki/BED) file, which will be analyzed using LOLA, GREAT, pycisTarget and ORA_GSEApy, RcisTarget (using the region-gene association provided by GREAT)
        - path to a query gene set as **.txt** file with one gene per line, which will be analyzed using ORA_GSEApy and RcisTarget
        - path to a query (preranked) gene-score table with a header as **.csv** file, where the first column consists of gene-symbols/names and the second of corresponding gene-scores (e.g., from differential expression analysis results), which will be analyzed using preranked_GSEApy
    - background_name: name of the background gene/region set (only required for region- and gene-sets, leave empty for gene-score tables)
    - background_path: path to the background/universe gene/region-set as .txt/.bed file (only required for region- and gene-sets, leave empty for gene-score tables)
    - group: enrichment results are aggregated and visualized per analysis and database based on this group variable (e.g., gene/region-sets resulting from the same analysis)

For **region-set inputs** (`features_path` or `background_path` ending in `.bed`),the workflow expects **[standard BED coordinates](https://en.wikipedia.org/wiki/BED)**.
- **BED is 0-based, start-inclusive, end-exclusive**, i.e. in interval notation: **`[start, end)`**.
- This means:
    - the first base of a chromosome is `0`
    - a region with `start=0` and `end=100` spans exactly 100 bases: positions `0` to `99`
- additionally, the workflow requires at least the first **3 BED columns**: `chrom`, `start`, `end`.

Set workflow-specific `resources` or command line arguments (CLI) in the workflow profile `workflow/profiles/default.config.yaml`, which supersedes global Snakemake profiles.
