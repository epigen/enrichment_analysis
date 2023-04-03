# Genomic Region Set & (Ranked) Gene Set Enrichment Analysis & Visualization Snakemake Workflow for Human and Mouse Genomes.

Given **human (hg19 or hg38) or mouse (mm9 or mm10)** based genomic region sets (i.e., region sets) and/or (ranked) gene sets of interest and respective background region/gene sets, the enrichment within the configured databases is determined using LOLA, GREAT, GSEApy (over-representation analysis (ORA) & preranked GSEA) and results saved as CSV files. Additionally, the most significant results are plotted for each region/gene set, database queried, and analysis performed. Finally, the results within the same "group" (e.g.,  stemming from the same DEA) are aggregated per database and analysis in summary CSV files and visualized using hierarchically clustered heatmaps and bubble plots. For collaboration, communication and documentation of results, methods and workflow information a detailed self-contained HTML report can be generated.

This workflow adheres to the module specifications of [MR. PARETO](https://github.com/epigen/mr.pareto), an effort to augment research by modularizing (biomedical) data science. For more details and modules check out the project's repository.

**If you use this workflow in a publication, don't forget to give credits to the authors by citing the URL of this (original) repository (and its DOI, see Zenodo badge above -> coming soon).**

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#methods)
  * [Features](#features)
  * [Results](#results)
  * [Usage](#usage)
  * [Configuration](#configuration)
  * [Examples](#examples)
  * [Links](#links)

# Authors
- [Stephan Reichl](https://github.com/sreichl)

# Software
This project wouldn't be possible without the following software and their dependencies:

| Software       | Reference (DOI)                                   |
| :------------: | :-----------------------------------------------: |
| Enrichr        | https://doi.org/10.1002/cpz1.90                   |
| ggplot2        | https://ggplot2.tidyverse.org/                    |
| GREAT          | https://doi.org/10.1371/journal.pcbi.1010378      |
| GSEA           | https://doi.org/10.1073/pnas.0506580102           |
| GSEApy         | https://doi.org/10.1093/bioinformatics/btac757    |
| LOLA           | https://doi.org/10.1093/bioinformatics/btv612     |
| pandas         | https://doi.org/10.5281/zenodo.3509134            |
| pheatmap       | https://cran.r-project.org/package=pheatmap       |
| rGREAT         | https://doi.org/10.1093/bioinformatics/btac745    |
| Snakemake      | https://doi.org/10.12688/f1000research.29032.2    |

# Methods
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (workflow/envs/\*.yaml files) or post execution (results_dir/envs/enrichment_analysis/\*.yaml files). Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g. [X].

The outlined analyses were performed using the programming languages R (ver) [ref] and Python (ver) [ref] unless stated otherwise. All approaches statistically correct their results using expressed/accessible background genomic region/gene sets from the respective analyses that yielded the query region/gene sets.

**Genomic region set enrichment analyses**

**LOLA.** Genomic region set enrichment analysis was performed using LOLA (ver) [ref], which uses Fisher’s exact test. The following databases were queried [lola_dbs].

**GREAT.** Genomic region set enrichment analysis was performed using GREAT [ref] implemented with rGREAT (ver) [ref]. The following databases were queried [great_dbs].

Furthermore, genomic regions (query- and background-sets) were mapped to genes using GREAT and then analyzed as gene-sets as described below for a complementary and extended perspective.

**Gene set enrichment analyses (GSEA)**

**Over-representation analysis (ORA).** Gene set ORA was performed using Enrichr [ref], which uses Fisher’s exact test (i.e., hypergeometric test), implemented with GSEApy's (ver) [ref] function _enrich_. The following databases were queried [enrichr_dbs][local_gmt_dbs][local_json_dbs].

**Preranked GSEA.** Preranked GSEA was performed using GSEA [ref], implemented with GSEApy's (ver) [ref] function _prerank_. The following databases were queried [enrichr_dbs][local_gmt_dbs][local_json_dbs].

**Aggregation**
The results of all queries belonging to the same analysis [group] were aggregated by method and database. Additionally, we filtered the results by retaining only the union of terms that were statistically significant (i.e. [adj_pvalue]<[adjp_th]) in at least one query.

**Visualization**
All analysis results were visualized in the same way.

For each query, method and database combination an enrichment dot plot was used to visualize the most important results.  The top [top_n] terms were ranked (along the y-axis) by the mean rank of statistical significance ([p_value]), effect-size ([effect_size]), and overlap ([overlap]) with the goal to make the results more balanced and interpretable. The significance (adjusted p-value) is denoted by the dot color, effect-size by the x-axis position, and overlap by the dot size.

The aggregated results per analysis [group], method and database combination were visualized using hierarchically clustered heatmaps and bubble plots. The union of the top [top_terms_n] most significant terms per query were determined and their effect-size and significance were visualized as hierarchically clustered heatmaps, and statistical significance ([adj_pvalue] < [adjp_th]) was denoted by \*. Furthermore, a hierarchically clustered bubble plot encoding both effect-size (color) and statistical significance (size) is provided, with statistical significance denoted by \*. All summary visualizations’ values were capped by [adjp_cap]/[or_cap]/[nes_cap] to avoid shifts in the coloring scheme caused by outliers.

**The analysis and visualizations described here were performed using a publicly available Snakemake (ver) [ref] workflow [ref - cite this workflow here].**


# Features
The three tools LOLA, GREAT and GSEApy (over-representation analysis (ORA) & preranked GSEA) are used for various enrichment analyses. Databases to be queried can be configured (see ./config/config.yaml). All approaches statistically correct their results using the provided background region/gene sets.
- enrichment analysis methods:
    - region-set
        - [LOLA](http://bioconductor.org/packages/release/bioc/html/LOLA.html): Genomic Locus Overlap Enrichment Analysis is run locally. Required (cached) databases, which are downloaded automatically during the first run. [Supported databases](https://databio.org/regiondb) depend on the genome (lola_dbs).
        - [GREAT](https://doi.org/10.1371/journal.pcbi.1010378) using [rGREAT](http://bioconductor.org/packages/release/bioc/html/rGREAT.html): Genomic Regions Enrichment of Annotations Tool is queried remotely (requires a working internet connection). [Supported databases](https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655440/Ontologies) depend on the genome (great_dbs).
            - query region sets with >500,000 regions are [not supported](https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655402/File+Size) and empty output files are generated to satisfy Snakemake
            - background region sets with >1,000,000 are [not supported](https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655402/File+Size) and the whole genome is used as background
    - gene-set over-representation analysis (ORA_GSEApy)
        - [GSEApy](https://gseapy.readthedocs.io/en/latest/) enrich() function performs Fisher’s exact test (i.e., hypergeoemtric test) and is run locally.
    - region-based gene-set over-representation analysis (ORA_GSEApy)
        - region-gene associations for each query and background region-set are obtained using GREAT.
        - they are used for a complementary ORA using GSEApy.
        - thereby an extended region-set enrichment perspective can be gained through association to genes by querying the same and/or more databases, that are not supported/provided by region-based tools. 
        - limitation: if the background region set exceeds GREAT's capacities (i.e., 1,000,000 regions), no background gene list is generated and background gene number (bg_n) of 20,000 is used in the ORA.
    - preranked gene-set enrichment analysis (preranked_GSEApy)
        - [GSEApy](https://gseapy.readthedocs.io/en/latest/) prerank() function performs [preranked GSEA](https://doi.org/10.1073/pnas.0506580102) and is run locally.
        - no duplicates allowed: only entries with the largest absolute score are kept.
- resources (databases) for both gene-based analyses are downloaded (Enrichr) or copied (local files) and saved as JSON files in /resources
    - all [Enrichr databases](https://maayanlab.cloud/Enrichr/#libraries) can be queried (enrichr_dbs).
    - local JSON database files can be queried (local_json_dbs).
    - local GMT database files (e.g., from [MSigDB](http://www.gsea-msigdb.org/gsea/msigdb)) can be queried (local_gmt_dbs).
- group aggregation of results per method and database
    - results of all queries belonging to the same group are aggregated per method (e.g., ORA_GSEApy) and database (e.g., GO_Biological_Process_2021) by concatenation and saved as a long-format table (CSV).
    - a filtered version taking the union of all statistically significant (i.e., adjusted p-value <{adjp_th}) terms per query is also saved as CSV file.
- visualization
    - region/gene-set specific enrichment dot plots are generated for each query, method and database combination
        - the top {top_n} terms are ranked (along the y-axis) by the mean rank of statistical significance ({p_value}), effect-size ({efect_size} e.g., log2(odds ratio) or normalized enrichemnt scores), and overlap ({overlap} e.g., coverage or support) with the goal to make the results more balanced and interpretable
        - significance (adjusted p-value) is presented by the dot color
        - effect-size is presented by the x-axis position
        - overlap is presented by the dot size
    - group summary/overview
        - the union of the top {top_terms_n} most significant terms per query, method, and database within a group is determined. 
        - their effect-size (effect) and statistical significance (adjp) are visualized as hierarchically clustered heatmaps, with statistical significance denoted by \* (PDF).
        - a hierarchically clustered bubble plot encoding both effect-size (color) and significance (size) is provided, with statistical significance denoted by \* (PNG and SVG).
        - all summary visualizations are configured to cap the values ({adjp_cap}/{or_cap}/{nes_cap}) to avoid shifts in the coloring scheme caused by outliers.

# Results
The result directory {result_path}/enrichment_analysis contains a folder for each region/gene-set {query} and {group}
- {query}/{method}/{database}/ containing:
    - result table (CSV file): {query}\_{database}.csv
    - enrichment dot plot (SVG and PNG): {query}\_{database}.{svg|png}
- {group}/{method}/{database}/ containing
    - aggregated result table (CSV file): {group}\_{database}\_all.csv
    - filtered aggregated result table (CSV file): {group}\_{database}\_sig.csv
    - hierarchically clustered heatmaps visualizing statistical significance and effect-sizes of the top {top_terms_n} terms (PDF): {group}\_{database}\_{adjp|effect}\_hm.pdf
    - hierarchically clustered bubble plot visualizing statistical significance and effect-sizes simultaneously (PNG and SVG):  {group}\_{database}\_summary.{svg|png}

# Usage
Here are some tips for the usage of this workflow:
- Run the analysis on every query gene/region set of interest (e.g., results of differential analyses) with the respective background genes/regions (e.g., all expressed genes in the data or consensus regions).
- generate the [Snakemake Report](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html)
- look through the overview plots of your dedicated groups and queried databases in the report
- dig deeper by looking at the 
    - aggregated result table underlying the summary/overview plot
    - enrichment plots for the individual query sets
- investigate interesting hits further by looking into the individual query result tables.

# Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# Examples
We provide four example queries:
- three are region-sets from a [LOLA Vignette](http://code.databio.org/LOLA/articles/usingLOLACore.html). Download the example data by following these [instructions](./.test/data/example_data_download_instructions.txt).
- one is a preranked gene-score set derived from the GDS289 [fgsea R package example data](https://github.com/ctlab/fgsea/blob/master/inst/extdata/GDS289.tsv) (score=-log10(p-value)\*sign(lfc)).

We provide two local example databases
- [MSigDB Human C2](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C2) in GMT format.
- A renamed copy of [Enrichr's WikiPathways 2019 Human](https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2019_Human) in the JSON format.

Follow these steps to run the complete analysis:
1. activate your snakemake conda environment
    ```console
    conda activate snakemake
    ```
2. enter the workflow directory
    ```console
    cd enrichment_analysis
    ```
3. run a snakemake dry-run (-n flag) using the provided configuration to check if everything is in order
    ```console
    snakemake -p --use-conda --configfile .test/config/example_enrichment_analysis_config.yaml -n
    ```
4. run the workflow
    ```console
    snakemake -p --use-conda --configfile .test/config/example_enrichment_analysis_config.yaml
    ```
5. generate report
    ```console
    snakemake --report .test/report.html --configfile .test/config/example_enrichment_analysis_config.yaml
    ```

# Links
- [GitHub Repository](https://github.com/epigen/enrichment_analysis/)
- [GitHub Page](https://epigen.github.io/enrichment_analysis/)
- [Zenodo Repository (coming soon)]()
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/enrichment_analysis)

# Resources
- Web versions of the used tools
    - [LOLA](http://lolaweb.databio.org/)
    - [GREAT](http://great.stanford.edu/public/html/index.php)
    - [Enrichr](https://maayanlab.cloud/Enrichr/)
- Gene-set databases
    - [GREAT databases](https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655440/Ontologies)
    - [Enrichr gene-set databases](https://maayanlab.cloud/Enrichr/#libraries)
    - [The Molecular Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb/)

# Publications
The following publications successfully used this module for their analyses.
- ...