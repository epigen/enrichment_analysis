# Genomic Region Set & (Ranked) Gene Set Enrichment Analysis & Visualization Snakemake Workflow for Human and Mouse Genomes.

Given **human (hg19 or hg38) or mouse (mm9 or mm10)** based genomic region sets (i.e., region sets) and/or (ranked) gene sets of interest and respective background region/gene sets, the enrichment within the configured databases is determined using LOLA, GREAT, GSEApy (over-represenation analysis (ORA) & preranked GSEA) and results saved as CSV files. Additionally, the most significant results are plotted for each region/gene set, database queried, and analysis performed. Finally, the results within the same "group" (e.g.,  stemming from the same DEA) are aggregated per database and analysis in summary CSV files and visualized using hierarchically clustered heatmaps and bubble-heatmap plots. For collaboration, communication and documentation of results, methods and workflow information a detailed self-contained HTML report can be generated.

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
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (.yaml file) or post execution. Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g. [X].

The outlined analyses were performed using the programming languages R (ver) [ref] and Python (ver) [ref] unless stated otherwise. All approaches statistically correct their results using expressed/accessible background genomic region/gene sets from the respective analyses that yielded the query region/gene sets.

**Genomic region set enrichment analyses**
**LOLA.** Genomic region set enrichment analysis was performed using LOLA (ver) [ref], which uses Fisher’s exact test. The following databases were queried [lola_dbs].

**GREAT.** Genomic region set enrichment analysis was performed using GREAT [ref] implemented with rGREAT (ver) [ref]. The following databases were queried [great_dbs].

Furthermore, genomic regions (query- and background-sets) were mapped to genes with GREAT and then analyzed as gene-sets as described below for a complementary/extended perspective.

**Gene set enrichment analyses**
GSEApy. Gene set enrichment analysis was performed using Enrichr [ref], which uses Fisher’s exact test, implemented with GSEApy (ver) [ref]. The following databases were queried [enrichr_dbs].

**Aggregation**
The results of all queries belonging to the same analysis [group] were aggregated by retaining all results of the union of all statistically significant [adjp_th] terms per query.

**Visualization**
Genomic region- and gene-set enrichment analysis results were visualized in the same way.

For each query set and database combination an enrichment dot plot was used to visualize the most important results.  The top [top_n] terms were ranked (along the y-axis) by the mean rank of statistical significance (p-value), effect-size (odds-ratio), and coverage/support (proportion of overlap) with the goal to make the results more balanced and interpretable. The significance (adjusted p-value) is denoted by the dot color, effect-size by the x-axis position, and coverage/support by the dot size.

The aggregated results per analysis [group] and database combination were visualized using hierarchically clustered heatmaps and bubble plots. The union of the top [top_terms_n] most significant terms per query were determined and their effect-size and significance were visualized as hierarchically clustered heatmaps (statistical significance of adjusted p-value < 0.05 was denoted by *). Furthermore, a hierarchically clustered overview bubble plot encoding both effect-size (size) and significance (color) is provided. All summary visualizations’ values were capped by [adjp_cap/or_cap] to avoid shifts in the coloring scheme caused by outliers.

**The analysis and visualizations described here were performed using a publicly available Snakemake (ver) [ref] workflow [ref - cite this workflow here].**


# Features
The three tools LOLA, GREAT and GSEApy are used for enrichment analysis. Databases to be queried can be configured (see ./config/config.yaml). All approaches statistically correct their results using the provided background region/gene sets.
- enrichment analysis:
    - region-set
        - [LOLA](http://lolaweb.databio.org/): Genomic Locus Overlap Enrichment Analysis is run locally. Required (cached) databases, which are downloaded automatically during the first run. [Supported databases](https://databio.org/regiondb) depend on the genome.
        - [GREAT](http://great.stanford.edu/public/html/index.php) Genomic Regions Enrichment of Annotations Tool is queried remotely (requires internet connection). [Supported databases](https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655440/Ontologies) depend on the genome.
    - gene-set
        - [GSEApy](https://gseapy.readthedocs.io/en/latest/) performs Fisher’s exact test and is run locally. Required databases are downloaded automatically during the first run (and saved as JSON files). All [Enrichr databases](https://maayanlab.cloud/Enrichr/#libraries) can be queried.
    - region-based gene-set
        - region-gene associations for each query and background region-set are obtained using GREAT
        - they are used for a complementary gene-set enrichment analysis using GSEApy
        - thereby an extended region-set enrichment perspective can be gained through association to genes by querying the same and/or more databases 
- group summary per queried database (aggregation of results)
    - results of all queries belonging to the same group are aggregated by taking the union of all statistically significant (adjp_th) terms per query
    - effect-size (or) and statistical significance (adjp) values are saved as CSV files
- visualization
    - region/gene-set specific enrichment dot plots are generated for each query and database combination
        - the top top_n terms are ranked (along the y-axis) by the mean rank of statistical significance (p_value), effect-size (odds_ratio), and coverage/support (overlap) with the goal to make the results more balanced and interpretable
        - significance (adjusted p-value) is presented by the dot color
        - effect-size is presented by the x-axis position
        - coverage/support is presented by the dot size
    - group summary
        - the union of the top (top_terms_n) most significant terms per query are determined and their effect-size (or) and significance (adjp) are visualized as hierarchically clustered heatmaps (statistical significance is denoted by \*).
        - a hierarchically clustered overview bubble plot encoding both effect-size (size) and significance (color) is provided as PNG and SVG.
        - all summary visualizations are configured to cap the values (adjp_cap/or_cap) to avoid shifts in the coloring scheme caused by outliers.

# Results
The result directory "enrichment_analysis" contains a folder for each query and group
- "query region/gene set" >> "tool" >> "database" containing:
    - CSV result file
    - enrichment dot plot (SVG and PNG)
- "group" >> "tool" >> "database" containing
    - CSV summary files for adjusted p-values and odds ratios
    - hierarchically clustered heatmaps visualizing the adjusted p-values and odds ratios of the top terms
    - summary bubble heatmap visualizing both (significance and effect size)

# Usage
Here are some tips for the usage of this workflow:
- Run the analysis on every query gene/region set of interest (eg results of differential analyses) with the respective background genes/regions (eg all expressed genes in the data or consensus regions).
- generate the [Snakemake Report](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html)
- look through the overview plots of your dedicated groups and queried databases in the report
- dig deeper by looking at the 
    - summary CSV file underlying the overview plot
    - enrichment plots for the individual query sets
- investigate interesting hits further by looking into the raw CSV output files of individual queries

# Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# Examples
--- COMING SOON ---

# Links
- [GitHub Repository](https://github.com/epigen/enrichment_analysis/)
- [GitHub Page](https://epigen.github.io/enrichment_analysis/)
- [Zenodo Repository (coming soon)]()
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/enrichment_analysis)
