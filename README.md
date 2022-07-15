# Gene Set & Genomic Region Set Enrichment Analysis & Visualization Snakemake Workflow for hg38 & mm10.

Given **hg38 or mm10** based gene sets and/or genomic region sets of interest and respective background gene/region sets, the enrichment within the configured databases is determined and results saved as .csv files. Additionally, the top 25 statistically significant results are plotted for each gene/region set and database queried. Finally, the results within the same "group" are aggregated per database and visualized using hierarchically clustered heatmaps and dotplots.

**If you use this workflow in a publication, don't forget to give credits to the authors by citing the URL of this (original) repository (and its DOI, see Zenodo badge above -> coming soon).**

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#methods)
  * [Features](#features)
  * [Usage](#usage)
  * [Configuration](#configuration)
  * [Examples](#examples)
  * [Links](#links)

# Authors
- [Stephan Reichl](https://github.com/sreichl)

# Software (TODO)
This project wouldn't be possible without the following software and their dependencies:

| Software       | Reference (DOI)                                   |
| :------------: | :-----------------------------------------------: |
| Enrichr        |  |
| ggplot2        | https://ggplot2.tidyverse.org/                    |
| GREAT          |        |
| GSEApy         |        |
| LOLA           |       |
| pandas         |       |
| pheatmap       | https://cran.r-project.org/package=pheatmap       |
| rGREAT         |       |
| Snakemake      | https://doi.org/10.12688/f1000research.29032.2    |

# Methods
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (.yaml file) or post execution. Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g. [X].

--- COMING SOON ---

# Features & Results (TODO: UPDATE)
the following steps are performed for each query region set (always with default settings):
- [LOLA](http://lolaweb.databio.org/) is run locally with LOLACore, LOLAJasper(Motifs) and LOLARoadmap(Epigenomics) (necessary cached databases are downloaded automatically)
- [GREAT](http://great.stanford.edu/public/html/index.php) is queried remotely (all available databases) and additionally used to determine region-gene associations (gene list is also saved)
- [Enrichr](https://maayanlab.cloud/Enrichr/) all available databases are queried in two different ways with the genes determined by GREAT:
    - R package [enrichR](https://cran.r-project.org/web/packages/enrichR/index.html): remotely  **without using a background gene/region set(!)**
    - python package [GSEApy](https://gseapy.readthedocs.io/en/latest/introduction.html): locally (databases are downloaded once and saved as a dictionary pickle file) **with a background gene set** obtained by querying GREAT with the background region set
- for each queried database a dotplot with the top 25 statistically significant hits (adjusted p-value < 0.05) is generated
    - the hits are ordered (along the y-axis) by the mean rank of p-value, odds ratio/foldchange, coverage/support with the goal to make the results more balanced and interpretable
    - p-value is presented by the dot color
    - odds ratio/foldchange is presented by the x-axis position
    - coverage/support is presented by the dot size

The results directory consists of one directory for each query region set, which contains a folder for each method:
- GREAT: for each queried database 2 files are generated
    - all results in a .tsv file
    - top 25 statistically significant hits as .svg dot plot file
    - additionally the following results are provided:
        - GREAT_job_description.txt file with the GREAT job query
        - GREAT_genes.txt file containing the list of genes associated with the query region set
        - GREAT_region_gene_assciations.pdf plot describing the region-gene association mapping
- GSEApy: for each queried database a folder is generated containing 3 files
    - all results in a .txt file (without odds ratio)
    - all results in a .csv file (including odds ratio)
    - top 25 statistically significant hits as .svg bar plot file
- LOLA: for each queried database (LOLACore, , LOLAJasper (Motifs) and LOLARoadmap (Epigenomics)) a folder is generated containing
    - all results in a .tsv file
    - top 25 statistically significant hits as .svg dot plot file
    - in case of LOLACore additional results are provided:
        - results from each collection individually as .tsv files
        - top 100 statistically significant hits grouped by cell type as .svg dot plot file

# Usage
Here are some tips for the usage of this workflow:
- Run the analysis on every query gene/region set of interest (eg results of differential analyses) with the respective background genes/regions (eg all expressed genes in the data or consensus regions).
- generate the report
- look through the overview plots in the report (includes a search function)
- dig deeper by looking at the enrichment plots for the individual query sets
- investigate interesting hits further by digging into the raw .csv output files

# Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# Examples
--- COMING SOON ---

# Links
- [GitHub Repository](https://github.com/epigen/enrichment_analysis/)
- [GitHub Page](https://epigen.github.io/enrichment_analysis/)
- [Zenodo Repository (coming soon)]()
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/enrichment_analysis)
