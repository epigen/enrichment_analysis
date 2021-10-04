# Genomic Region Enrichment Analysis Snakemake Workflow using LOLA, GREAT and Enrichr for hg38 & mm10.

Given **hg38 or mm10** based genomic region sets of interest and a background region set (as .bed files), the enrichment within mutliple databases is determined and results saved as .tsv/.csv files. Additionally, the top 25 statistically significant results are plotted for each database queried.

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

# Features
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

# Recommended Usage
1. Run the analysis on every query region set of interest (eg results of differential analyses) with the respective background regions (eg filtered consensus regions)
2. generate the report
3. look through the generated plots in the report (includes a search function)
4. investigate interesting hits further by digging into the raw .tsv output files

# Installation (<10 minutes)
1. install snakemake, which requires conda & mamba, according to the [documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
2. clone/download this repository (eg git clone https://github.com/sreichl/genomic_region_enrichment.git)

All software/package dependencies are installed and managed automatically via Snakemake and conda.

# Configuration
2 configuration files are needed. Always use absolute paths.
- general configuration (config/config.yaml):
    - region_annotation: path to the region annotation .csv file (see below)
    - results_dir: path to the directory for the enrichment results (results directory)
    - genome: genome reference that was used ('hg38' or 'mm10')
    - parameters for cluster execution: partition, memory, threads (if in doubt try the default values)
    - an example is provided in the repository
- region annotation: CSV file with 4 columns (mandatory)
    - name: name of the query region set
    - regions_bed: path to the query region set .bed file
    - background_name: name of the background region set
    - background_bed: path to the background/universe region set .bed file

# Execution
## 1. Change working directory & activate conda environment
Execute always from within top level of the pipeline directory (ie genomic_region_enrichment/).
Snakemake commands only work from within the snakemake conda environment.
```
cd genomic_region_enrichment
conda activate snakemake
```
## 2. Execute a dry-run
command for a dry-run with option -n (-p makes Snakemake print the resulting shell command for illustration)
```
snakemake -p -n
```
## 3. Execute workflow local or on a cluster
### 3a. Local execution
command for execution with one core
```
snakemake -p -j1 --use-conda
```
### 3b. Cluster execution
command for **vanilla cluster execution** on cluster engines that support shell scripts and have access to a common filesystem, (e.g. the Sun Grid Engine), more info in the [Snakemake Cluster Execution documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)
```
snakemake -p --use-conda --cluster qsub -j 32
```

command for **cluster execution by using --profile**, submits every task as separate job with dependencies
```
snakemake -p --use-conda --profile config/slurm.cemm
```
the profile for CeMM's slurm environment is provided in the config/ directory, of note: 
- the number of jobs in the slurm.cemm/config.yaml should be set as high as necessary, because arrayed job subsmission does not work (yet) and the scheduler (eg SLURM) should take care of the priorization
- jobs which dependencies can never be fulfilled are automatically removed from the queue

If you are using another setup get your cluster execution profile here: [The Snakemake-Profiles project](https://github.com/snakemake-profiles/doc)

# Use as module in another Snakemake workflow (soon)
- [https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-modules)
- [https://slides.com/johanneskoester/snakemake-6#/8](https://slides.com/johanneskoester/snakemake-6#/8)

# Report
command for report generation, executed from within the directory genomic_region_enrichment/ (can take a few minutes, depending on the number of tested region sets)
```
snakemake --report /absolute/path/to/report.zip
```

The command creates a self contained HTML based report in a zip archive containing the following sections:
- Workflow: interactive rulegraph to recapitulate individual steps, used software and conrete code (reproducibility)
- Statistics: duration and timing of individual steps
- Configuration: used general configuration (accountability)
- Results
    - Configuration: the used region annotation configuration file
    - a section for each query region set (row in the region annotation): each enrichment analysis method (LOLA, GREAT, Enrichr/GSEApy) has its own searchable sub-section within the query region set section consisting of all generated dot plots of significant hits

# Results
The results directory consists of one directory for each query region set, which contains a folder for each method:
- Enrichr: for each queried database 2 files are generated
    - all results in a .tsv file
    - top 25 statistically significant hits as .svg dot plot file
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


# Examples (coming soon)

# Tips, FAQ & Troubleshooting
- always first perform a dry-run with option -n
- the pipeline only checks if the result directories for each method within the query region set result directories exist. Hence, if you want to re-run the analysis, you have to rename/remove these directories beforehand or use the flag --forceall.
- in case the pipeline crashes, you manually canceled your jobs or when snakemake tries to "resume.. resubmit.." jobs, then remove the .snakemake/incomplete directory!
- if you commit a lot of jobs eg via slurm (>500) this might take some time (ie 1s/job commit)
