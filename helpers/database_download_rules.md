# Download Databases using Snakemake Rules

The `enrichment_analysis` workflow requires various resources for different enrichment tools (e.g., GMT files for GSEApy, `.tbl` and `.feather` files for cisTarget, region directories for LOLA). Automating their download via Snakemake ensures immediate reproducibility, version control, and seamless execution without manual file wrangling. If an analysis needs a database, Snakemake automatically fetches it on demand.

Here we share templates of Snakemake rules for downloading databases used by the `enrichment_analysis` workflow. 

## Usage

We are leveraging Snakemake's capability to automatically recognize when a required resource, specified in the configuration file, represents a missing input and triggers the appropriate download rule before continuing the downstream analysis.

1. Decide which resources you want to use for your enrichment analyses
2. Copy the corresponding rule templates into your workflow
    - We recommend a dedicated Snakefile for resources e.g. `workflow/rules/resources.smk`
3. Adapt the templates to your needs, specifically choose meaningful output paths and filenames
4. Put that exact output path into the corresponding configuration file of your enrichment analysis (see example below)

> [!IMPORTANT] 
> Make sure the output path in the configuration file matches the composed output path of the corresponding download rule.

## Single file downloads (ORA GSEA, preranked GSEA, GREAT, cisTarget)

The simplest and most robust approach for single file resources (like GMTs for over-representation analysis, or `.tbl` and `.feather` files for cisTarget) is to explicitly map the desired local filename to its direct download URL in a Python dictionary. 

### Rule template

```python
# Map the local resource basename to its exact remote download URL.
enrichment_analysis_resource_urls = {
    # Enrichr
    "KEGG_2021_Human.gmt": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human",
    # MSigDB
    "c2.cgp.v2023.2.Hs.symbols.gmt": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c2.cgp.v2023.2.Hs.symbols.gmt",
    # cisTarget motif annotation
    "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
    # cisTarget rankings
    "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather": "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather",
}

rule download_enrichment_analysis_resources:
    output:
        "resources/enrichment_analysis/{resource}"
    params:
        url=lambda wildcards: enrichment_analysis_resource_urls[wildcards.resource],
        extra="--file-allocation none --retry-wait 5 --console-log-level warn --log-level notice"
    threads: 4
    resources:
        mem_mb=1000
    log:
        "logs/aria2c/download_enrichment_analysis_resource_{resource}.log"
    wrapper:
        "v7.5.0/utils/aria2c"
```

### Corresponding configuration example

With the rule above, your configuration file points to the paths rendered by substituting `{resource}` with the dictionary keys:

```yaml
local_databases:
    KEGG_2021_Human: "resources/enrichment_analysis/KEGG_2021_Human.gmt"
    MSigDB_C1_ALL: "resources/enrichment_analysis/c2.cgp.v2023.2.Hs.symbols.gmt"
...
pycistarget_parameters:
    path_to_motif_annotations: "resources/enrichment_analysis/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    databases:
        hg38_screen_v10clust: "resources/enrichment_analysis/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
```

## Custom JSON databases from gene lists

Oftentimes, you have specific lists of genes derived from orthogonal analyses or literature that do not exist in public resources like Enrichr. You can build custom databases from simple text files, actively converting them into a `.json` database. The `enrichment_analysis` module natively consumes these `.json` databases and effortlessly converts them to `.gmt`.

This approach keeps your custom gene lists organized individually as simple text files and dynamically bundles them when required.

### Rule template

```python
from pathlib import Path
import json


rule create_custom_databases:
    input:
        Special_genes=[
            "config/genes/XYZ.txt",
            "config/genes/ABC.txt",
        ],
        Super_genes=[
            "config/genes/XYZ.txt",
            "config/genes/ABC.txt",
        ],
    output:
        Special_db="resources/custom_databases/Special.json",
        Super_db="resources/custom_databases/Super.json",
    run:
        for key, gene_paths in input.items():
            entries = {}
            for gene_file in map(Path, gene_paths):
                with gene_file.open("r", encoding="utf-8") as handle:
                    # The filename (without extension) becomes the term name
                    entries[gene_file.stem] = [
                        line.strip()
                        for line in handle
                        if line.strip()
                    ]

            destination = getattr(output, key.replace("_genes", "_db"))
            with Path(destination).open("w", encoding="utf-8") as out:
                json.dump(entries, out, indent=4)
```

### Input example

`config/genes/XYZ.txt` *(one gene symbol per line)*
```text
STAT1
IRF1
CXCL10
GBP1
```

### Output example

`Special.json`
```json
{
    "XYZ": ["STAT1", "IRF1", "CXCL10", "GBP1"],
    "ABC": ["GATA3", "IL7R", "LTB", "MAL"]
}
```

### Corresponding configuration example

With the rule above, your configuration file points to the paths of the custom databases:

```yaml
local_databases:
    Special: "resources/custom_databases/Special.json"
    Super: "resources/custom_databases/Super.json"
```

## Directory downloads (LOLA)

LOLA databases are used for genomic region set enrichment analysis. Unlike standard single file databases, LOLA requires an extracted directory structure that contains the indexed genomes (e.g. `resources/LOLACore/hg38`). 

Since different LOLA database archives (e.g. from `databio.org`) unpack into varying internal folder hierarchies, it's safer and more explicit to write a dedicated download and copy rule for the specific archive you are using.

### Rule template

This example demonstrates how to explicitly download and extract the commonly used `LOLACoreCaches_180412.tgz`.

```python
rule download_lola_database:
    output:
        directory("resources/LOLACore/hg38")
    params:
        url="http://big.databio.org/regiondb/LOLACoreCaches_180412.tgz",
        archive="resources/LOLACoreCaches_180412.tgz"
    threads: 1
    resources:
        mem_mb=2000
    conda:
        "../envs/wget.yaml" # make sure to copy this to your workflow/envs folder
    log:
        "logs/wget/download_lola_database.log"
    shell:
        """
        wget -c -O {params.archive} '{params.url}' > {log} 2>&1
        tar -xzf {params.archive} -C resources >> {log} 2>&1
        mv resources/nm/t1/resources/regions/LOLACore resources/ >> {log} 2>&1
        rm -rf resources/nm >> {log} 2>&1
        """
```

### Corresponding configuration example

Because the rule outputs directories matched to the genomes, downstream you configure the directory path. 

```yaml
lola_databases:
    LOLACore: "resources/LOLACore/hg38"
```
