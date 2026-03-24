# Database Download Rules

This document collects reusable Snakemake rules for downloading common databases used by the `enrichment_analysis` workflow.


## Minimal usage pattern

The easiest way to use these rules is to place them into a small helper Snakefile and run that helper file directly with Snakemake.

Recommended workflow:

1. create a file such as `workflow/DownloadDatabasesSnakefile`
2. copy the download rules you want from this document into that file
3. add a `rule all` listing the database outputs you want to create
4. run Snakemake with `-s workflow/DownloadDatabasesSnakefile`

Minimal helper Snakefile example:

```python
rule all:
    input:
        "resources/Enrichr_databases/GO_Biological_Process_2025.gmt"


rule download_enrichr_database:
    output:
        "resources/Enrichr_databases/{database}.gmt"
    params:
        url=lambda wildcards: (
            "https://maayanlab.cloud/Enrichr/geneSetLibrary"
            "?mode=text&libraryName={}".format(wildcards.database)
        )
    threads: 1
    resources:
        mem_mb=1000
    log:
        "logs/wget/download_enrichr_database_{database}.log"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -O {output} '{params.url}' > {log} 2>&1
        """
```

Run the helper file with:

```sh
snakemake -s workflow/DownloadDatabasesSnakefile --cores 1
```

Snakemake will then execute the rules needed to create the files listed in `rule all`.

If users prefer, they can also run one explicit target from the same helper file:

```sh
snakemake -s workflow/DownloadDatabasesSnakefile --cores 1 \
  resources/Enrichr_databases/GO_Biological_Process_2025.gmt
```

## 1. GMT gene set databases: Enrichr and MSigDB

### What it is for

Enrichr and MSigDB both provide GMT databases for gene set enrichment analysis. In this workflow, these files are typically used through `local_databases` for:

- `ORA_GSEApy`
- `preranked_GSEApy`
- local database queries associated with `rGREAT`

### Enrichr example

```python
rule download_enrichr_database:
    output:
        "resources/Enrichr_databases/{database}.gmt"
    params:
        url=lambda wildcards: (
            "https://maayanlab.cloud/Enrichr/geneSetLibrary"
            "?mode=text&libraryName={}".format(wildcards.database)
        )
    threads: 1
    resources:
        mem_mb=1000
    log:
        "logs/wget/download_enrichr_database_{database}.log"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -O {output} '{params.url}' > {log} 2>&1
        """
```

Example targets:

- `resources/Enrichr_databases/GO_Biological_Process_2025.gmt`
- `resources/Enrichr_databases/KEGG_2021_Human.gmt`
- `resources/Enrichr_databases/Reactome_2022.gmt`
- `resources/Enrichr_databases/Azimuth_2023.gmt`

Notes:

- The wildcard value must match the Enrichr library name exactly.
- Enrichr is one of the easiest sources for quickly obtaining GMT files.
- Source: <https://maayanlab.cloud/Enrichr/#libraries>

### MSigDB example

```python
rule download_msigdb_gmt:
    output:
        "resources/MSigDB/{collection}.gmt"
    params:
        url=lambda wildcards: (
            "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
            "2024.1.Hs/{collection}.gmt"
        ).format(collection=wildcards.collection)
    threads: 1
    resources:
        mem_mb=1000
    log:
        "logs/wget/download_msigdb_gmt_{collection}.log"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -O {output} '{params.url}' > {log} 2>&1
        """
```

Example targets:

- `resources/MSigDB/h.all.v2024.1.Hs.symbols.gmt`
- `resources/MSigDB/c2.cgp.v2024.1.Hs.symbols.gmt`

Notes:

- MSigDB release paths may change over time.
- Some MSigDB downloads may require accepting license terms or using an authenticated session.
- It is usually best to keep the full original collection name in the filename.
- Source: <https://www.gsea-msigdb.org/gsea/msigdb>

### Custom enrichment databases from gene lists

#### What it is for

If you want to build their own database from simple gene lists instead of downloading a public GMT file. This is useful for:

- project-specific marker sets
- curated literature gene sets
- pathway collections assembled inside the lab
- small custom databases for ORA or preranked GSEA

The rule below collects multiple plain-text gene lists and writes them into JSON databases that can be used as local enrichment resources.

#### Recommended rule

```python
from pathlib import Path
import json


#### create custom enrichment database from gene lists ####
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
                    entries[gene_file.stem] = [
                        line.strip()
                        for line in handle
                        if line.strip()
                    ]

            destination = getattr(output, key.replace("_genes", "_db"))
            with Path(destination).open("w", encoding="utf-8") as handle:
                json.dump(entries, handle, indent=4)
```

#### Expected input format

Each input text file should contain one gene symbol per line, for example:

```text
STAT1
IRF1
CXCL10
GBP1
```

#### Output structure

The generated JSON file will look like this:

```json
{
    "XYZ": ["STAT1", "IRF1", "CXCL10", "GBP1"],
    "ABC": ["GATA3", "IL7R", "LTB", "MAL"]
}
```

#### Notes

- The JSON keys are taken from the input filenames via `gene_file.stem`.
- Empty lines are ignored.
- Gene symbols should be standardized before building the database.
- This is a convenient way to create local databases for `local_databases`.

## 2. LOLA region databases

### What it is for

LOLA databases are used for genomic region set enrichment analysis with the `LOLA` part of the workflow.

Unlike GMT-based databases, LOLA expects a directory structure rather than a single file. The configured path should point to the extracted genome folder, for example:

- `resources/LOLACore/hg38`

The example below uses the `LOLACoreCaches_180412.tgz` archive, which is also used in the project README as a working example. This archive contains the prebuilt `LOLACore` region databases and, after extraction and moving the files into place, provides genome-specific folders such as:

- `resources/LOLACore/hg19`
- `resources/LOLACore/hg38`
- `resources/LOLACore/mm10`

These genome folders are the actual databases consumed by the workflow.


### Recommended rule

```python
rule download_lola_database:
    output:
        directory("resources/LOLACore")
    params:
        url="http://big.databio.org/regiondb/LOLACoreCaches_180412.tgz",
        archive="resources/LOLACoreCaches_180412.tgz"
    threads: 1
    resources:
        mem_mb=2000
    log:
        "logs/wget/download_lola_database.log"
    shell:
        """
        mkdir -p resources/LOLACore $(dirname {log})
        wget -O {params.archive} '{params.url}' > {log} 2>&1
        tar -xzf {params.archive} -C resources >> {log} 2>&1
        mv resources/nm/t1/resources/regions/LOLACore resources/ >> {log} 2>&1
        rm -rf resources/nm >> {log} 2>&1
        test -d {output}
        """
```

### Notes

- This example downloads the cached `LOLACore` archive from the LOLA RegionDB resources.
- The rule downloads the complete multi-genome LOLACore archive and extracts it into `resources/LOLACore/`.
- The useful final outputs are the extracted genome folders inside `resources/LOLACore/`, such as `resources/LOLACore/hg38`.
- Source: <https://databio.org/regiondb>

## 3. cisTarget rankings: region-based and gene-based

### What it is for

These ranking databases are used for motif enrichment with both `pycisTarget` (region-based) and `RcisTarget` (gene-based). They work along a motif annotation table.
The ranking database tells `pycisTarget` how motifs rank across regions, while the motif annotation table is needed afterward to map enriched motifs back to transcription factors.


If users need a different organism or annotation flavor, they should choose an alternative table from:

- <https://resources.aertslab.org/cistarget/motif2tf/>

A cisTarget ranking database is a `*.feather` file containing precomputed motif rankings for one genome assembly and one database build. In practice:

- `pycisTarget` uses region-based rankings
- `RcisTarget` uses gene-based rankings

The filenames are informative. For example:

- `hg38_screen_v10_clust.regions_vs_motifs.rankings.feather`
- `hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`

These names encode information such as:

- genome assembly such as `hg38`
- search space or region universe such as `screen` or `500bp_up_100bp_down_full_tx`
- cisTarget database version such as `v10_clust`
- whether the file is `regions_vs_motifs` or `genes_vs_motifs`

General notes:

- These files can be large, so resumable download is recommended.
- In practice, it is simplest to download the motif annotation table with a separate rule and let the ranking database rules focus only on the `*.feather` files.
- Users who need a different motif annotation table can choose one from <https://resources.aertslab.org/cistarget/motif2tf/>.
- Source: <https://resources.aertslab.org/cistarget/>

### Shared motif annotation rule

This rule downloads the tested human motif annotation table used by both example ranking databases below.
- `resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`


```python
rule download_cistarget_motif_annotations:
    output:
        "resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    params:
        url=(
            "https://resources.aertslab.org/cistarget/motif2tf/"
            "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
        )
    threads: 1
    resources:
        mem_mb=1000
    log:
        "logs/wget/download_cistarget_motif_annotations.log"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -O {output} '{params.url}' > {log} 2>&1
        """
```

### Region-based pycisTarget example

This example uses one verified human `hg38` region-based ranking database.

```python
rule download_pycistarget_database:
    input:
        motif_annotations="resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    output:
        "resources/cistarget/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
    params:
        rankings_url=(
            "https://resources.aertslab.org/cistarget/databases/"
            "homo_sapiens/hg38/screen/mc_v10_clust/region_based/"
            "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
        )
    threads: 1
    resources:
        mem_mb=4000
    log:
        "logs/wget/download_pycistarget_hg38_screen_v10_clust.log"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -c -O {output} '{params.rankings_url}' > {log} 2>&1
        """
```

Example target:

```sh
snakemake --cores 1 \
  resources/cistarget/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather
```



### Gene-based RcisTarget example

#### Recommended rule

This example uses one verified human `hg38` gene-based ranking database.

```python
rule download_rcistarget_database:
    input:
        motif_annotations="resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    output:
        "resources/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
    params:
        rankings_url=(
            "https://resources.aertslab.org/cistarget/databases/"
            "homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/"
            "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
        )
    threads: 1
    resources:
        mem_mb=4000
    log:
        "logs/wget/download_rcistarget_hg38_500bp_up_100bp_down_full_tx_v10_clust.log"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -c -O {output} '{params.rankings_url}' > {log} 2>&1
        """
```

Example target:

```sh
snakemake --cores 1 \
  resources/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
```



## Suggested helper Snakefile

If users want to keep download rules separate from the main workflow, a standalone helper Snakefile can be used. This is the pattern we tested in practice.

```python
rule all:
    input:
        "resources/Enrichr_databases/GO_Biological_Process_2025.gmt",
        "resources/MSigDB/h.all.v2024.1.Hs.symbols.gmt",
        "resources/LOLACore",
        "resources/cistarget/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather",
        "resources/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "resources/custom_databases/Special.json",
        "resources/custom_databases/Super.json"


# Then place the rule definitions from this document below, for example:
# - download_enrichr_database
# - download_msigdb_gmt
# - create_custom_databases
# - download_lola_database
# - download_cistarget_motif_annotations
# - download_pycistarget_database
# - download_rcistarget_database
```

Example usage:

```sh
snakemake -s workflow/DownloadDatabasesSnakefile -n
snakemake -s workflow/DownloadDatabasesSnakefile --cores 1
```
