# Database Download Rules

This document collects reusable Snakemake rules for downloading common databases used by the `enrichment_analysis` workflow.

The goal is to give users:

- standardized rule templates they can copy into their own Snakefile
- a clear target file layout under `resources/`
- short explanations of what each database is used for
- example output files that already match this workflow's conventions

These examples are based on the helper rules in [`workflow/rules/help.smk`](./workflow/rules/help.smk).

## General principles

The rules below follow a few conventions:

- write downloaded files into `resources/`
- write logs into `logs/wget/`
- create parent directories automatically
- use predictable output file names
- use `wget -c` for large files so interrupted downloads can resume
- keep the rule output aligned with paths that can be referenced directly in `config.yaml`

In general:

- `local_databases` expects local GMT or JSON databases
- `lola_databases` expects the extracted LOLA genome folder
- `pycistarget_parameters:databases` expects region-based ranking databases
- `rcistarget_parameters:databases` expects gene-based ranking databases
- `path_to_motif_annotations` expects the cisTarget motif annotation table

## Minimal usage pattern

Users can either:

1. copy one or more rules into their own Snakefile
2. place them into a helper file and include it with Snakemake

Example include statement:

```python
include: os.path.join("workflow", "rules", "help.smk")
```

Example target execution:

```sh
snakemake --cores 1 resources/Enrichr_databases/GO_Biological_Process_2025.gmt
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
- `resources/Enrichr_databases/Azimuth_Cell_Types_2023.gmt`

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

Sometimes users want to build their own database from simple gene lists instead of downloading a public GMT file. This is useful for:

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
        directory("resources/LOLACore/{genome}")
    params:
        url="http://big.databio.org/regiondb/LOLACoreCaches_180412.tgz",
        archive="resources/LOLACoreCaches_180412.tgz"
    threads: 1
    resources:
        mem_mb=2000
    log:
        "logs/wget/download_lola_database_{genome}.log"
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
- The useful final outputs are the extracted genome folders such as `resources/LOLACore/hg38`, not the archive itself.
- The rule is written so users can request a specific genome folder as Snakemake target, for example `resources/LOLACore/hg38`.
- Source: <https://databio.org/regiondb>

## 3. cisTarget rankings: region-based and gene-based

### What it is for

These ranking databases are used for motif enrichment with both `pycisTarget` (region-based) and `RcisTarget` (gene-based).

The two examples below each download a ranking database plus the matching motif annotation table (`motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`) used by the human `hg38` examples.

The file downloaded here is a `*.feather` ranking database from the Aerts lab cisTarget resources. In practice, this file contains precomputed motif rankings over a large collection of genomic regions for one genome assembly and one database build. `pycisTarget` uses this ranking database to test whether the motifs associated with the user's input regions are enriched compared with the background represented in the cisTarget resource.

For the example below, the downloaded file is:

- `hg38_screen_v10_clust.regions_vs_motifs.rankings.feather`

This name already tells the user a lot about the database:

- `hg38`: genome assembly
- `screen`: the region universe used by the resource
- `v10_clust`: cisTarget database version
- `regions_vs_motifs`: region-based motif ranking database
- `.feather`: binary table format used for efficient loading

### Region-based pycisTarget example

If you want to make this easier to reuse, the rule can be written with separate wildcards for species, genome, region collection, and database version.

```python
rule download_pycistarget_database:
    output:
        rankings=(
            "resources/cistarget/"
            "{genome}_{region_set}_{db_label}.regions_vs_motifs.rankings.feather"
        ),
        motif_annotations="resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    params:
        rankings_url=lambda wildcards: (
            "https://resources.aertslab.org/cistarget/databases/"
            "{species_dir}/{genome}/{region_set}/{db_dir}/region_based/"
            "{genome}_{region_set}_{db_label}.regions_vs_motifs.rankings.feather"
        ).format(
            species_dir={
                "human": "homo_sapiens",
                "mouse": "mus_musculus",
            }[wildcards.species],
            genome=wildcards.genome,
            region_set=wildcards.region_set,
            db_dir=wildcards.db_dir,
            db_label=wildcards.db_label,
        ),
        motif_annotations_url=(
            "https://resources.aertslab.org/cistarget/motif2tf/"
            "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
        )
    threads: 1
    resources:
        mem_mb=4000
    log:
        (
            "logs/wget/download_pycistarget_"
            "{species}_{genome}_{region_set}_{db_label}.log"
        )
    shell:
        """
        mkdir -p $(dirname {output.rankings}) $(dirname {output.motif_annotations}) $(dirname {log})
        wget -c -O {output.rankings} '{params.rankings_url}' > {log} 2>&1
        wget -O {output.motif_annotations} '{params.motif_annotations_url}' >> {log} 2>&1
        """
```

Example target:

```sh
snakemake --cores 1 \
  resources/cistarget/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather
```

This target corresponds to one specific official pycisTarget ranking database:

- human
- genome assembly `hg38`
- `screen` region collection
- version label `v10_clust`
- region-based motif rankings

In practice, this ranking database is usually downloaded together with the matching motif annotation table:

- `resources/cistarget/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather`
- `resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`

This pairing is useful because the ranking database tells `pycisTarget` how motifs rank across regions, while the motif annotation table is needed afterward to map enriched motifs back to transcription factors.

Other database examples that can be represented with the same modular pattern:

- `resources/cistarget/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather`
- `resources/cistarget/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.regions_vs_motifs.rankings.feather`
- `resources/cistarget/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather`

These examples illustrate how the same rule can cover:

- different genomes such as `hg38` or `mm10`
- different region universes such as `screen` or `10kbp_up_10kbp_down_full_tx`
- the same cisTarget database version, here `v10_clust`

Before using this modular rule for other combinations, users should still verify the exact folder naming on the cisTarget resource site because the directory structure must match the remote server exactly.

### Notes

- These files can be large, so resumable download is recommended.
- The exact database should match the genome assembly used in the analysis.
- Source: <https://resources.aertslab.org/cistarget/>

### Gene-based RcisTarget example

These ranking databases are used for gene-based motif enrichment with `RcisTarget`.

The file downloaded here is a `*.feather` ranking database from the Aerts lab cisTarget resources. In practice, this file contains precomputed motif rankings over a large collection of genes for one genome assembly and one database build. `RcisTarget` uses this ranking database to test whether motifs are enriched in the input gene set.

For the example below, the downloaded file is:

- `hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`

This name already tells the user a lot about the database:

- `hg38`: genome assembly
- `500bp_up_100bp_down_full_tx`: gene-centered search space used by the resource
- `v10_clust`: cisTarget database version
- `genes_vs_motifs`: gene-based motif ranking database
- `.feather`: binary table format used for efficient loading

#### Recommended rule

If you want to make this easier to reuse, the rule can be written with separate wildcards for species, genome, gene space, and database version.

```python
rule download_rcistarget_database:
    output:
        rankings=(
            "resources/cistarget/"
            "{genome}_{gene_space}_{db_label}.genes_vs_motifs.rankings.feather"
        ),
        motif_annotations="resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    params:
        rankings_url=lambda wildcards: (
            "https://resources.aertslab.org/cistarget/databases/"
            "{species_dir}/{genome}/{annotation}/{db_dir}/gene_based/"
            "{genome}_{gene_space}_{db_label}.genes_vs_motifs.rankings.feather"
        ).format(
            species_dir={
                "human": "homo_sapiens",
                "mouse": "mus_musculus",
            }[wildcards.species],
            genome=wildcards.genome,
            annotation=wildcards.annotation,
            db_dir=wildcards.db_dir,
            gene_space=wildcards.gene_space,
            db_label=wildcards.db_label,
        ),
        motif_annotations_url=(
            "https://resources.aertslab.org/cistarget/motif2tf/"
            "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
        )
    threads: 1
    resources:
        mem_mb=4000
    log:
        (
            "logs/wget/download_rcistarget_"
            "{species}_{genome}_{gene_space}_{db_label}.log"
        )
    shell:
        """
        mkdir -p $(dirname {output.rankings}) $(dirname {output.motif_annotations}) $(dirname {log})
        wget -c -O {output.rankings} '{params.rankings_url}' > {log} 2>&1
        wget -O {output.motif_annotations} '{params.motif_annotations_url}' >> {log} 2>&1
        """
```

Example target:

```sh
snakemake --cores 1 \
  resources/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
```

This target corresponds to one specific official RcisTarget ranking database:

- human
- genome assembly `hg38`
- annotation path `refseq_r80`
- gene-centered search space `500bp_up_100bp_down_full_tx`
- version label `v10_clust`
- gene-based motif rankings

In practice, this ranking database is usually downloaded together with the matching motif annotation table:

- `resources/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
- `resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`

This pairing is useful because the ranking database is used for the motif enrichment statistics, while the motif annotation table is needed to translate motif hits into TF annotations.

Other database examples that can be represented with the same modular pattern:

- `resources/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
- `resources/cistarget/hg38_10kb_up_10kb_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
- `resources/cistarget/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`

These examples illustrate how the same rule can cover:

- different genomes such as `hg38` or `mm10`
- different gene-centered search spaces such as `500bp_up_100bp_down_full_tx` or `10kb_up_10kb_down_full_tx`
- different annotation paths when needed
- the same cisTarget database version, here `v10_clust`

Before using this modular rule for other combinations, users should still verify the exact folder naming on the cisTarget resource site because the directory structure must match the remote server exactly.

### Notes

- Use a ranking file that matches the genome build and annotation model relevant to the experiment.
- As with pycisTarget, resumable download is strongly recommended.
- Source: <https://resources.aertslab.org/cistarget/>

## Example target set

The helper entry point in [`workflow/HelpSnakefile`](./workflow/HelpSnakefile) currently uses these example outputs:

- `resources/Enrichr_databases/GO_Biological_Process_2025.gmt`
- `resources/MSigDB/h.all.v2024.1.Hs.symbols.gmt`
- `resources/cistarget/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather`
- `resources/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
- `resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`

These make a good starting point for testing the helper rules.

## Suggested helper Snakefile

If users want to keep download rules separate from the main workflow, a small helper Snakefile can be used:

```python
import os


rule all:
    input:
        "resources/Enrichr_databases/GO_Biological_Process_2025.gmt",
        "resources/MSigDB/h.all.v2024.1.Hs.symbols.gmt",
        "resources/cistarget/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather",
        "resources/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",


include: os.path.join(workflow.basedir, "rules", "help.smk")
```

Example usage:

```sh
snakemake -s workflow/HelpSnakefile -n
snakemake -s workflow/HelpSnakefile --cores 1
```
