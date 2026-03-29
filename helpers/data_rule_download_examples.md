# Database Download Rules

This document collects reusable Snakemake rules and setup patterns for downloading databases used by the `enrichment_analysis` workflow.

## Intended workflow

The intended workflow is always the same:

1. decide which external resource you want to use
2. decide where you want it to live locally under `resources/`
3. put that final local path into `config/config.yaml`
4. copy the matching rule or shell snippet into `workflow/rules/resources.smk`
5. run the normal analysis workflow

The important idea is that the config stores the final local path used downstream by the analysis rules. If the file or folder is missing, Snakemake matches that path to a rule, creates it, and then continues with the analysis.


## Minimal usage pattern

The working pattern in this repo is:

- `config/config.yaml` contains local resource paths
- `workflow/rules/resources.smk` contains rules that know how to create those paths
- downstream analysis rules use those configured paths as inputs

Example:

```yaml
local_databases:
    KEGG_2021_Human: "resources/KEGG_2021_Human.gmt"
    MSigDB_H_C2_CGP: "resources/c2.cgp.v2023.2.Hs.symbols.gmt"
```

If the files above do not exist yet, Snakemake can create them through the matching download rule.

Then run the analysis normally:

```sh
snakemake --cores 4 --use-conda
```


## 1. GMT gene set databases: Enrichr and MSigDB

### What it is for

Enrichr and MSigDB both provide GMT files for gene set enrichment.

In this workflow, these resources are typically used through `local_databases` for:

- `ORA_GSEApy`
- `preranked_GSEApy`
- GMT-backed gene set enrichment through `GREAT`

The following download rule manages two common sources of GMT files:
- Enrichr: <https://maayanlab.cloud/Enrichr/#libraries>
- MSigDB: <https://data.broadinstitute.org/gsea-msigdb/msigdb/release/>

### Download and config-path logic

The user first chooses the remote GMT they want:

- for Enrichr, choose the library name
- for MSigDB, choose the exact `.gmt` filename listed in the Broad release directory

Then the user writes the desired final local path in `config.yaml`.
The rule then matches the basename and decides which remote URL to build.
As examples of what to write in the config :

```yaml
local_databases:
    KEGG_2021_Human: "resources/KEGG_2021_Human.gmt"
local_databases:
    MSigDB_C1_ALL: "esources/c2.cgp.v2023.2.Hs.symbols.gmt"
```

### Recommended rule

In the working workflow, Enrichr and MSigDB are handled by one single GMT rule:

```python
import re


MSIGDB_GMT_PATTERN = re.compile(
    r"^[^.]+\.[^.]+\.v(?P<year>\d{4})\.(?P<minor>\d+)\.(?P<species>Hs|Mm)\.(symbols|entrez)$"
)


def get_gmt_download_url(wildcards):
    database = wildcards.database
    match = MSIGDB_GMT_PATTERN.match(database)
    if match:
        release = "{year}.{minor}.{species}".format(**match.groupdict())
        return (
            "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
            "{release}/{database}.gmt"
        ).format(release=release, database=database)

    return (
        "https://maayanlab.cloud/Enrichr/geneSetLibrary"
        "?mode=text&libraryName={}".format(database)
    )


rule download_gmt_database:
    wildcard_constraints:
        database="[^/]+"
    output:
        "resources/{database}.gmt"
    params:
        url=get_gmt_download_url,
        extra="--file-allocation none --retry-wait 5 --console-log-level warn --log-level notice"
    threads: 4
    resources:
        mem_mb=1000
    log:
        "logs/wget/download_gmt_database_{database}.log"
    wrapper:
        "v7.5.0/utils/aria2c"
```

For network downloads, this example uses the Snakemake `aria2c` wrapper because it is often faster and more robust than plain `wget` or `curl`, especially for larger resource files.

### How the MSigDB vs Enrichr matching works

This is the main matching logic:

- if the basename looks like an MSigDB filename, download from MSigDB
- otherwise, treat it as an Enrichr library name

The MSigDB pattern is:

- `<collection>.<subcategory>.vYYYY.N.<Hs|Mm>.<symbols|entrez>`

Examples that match MSigDB:

- `c2.cgp.v2023.2.Hs.symbols`
- `c1.all.v2025.1.Hs.entrez`
- `h.all.v2023.2.Hs.symbols`

Examples that do not match this pattern and therefore fall through to Enrichr:

- `KEGG_2021_Human`
- `GO_Biological_Process_2025`
- `Reactome_2022`


### What users write in `config.yaml`

Examples:

Notes:

- for Enrichr, the basename must match the Enrichr library name exactly
- for MSigDB, the basename must be the real MSigDB filename without the `.gmt` stripped off
- for `GREAT` and `GSEApy`, `symbols` files are usually safer than `entrez` files unless you know your downstream ID type matches



## 2. LOLA region databases

### What it is for

LOLA databases are used for genomic region set enrichment analysis with the `LOLA` part of the workflow.

### Download and config-path logic

Unlike GMT-based databases, LOLA expects a directory structure rather than a single file. The configured path should point to the extracted genome folder, for example:

- `resources/LOLACore/hg38`

The example below uses the `LOLACoreCaches_180412.tgz` archive, which is also used in the project README as a working example. This archive contains the prebuilt `LOLACore` region databases and, after extraction and moving the files into place, provides genome-specific folders such as:

- `resources/LOLACore/hg19`
- `resources/LOLACore/hg38`
- `resources/LOLACore/mm10`

These genome folders are the actual databases consumed by the workflow. The config can look like this

```yaml
lola_databases:
    LOLACore: "resources/LOLACore/hg38"
```

### Important limitation

This rule is intentionally specific to the tested `LOLACoreCaches_180412.tgz` archive (hardcoded in the rule).

Other LOLA `tgz` files can unpack with different internal paths, so users should not assume the same `mv resources/nm/t1/resources/regions/LOLACore ...` command will work for every archive. If they want another LOLA archive, they need to inspect that tarball and manually adapt the extraction command accordingly.

Source:

- <https://databio.org/regiondb>

### Recommended rule

```python
rule download_lola_database:
    output:
        directory("resources/LOLACore/hg19")
    params:
        url="http://big.databio.org/regiondb/LOLACoreCaches_180412.tgz",
        archive="resources/LOLACoreCaches_180412.tgz"
    threads: 1
    resources:
        mem_mb=2000
    conda:
        "../envs/download.yaml"
    log:
        "logs/wget/download_lola_database.log"
    shell:
        """
        wget -O {params.archive} '{params.url}' > {log} 2>&1
        tar -xzf {params.archive} -C resources >> {log} 2>&1
        mv resources/nm/t1/resources/regions/LOLACore resources/ >> {log} 2>&1
        rm -rf resources/nm >> {log} 2>&1
        """

```




## 3. cisTarget motif annotation tables

### What it is for

The motif annotation table is a `*.tbl` file that maps motifs back to transcription factors. It is needed by both:

- `pycisTarget`
- `RcisTarget`

### Download and config-path logic

The user first chooses the exact `tbl` filename from:

- <https://resources.aertslab.org/cistarget/motif2tf/>

Then the user writes the final local path in `config.yaml`, for example:
```yaml
pycistarget_parameters:
    path_to_motif_annotations: "resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

rcistarget_parameters:
    motifAnnot: "resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
```

### Recommended rule

```python
def get_cistarget_motif2tf_url(wildcards):
    return "https://resources.aertslab.org/cistarget/motif2tf/{}".format(
        wildcards.annotation_file
    )


rule download_cistarget_motif_annotations:
    wildcard_constraints:
        annotation_file="[^/]+\\.tbl"
    output:
        "resources/cistarget/{annotation_file}"
    params:
        url=get_cistarget_motif2tf_url,
        extra="--file-allocation none --retry-wait 5 --console-log-level warn --log-level notice"
    threads: 4
    resources:
        mem_mb=1000
    log:
        "logs/wget/download_cistarget_motif_annotations_{annotation_file}.log"
    wrapper:
        "v7.5.0/utils/aria2c"
```

For cisTarget downloads, the `aria2c` wrapper is used here because these resources can be large and `aria2` is often faster and more robust than plain `wget`.

## 4. cisTarget ranking databases

### What it is for

These ranking databases are `*.feather` files used by:

- `pycisTarget` for region-based motif enrichment
- `RcisTarget` for gene-based motif enrichment


### Download and config-path logic

First, identify the correct cisTarget database path components on the AertsLab resource server: <https://resources.aertslab.org/cistarget/>

- `species`
- `genome assembly`
- `genome annotation`
- `database version`
- resource subdir : `region_based` or `gene_based`
- the actual feather `filename`

Then, encode those pieces directly in the local filename and writes the final local path in `config.yaml`.

The local filename itself needs to be structured, so that the rule can create an exact matching to the corresponding URL.

### Filename convention

The local filename is:

```text
species.genome_assembly.annotation.database_version.resource_subdir.resource_name.feather
```


Everything after that is treated as the resource-specific name.

Examples:
Examples:

```yaml
pycistarget_parameters:
    databases:
        hg38_screen_v10clust: "resources/cistarget/homo_sapiens.hg38.screen.mc_v10_clust.region_based.hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"

rcistarget_parameters:
    databases:
        hg38_500bp_up_100bp_down_v10clust: "resources/cistarget/homo_sapiens.hg38.refseq_r80.mc_v10_clust.gene_based.hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

        mm9: "resources/cistarget/mus_musculus.mm9.refseq_r45.mc9nr.gene_based.mm9-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather"

        drosophilia: "drosophila_melanogaster.dm6.flybase_r6.02.tc_v1.gene_based.encode_modERN_20190621__ChIP_seq.drosophila_melanogaster.dm6.gene_based.max.genes_vs_tracks.rankings.feather"
```

- `homo_sapiens.hg38.screen.mc_v10_clust.region_based.hg38_screen_v10_clust.regions_vs_motifs.rankings.feather`
- `homo_sapiens.hg38.refseq_r80.mc_v10_clust.gene_based.hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
- `mus_musculus.mm9.refseq_r45.mc9nr.gene_based.mm9-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather`
- 

### Recommended rule

```python
def parse_cistarget_filename(database):
    if not database.endswith(".feather"):
        raise ValueError("cisTarget filename must end with .feather: {}".format(database))

    stem = database[:-8]
    parts = stem.split(".")
    if len(parts) < 6:
        raise ValueError(
            "cisTarget filename must have at least 6 dot-separated fields: {}".format(database)
        )

    species, genome_assembly, annotation, database_version, resource_subdir = parts[:5]
    resource_name = ".".join(parts[5:])

    if resource_subdir not in {"region_based", "gene_based"}:
        raise ValueError(
            "cisTarget filename must use region_based or gene_based as the fifth field: {}".format(database)
        )

    return {
        "species": species,
        "genome_assembly": genome_assembly,
        "annotation": annotation,
        "database_version": database_version,
        "resource_subdir": resource_subdir,
        "resource_name": resource_name,
    }


def get_cistarget_download_url(wildcards):
    parsed = parse_cistarget_filename(wildcards.database)
    filename = parsed["resource_name"] + ".feather"

    return (
        "https://resources.aertslab.org/cistarget/databases/"
        "{species}/{genome_assembly}/{annotation}/{database_version}/{resource_subdir}/{filename}"
    ).format(
        species=parsed["species"],
        genome_assembly=parsed["genome_assembly"],
        annotation=parsed["annotation"],
        database_version=parsed["database_version"],
        resource_subdir=parsed["resource_subdir"],
        filename=filename,
    )


rule download_pycistarget_database:
    wildcard_constraints:
        database="[^/]+\\.feather"
    output:
        "resources/cistarget/{database}"
    params:
        url=get_cistarget_download_url,
        extra="--file-allocation none --retry-wait 5 --console-log-level warn --log-level notice"
    threads: 4
    resources:
        mem_mb=4000
    log:
        "logs/wget/download_pycistarget_{database}.log"
    wrapper:
        "v7.5.0/utils/aria2c"
```

Again, the cisTarget `feather` downloads use the `aria2c` wrapper because these files are often large enough that a parallel downloader is helpful.



## 5. Custom enrichment databases from gene lists

### What it is for

This is useful when you want to build a small custom database from your own text files rather than download a public GMT.

### Recommended rule

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
                    entries[gene_file.stem] = [
                        line.strip()
                        for line in handle
                        if line.strip()
                    ]

            destination = getattr(output, key.replace("_genes", "_db"))
            with Path(destination).open("w", encoding="utf-8") as handle:
                json.dump(entries, handle, indent=4)
```

### Expected input format

```text
STAT1
IRF1
CXCL10
GBP1
```

### Output structure

```json
{
    "XYZ": ["STAT1", "IRF1", "CXCL10", "GBP1"],
    "ABC": ["GATA3", "IL7R", "LTB", "MAL"]
}
```
