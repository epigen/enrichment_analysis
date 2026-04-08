# Test Workflow

## CI Setup

The test job in [.github/workflows/cy.yaml](../.github/workflows/cy.yaml) first runs:

```bash
bash test/setup_test_resources.sh
```

That script restores:

- [test/data/](data/) from [test/compressed_resources/test_data.tar.gz](compressed_resources/test_data.tar.gz)
- [test/resources/enrichment_analysis](resources/enrichment_analysis) from [test/compressed_resources/test_resources.tar.gz](compressed_resources/test_resources.tar.gz)
- [test/resources/LOLACore](resources/LOLACore) by downloading `LOLACoreCaches_180412.tgz` online.

After setup, CI runs the workflow through `snakemake/snakemake-github-action` using:

```bash
bash test/setup_test_resources.sh
snakemake --snakefile workflow/Snakefile --cores all --sdm conda --show-failed-logs
```

The active config and annotation are:

- [config/config.yaml](../config/config.yaml)
- [config/annotation.csv](../config/annotation.csv)

## Input Test Data

The test data are derived from MrBiomics intermediate outputs of an analysis based on the hematopoietic dataset from [Corces et al. (2016)](https://www.nature.com/articles/ng.3646). Refer to the MrBiomics [ATAC-seq recipe](../../MrBiomics.wiki/ATAC‐seq-Analysis-Recipe.md) and [RNA-seq recipe](../../MrBiomics.wiki/RNA‐seq-Analysis-Recipe.md) for more context. The dataset provides ATAC-seq and RNA-seq data for several hematopoietic cell lines. We focus only on B cells and erythroid cells here to make the test minimal. The input data are separated in two groups (see [config/annotation.csv](../config/annotation.csv)):

- `ATAC`: ATAC-seq reads were processed into a consensus-region count matrix, differential accessibility was called with `dea_limma`, and regions with increased accessibility in the selected lineage were exported as BED files for region enrichment.
- `RNA`: RNA-seq reads were processed into a gene count matrix, differential expression was called with `dea_limma`, and selected feature-score tables were kept as ranked gene-list inputs.

The active annotation defines four feature sets: `Bcell_open_regions`, `Ery_open_regions`, `Bcell_ranked`, and `Ery_ranked`.

## Test Resources

The resources restored by [test/setup_test_resources.sh](setup_test_resources.sh) all live under [test/resources/](resources/):

- `Azimuth_2023.json`: local JSON gene-set database
- `ReactomePathways.gmt`: local GMT gene-set database
- `600regions_test.regions_vs_motifs.rankings.feather`: synthetic pycisTarget ranking database
- `5000genes_test.genes_vs_motifs.rankings.feather`: synthetic RcisTarget ranking database
- `89_motifs_test.tbl`: synthetic shared motif annotation table
- `LOLACore`: LOLA database

> [!IMPORTANT]
> The cisTarget files are synthetic test resources. They were shaped to trigger stable code paths and predictable enrichment behavior, but they are not biologically meaningful and must not be used for real analysis.

## Coverage

The integration test covers:

- region-set enrichment with `GREAT` and `LOLA`
- region motif enrichment with `pycisTarget`
- ranked gene-list enrichment with `GSEApy`
- `GREAT` region-to-gene association
- downstream reuse of `GREAT`-derived gene lists in `ORA_GSEApy` and `RcisTarget`
- local JSON-to-GMT conversion and direct GMT database use
- result aggregation, plots, and exported config/env files
