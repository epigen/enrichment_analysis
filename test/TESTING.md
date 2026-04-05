# Test Workflow Notes

This repository includes a small integration-style test for the enrichment workflow. The goal is to exercise the normal input modes, the main analysis branches, and the reporting layer with a compact set of fixtures that can run in CI.

## What the test runs

The CI test in [.github/workflows/cy.yaml](/home/stoll/work/enrichment_analysis/.github/workflows/cy.yaml) runs:

```bash
snakemake --snakefile workflow/Snakefile --cores all --sdm conda --show-failed-logs
```

That command uses the default config loaded by [workflow/Snakefile](/home/stoll/work/enrichment_analysis/workflow/Snakefile), which is [config/config.yaml](/home/stoll/work/enrichment_analysis/config/config.yaml).

## Which config and annotation are used

The default test config points to:

- [config/config.yaml](/home/stoll/work/enrichment_analysis/config/config.yaml)
- [test/config/corces_minimal_enrichment_analysis_annotation.csv](/home/stoll/work/enrichment_analysis/test/config/corces_minimal_enrichment_analysis_annotation.csv)

The annotation file defines six feature sets:

- two region sets: `Bcell_regions`, `Ery_regions`
- two direct gene sets: `Bcell_geneSet`, `Ery_geneSet`
- two ranked gene lists: `Bcell_ranked`, `Ery_ranked`

These are grouped for summary plotting into:

- `corcesatac_regions`
- `corcesrna_gene_sets`
- `corcesrna_preranked`

## Where the test inputs come from

The fixtures under [test/data](/home/stoll/work/enrichment_analysis/test/data) are reduced copies of real MrBiomics intermediate outputs derived from the Corces example datasets used in the broader MrBiomics project, based on the hematopoietic dataset from [Corces et al. (2016)](https://www.nature.com/articles/ng.3646).

The upstream provenance is:

- `CorcesATAC`
  - corresponds to the MrBiomics ATAC-seq workflow output for the Corces hematopoietic ATAC-seq example
  - upstream flow: raw ATAC-seq reads -> `atacseq_pipeline` -> consensus-region count matrix -> `dea_limma` one-vs-all differential accessibility -> region lists -> conversion to `BED` for enrichment analysis
  - recipe documentation: [ATAC‐seq-Analysis-Recipe.md](/home/stoll/work/MrBiomics.wiki/ATAC‐seq-Analysis-Recipe.md)
- `CorcesRNA`
  - corresponds to the MrBiomics RNA-seq workflow output for the Corces hematopoietic RNA-seq example
  - upstream flow: raw RNA-seq reads -> `rnaseq_pipeline` -> gene count matrix -> `dea_limma` one-vs-all differential expression -> gene lists and ranked gene tables for enrichment analysis
  - recipe documentation: [RNA‐seq-Analysis-Recipe.md](/home/stoll/work/MrBiomics.wiki/RNA‐seq-Analysis-Recipe.md)

The current test keeps the fixture small by using two cell types, B cells and erythroid cells, while still covering every normal input class of the module.

For each lineage, the active inputs are:

- [Bcell_up_features.bed](/home/stoll/work/enrichment_analysis/test/data/CorcesATAC/Bcell_up_features.bed) and [Ery_up_features.bed](/home/stoll/work/enrichment_analysis/test/data/CorcesATAC/Ery_up_features.bed)
  - genomic consensus regions that are differentially accessible in the selected lineage
  - used as the region-set queries for `GREAT`, `LOLA`, and `pycisTarget`
- [Bcell_up_features_annot.txt](/home/stoll/work/enrichment_analysis/test/data/CorcesRNA/Bcell_up_features_annot.txt) and [Ery_up_features_annot.txt](/home/stoll/work/enrichment_analysis/test/data/CorcesRNA/Ery_up_features_annot.txt)
  - plain gene-symbol lists from the differential expression output
  - used as direct gene-set inputs for `ORA_GSEApy` and `RcisTarget`
- [Bcell_featureScores_annot.csv](/home/stoll/work/enrichment_analysis/test/data/CorcesRNA/Bcell_featureScores_annot.csv) and [Ery_featureScores_annot.csv](/home/stoll/work/enrichment_analysis/test/data/CorcesRNA/Ery_featureScores_annot.csv)
  - ranked gene tables from the differential expression output
  - used as preranked inputs for `GSEApy`

The shared background inputs are:

- [ALL_features.bed](/home/stoll/work/enrichment_analysis/test/data/CorcesATAC/ALL_features.bed)
  - background universe for the region-set analyses
- [ALL_features_annot.txt](/home/stoll/work/enrichment_analysis/test/data/CorcesRNA/ALL_features_annot.txt)
  - background universe for the direct gene-set analyses

Ranked `.csv` inputs do not use an explicit background file.

## Which input resources are used

All test resources are grouped under [test/resources](/home/stoll/work/enrichment_analysis/test/resources).

The current default config actively uses:

- [test/resources/enrichment_analysis/Azimuth_2023.json](/home/stoll/work/enrichment_analysis/test/resources/enrichment_analysis/Azimuth_2023.json)
  - a JSON dictionary of gene sets
  - used as a local database for `GREAT`, `ORA_GSEApy`, and `preranked_GSEApy`
  - converted by the workflow to GMT before the enrichment tools consume it
- [test/resources/enrichment_analysis/ReactomePathways.gmt](/home/stoll/work/enrichment_analysis/test/resources/enrichment_analysis/ReactomePathways.gmt)
  - a second local database already stored in GMT format
  - used by the same `GREAT`, `ORA_GSEApy`, and `preranked_GSEApy` branches
- [test/resources/LOLACore/hg38](/home/stoll/work/enrichment_analysis/test/resources/LOLACore/hg38)
  - local LOLA region-set database
- [test/resources/pycistarget_toy/toy_600regions.regions_vs_motifs.rankings.feather](/home/stoll/work/enrichment_analysis/test/resources/pycistarget_toy/toy_600regions.regions_vs_motifs.rankings.feather)
  - toy motif-ranking database for `pycisTarget`
- [test/resources/rcistarget_toy/toy_10genes.genes_vs_motifs.rankings.feather](/home/stoll/work/enrichment_analysis/test/resources/rcistarget_toy/toy_10genes.genes_vs_motifs.rankings.feather)
  - toy motif-ranking database for `RcisTarget`
- [test/resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl](/home/stoll/work/enrichment_analysis/test/resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl)
  - shared motif annotation table used by both `pycisTarget` and `RcisTarget`

## How the toy feather files were constructed

The two cisTarget feather files were created specifically for this test setup. They are not intended to be realistic biological databases; they are small motif-ranking matrices designed to give stable and interpretable test behavior.

Each feather is structured as:

- `pycisTarget`: motif rows x genomic region columns
- `RcisTarget`: motif rows x gene columns

Lower numeric ranks mean stronger signal for that motif.

The construction strategy was:

- choose a few motif rows that should become enriched for the selected test inputs
- place the relevant B-cell or erythroid regions or genes at the best ranks for those rows
- keep additional motif rows as background so enrichment is computed against a non-trivial motif universe

The current toy databases were tuned so that:

- [test/resources/pycistarget_toy/toy_600regions.regions_vs_motifs.rankings.feather](/home/stoll/work/enrichment_analysis/test/resources/pycistarget_toy/toy_600regions.regions_vs_motifs.rankings.feather)
  - yields non-empty `pycisTarget` results for the B-cell region set
  - yields an empty main result table for the erythroid region set at [Ery_regions_hg38_screen_v10clust.csv](/home/stoll/work/enrichment_analysis/test/results/enrichment_analysis/Ery_regions/pycisTarget/hg38_screen_v10clust/Ery_regions_hg38_screen_v10clust.csv)
  - this intentionally exercises the plotting and summary code on a branch with no enriched motifs
- [test/resources/rcistarget_toy/toy_10genes.genes_vs_motifs.rankings.feather](/home/stoll/work/enrichment_analysis/test/resources/rcistarget_toy/toy_10genes.genes_vs_motifs.rankings.feather)
  - yields non-empty `RcisTarget` results for the direct gene-set inputs
  - also yields non-empty `RcisTarget` results for the region-derived gene lists produced by `GREAT`

In other words, the toy cisTarget resources were constructed so the test covers both successful motif-enrichment branches and one known empty-output branch.

## What is tested

The current test is a minimal integration test of the normal input modes and main analysis branches of the module.

Functionally, it tests:

- region-set inputs from `.bed` files
  - region enrichment with `GREAT`
  - region enrichment with `LOLA`
  - direct region motif enrichment with `pycisTarget`
- direct gene-set inputs from `.txt` files
  - gene-set enrichment with `ORA_GSEApy`
  - gene motif enrichment with `RcisTarget`
- ranked gene-list inputs from `.csv` files
  - preranked enrichment with `GSEApy`
- region-to-gene downstream logic
  - `GREAT`-based region-to-gene association
  - reuse of the derived gene lists in downstream `ORA_GSEApy` and `RcisTarget`
- local database handling
  - JSON local database conversion to GMT before enrichment
  - direct use of an already-GMT local database
- result aggregation and reporting
  - per-feature-set result tables and plots
  - group-level summary plots
  - export of the effective config and environment files

So in terms of normal supported input classes, the test covers:

- `.bed`
- `.txt`
- `.csv`

It does not focus on malformed inputs or failure-path edge cases. It is minimal for functional workflow coverage, not for exhaustive input validation.
