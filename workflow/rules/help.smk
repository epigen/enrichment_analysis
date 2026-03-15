# Helper examples for downloading enrichment resources.
#
# This file is intentionally not included from workflow/Snakefile by default.
# If you want to make these rules available, add:
# include: os.path.join("rules", "help.smk")
#
# These are example rules users can adapt for their own project and config.
# Most rules write into the same "resources/" folder structure expected by
# config/config.yaml.
#
# Typical matching config entries:
# local_databases:
#     Hallmark: "/abs/path/to/resources/MSigDB/h.all.v2024.1.Hs.symbols.gmt"
#     GO_BP_Enrichr: "/abs/path/to/resources/Enrichr_databases/GO_Biological_Process_2025.gmt"
# lola_databases:
#     LOLACore: "/abs/path/to/resources/LOLACore/hg38"
# pycistarget_parameters:
#     databases:
#         hg38_screen_v10clust: "/abs/path/to/resources/cistarget/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
#     path_to_motif_annotations: "/abs/path/to/resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
# rcistarget_parameters:
#     databases:
#         hg38_500bp_up_100bp_down_v10clust: "/abs/path/to/resources/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
#
# Notes:
# - Enrichr libraries can be downloaded directly as GMT text.
# - MSigDB downloads may require accepting the license terms or using a logged-in session.
# - cisTarget databases are large; resumable downloads are strongly recommended.
# - LOLA expects the path to the extracted genome folder (for example:
#   resources/LOLACore/hg38), not the tarball itself.


#### Download Enrichr databases ####
rule download_Enrichr_databases:
    output:
        "resources/Enrichr_databases/{database}.gmt",
    params:
        url=lambda w: (
            "https://maayanlab.cloud/Enrichr/geneSetLibrary"
            "?mode=text&libraryName={}".format(w.database)
        ),
    resources:
        mem_mb="1000",
    threads: 1
    log:
        "logs/wget/download_Enrichr_databases_{database}.log",
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -O {output} '{params.url}' > {log} 2>&1
        """


#### Example MSigDB template ####
# Replace the URL template below with the exact collection and version you need.
# Example output already matches your current local Hallmark file name.
#
# The Broad/MSigDB release paths change over time, and some downloads may require
# an authenticated session after accepting the MSigDB terms.
rule download_MSigDB_gmt_example:
    output:
        "resources/MSigDB/{collection}.gmt",
    params:
        url=lambda w: (
            "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
            "{release}/{collection}.gmt"
        ).format(
            release=config.get("msigdb_release", "2024.1.Hs"),
            collection=w.collection,
        ),
    resources:
        mem_mb="1000",
    threads: 1
    log:
        "logs/wget/download_MSigDB_gmt_{collection}.log",
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -O {output} '{params.url}' > {log} 2>&1
        """


#### Example LOLA template ####
# Fill in the exact archive URL from https://databio.org/regiondb
# and make sure the extraction layout produces the genome folder expected by
# config["lola_databases"], e.g. resources/LOLACore/hg38
rule download_LOLA_database_example:
    output:
        directory("resources/LOLACore/{genome}"),
    params:
        url=lambda w: config.get(
            "lola_example_urls",
            {},
        ).get(
            w.genome,
            "https://databio.org/regiondb",
        ),
        archive=lambda w: os.path.join("resources", "LOLACore", "{}.tgz".format(w.genome)),
    resources:
        mem_mb="2000",
    threads: 1
    log:
        "logs/wget/download_LOLA_database_{genome}.log",
    shell:
        """
        mkdir -p resources/LOLACore $(dirname {log})
        wget -O {params.archive} '{params.url}' > {log} 2>&1
        tar -xzf {params.archive} -C resources/LOLACore >> {log} 2>&1
        test -d {output}
        """


#### Official cisTarget examples ####
# Verified against the current Aerts lab resource pages for hg38 human databases.
# These files are large; wget -c allows resuming interrupted downloads.

rule download_pycistarget_hg38_screen_v10clust:
    output:
        "resources/cistarget/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather",
    params:
        url=(
            "https://resources.aertslab.org/cistarget/databases/"
            "homo_sapiens/hg38/screen/mc_v10_clust/region_based/"
            "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
        ),
    resources:
        mem_mb="4000",
    threads: 1
    log:
        "logs/wget/download_pycistarget_hg38_screen_v10clust.log",
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -c -O {output} '{params.url}' > {log} 2>&1
        """


rule download_rcistarget_hg38_500bp_up_100bp_down_v10clust:
    output:
        "resources/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    params:
        url=(
            "https://resources.aertslab.org/cistarget/databases/"
            "homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/"
            "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
        ),
    resources:
        mem_mb="4000",
    threads: 1
    log:
        "logs/wget/download_rcistarget_hg38_500bp_up_100bp_down_v10clust.log",
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -c -O {output} '{params.url}' > {log} 2>&1
        """


rule download_cistarget_hgnc_motif_annotations:
    output:
        "resources/cistarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
    params:
        url=(
            "https://resources.aertslab.org/cistarget/motif2tf/"
            "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
        ),
    resources:
        mem_mb="1000",
    threads: 1
    log:
        "logs/wget/download_cistarget_hgnc_motif_annotations.log",
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -O {output} '{params.url}' > {log} 2>&1
        """
