
# performs region enrichment analysis using LOLA
rule region_enrichment_analysis_LOLA:
    input:
        regions = get_region_path,
        background = get_background_region_path,
        database = get_lola_db_path,
    output:
        result = os.path.join(result_path,'{region_set}','LOLA','{database}','{region_set}_{database}.csv'),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/region_enrichment_analysis.yaml",
    log:
        "logs/rules/region_enrichment_analysis_LOLA_{region_set}_{database}.log"
    script:
        "../scripts/region_enrichment_analysis_LOLA.R"
        
# performs region enrichment analysis using GREAT
rule region_enrichment_analysis_GREAT:
    input:
        regions = get_region_path,
        background = get_background_region_path,
        database = os.path.join("resources", config["project_name"], "{database}.gmt"),
    output:
        result = os.path.join(result_path,'{region_set}','GREAT','{database}','{region_set}_{database}.csv'),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/region_enrichment_analysis.yaml",
    log:
        "logs/rules/region_enrichment_analysis_GREAT_{region_set}_{database}.log"
    script:
        "../scripts/region_enrichment_analysis_GREAT.R"
        
# region-gene association using GREAT for downstream gene-base analysis of genomic regions
rule region_gene_association_GREAT:
    input:
        regions = get_region_path,
        database = os.path.join("resources", config["project_name"],"{}.gmt".format(next(iter(database_dict)))), #get_first_database,
    output:
        genes = os.path.join(result_path,'{region_set}','GREAT','genes.txt'),
        associations_table = os.path.join(result_path,'{region_set}','GREAT','region_gene_associations.csv'),
        associations_plot = os.path.join(result_path,'{region_set}','GREAT','region_gene_associations.pdf'),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/region_enrichment_analysis.yaml",
    log:
        "logs/rules/region_gene_association_GREAT_{region_set}.log"
    script:
        "../scripts/region_gene_association_GREAT.R"

# performs region TFBS motif enrichment analysis using pycisTarget
rule region_motif_enrichment_analysis_pycisTarget:
    input:
        regions = get_region_path,
        ctx_db = get_pycistarget_db_path,
        motif2tf = config["pycistarget_parameters"]["path_to_motif_annotations"],
    output:
        motif_hdf5 = os.path.join(result_path,'{region_set}','pycisTarget','{database}','motif_enrichment_cistarget_{region_set}.hdf5'),
        motif_html = os.path.join(result_path,'{region_set}','pycisTarget','{database}','motif_enrichment_cistarget_{region_set}.html'),
    params:
        fraction_overlap_w_cistarget_database = config["pycistarget_parameters"]["fraction_overlap_w_cistarget_database"],
        auc_threshold = config["pycistarget_parameters"]["auc_threshold"],
        nes_threshold = config["pycistarget_parameters"]["nes_threshold"],
        rank_threshold =  config["pycistarget_parameters"]["rank_threshold"],
        annotation_version = config["pycistarget_parameters"]["annotation_version"],
        annotations_to_use = config["pycistarget_parameters"]["annotations_to_use"],
        motif_similarity_fdr = config["pycistarget_parameters"]["motif_similarity_fdr"],
        orthologous_identity_threshold = config["pycistarget_parameters"]["orthologous_identity_threshold"],
        species = 'homo_sapiens' if config["genome"] in ["hg19", "hg38"] else 'mus_musculus' if config["genome"] in ["mm9", "mm11"] else None,
    threads: 10 * config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/pycisTarget.yaml",
    log:
        "logs/rules/region_enrichment_analysis_pycisTarget_{region_set}_{database}.log"
    shell:
        """
        {{
            pycistarget cistarget \
                --cistarget_db_fname {input.ctx_db} \
                --bed_fname {input.regions} \
                --output_folder $(dirname {output.motif_hdf5}) \
                --fr_overlap_w_ctx_db {params.fraction_overlap_w_cistarget_database} \
                --auc_threshold {params.auc_threshold} \
                --nes_threshold {params.nes_threshold} \
                --rank_threshold {params.rank_threshold} \
                --path_to_motif_annotations {input.motif2tf} \
                --annotation_version {params.annotation_version} \
                --annotations_to_use {params.annotations_to_use} \
                --motif_similarity_fdr {params.motif_similarity_fdr} \
                --orthologous_identity_threshold {params.orthologous_identity_threshold} \
                --species {params.species} \
                --name {wildcards.region_set} \
                --output_mode 'hdf5' \
                --write_html
        }} || {{ 
            echo "An error occurred during the region TFBS motif enrichment analysis using pycisTarget"; touch {output.motif_hdf5} {output.motif_html}; exit 0;
        }}
        """

# postprocess results from pycisTarget
rule process_results_pycisTarget:
    input:
        motif_hdf5 = os.path.join(result_path,'{region_set}','pycisTarget','{database}','motif_enrichment_cistarget_{region_set}.hdf5'),
    output:
        motif_csv = os.path.join(result_path,'{region_set}','pycisTarget','{database}','{region_set}_{database}.csv'),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/pycisTarget.yaml",
    log:
        "logs/rules/process_results_pycisTarget_{region_set}_{database}.log"
    script:
        "../scripts/process_results_pycisTarget.py"
                
# performs gene over-represenation analysis (ORA) using GSEApy
rule gene_ORA_GSEApy:
    input:
        query_genes=get_gene_path,
        background_genes=get_background_gene_path,
        database = os.path.join("resources", config["project_name"], "{db}.gmt"),
    output:
        result_file = os.path.join(result_path,'{gene_set}','ORA_GSEApy','{db}','{gene_set}_{db}.csv'),
    params:
        database = lambda w: "{}".format(w.db),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/gene_enrichment_analysis.yaml",
    log:
        "logs/rules/gene_ORA_GSEApy_{gene_set}_{db}.log"
    script:
        "../scripts/gene_ORA_GSEApy.py"
        
# performs gene preranked GSEA and generate plots using GSEApy
rule gene_preranked_GSEApy:
    input:
        query_genes=get_rnk_path,
        database = os.path.join("resources", config["project_name"], "{db}.gmt"),
    output:
        result_file = os.path.join(result_path,'{gene_set}','preranked_GSEApy','{db}','{gene_set}_{db}.csv'),
    params:
        database = lambda w: "{}".format(w.db),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/gene_enrichment_analysis.yaml",
    log:
        "logs/rules/gene_preranked_GSEApy_{gene_set}_{db}.log"
    script:
        "../scripts/gene_preranked_GSEApy.py"

# performs TFBS motif enrichment analysis on genes using RcisTarget
rule gene_motif_enrichment_analysis_RcisTarget:
    input:
        genes=get_gene_path,
        background_genes=get_background_gene_path,
        database = get_rcistarget_db_path,
        motif2tf = config["rcistarget_parameters"]["motifAnnot"],
    output:
        result = os.path.join(result_path,'{gene_set}','RcisTarget','{database}','{gene_set}_{database}.csv'),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/RcisTarget.yaml",
    log:
        "logs/rules/gene_motif_enrichment_analysis_RcisTarget_{gene_set}_{database}.log"
    script:
        "../scripts/gene_enrichment_analysis_RcisTarget.R"

# plot enrichment results
rule plot_enrichment_result:
    input:
        enrichment_result=os.path.join(result_path,'{feature_set}','{tool}','{db}','{feature_set}_{db}.csv'),
    output:
        enrichment_plot=report(os.path.join(result_path,'{feature_set}','{tool}','{db}','{feature_set}_{db}.png'),
                             caption="../report/enrichment_plot.rst", 
                             category="{}_{}".format(config["project_name"], module_name),
                             subcategory="{feature_set}",
                               labels={
                                  "name": "{tool}",
                                  "type": "enrichment plot",
                                  "misc": "{db}",
                              }),
    params:
        utils_path=workflow.source_path("../scripts/utils.R")
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/visualization.yaml",
    log:
        "logs/rules/plot_enrichment_result_{tool}_{feature_set}_{db}.log"
    script:
        "../scripts/enrichment_plot.R"
