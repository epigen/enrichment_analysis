
# performs region enrichment analysis using LOLA
rule region_enrichment_analysis_LOLA:
    input:
        regions = get_region_path,
        background = get_background_region_path,
        lola_resources = rules.load_lola_resources.output.lola_resources
#         lola_resources = os.path.abspath(os.path.join("resources", config["project_name"], "LOLA")),
    output:
        results = expand(os.path.join(result_path,'{{region_set}}','LOLA','{db}','{{region_set}}_{db}.csv'),db=config["lola_dbs"]),
    params:
        region_set = lambda w: "{}".format(w.region_set),
        result_path = result_path,
        partition=config.get("partition"),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/region_enrichment_analysis.yaml",
    log:
        "logs/rules/region_enrichment_analysis_LOLA_{region_set}.log"
    script:
        "../scripts/region_enrichment_analysis_LOLA.R"
        
# performs region enrichment analysis using GREAT
rule region_enrichment_analysis_GREAT:
    input:
        regions=get_region_path,
        background=get_background_region_path,
    output:
        results = expand(os.path.join(result_path,'{{region_set}}','GREAT','{db}','{{region_set}}_{db}.csv'),db=great_dbs),
        genes = os.path.join(result_path,'{region_set}','GREAT','genes.txt'),
    params:
        region_set = lambda w: "{}".format(w.region_set),
        partition=config.get("partition"),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/region_enrichment_analysis.yaml",
    log:
        "logs/rules/region_enrichment_analysis_GREAT_{region_set}.log"
    script:
        "../scripts/region_enrichment_analysis_GREAT.R"
                
# performs gene over-represenation analysis (ORA) using GSEApy
rule gene_ORA_GSEApy:
    input:
        query_genes=get_gene_path,
        background_genes=get_background_gene_path,
        database = os.path.join("resources", config["project_name"], "{db}.json"),
    output:
        result_file = os.path.join(result_path,'{gene_set}','ORA_GSEApy','{db}','{gene_set}_{db}.csv'),
    params:
        database = lambda w: "{}".format(w.db),
        partition=config.get("partition"),
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
        database = os.path.join("resources", config["project_name"], "{db}.json"),
    output:
        result_file = os.path.join(result_path,'{gene_set}','preranked_GSEApy','{db}','{gene_set}_{db}.csv'),
    params:
        database = lambda w: "{}".format(w.db),
        partition=config.get("partition"),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/gene_enrichment_analysis.yaml",
    log:
        "logs/rules/gene_preranked_GSEApy_{gene_set}_{db}.log"
    script:
        "../scripts/gene_preranked_GSEApy.py"
        
        
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
        tool = lambda w: "{}".format(w.tool),
        feature_set = lambda w: "{}".format(w.feature_set),
        database = lambda w: "{}".format(w.db),
        partition=config.get("partition"),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/visualization.yaml",
    log:
        "logs/rules/plot_enrichment_result_{tool}_{feature_set}_{db}.log"
    script:
        "../scripts/enrichment_plot.R"