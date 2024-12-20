
# aggregate all results of the same group per database
rule aggregate:
    input:
        enrichment_results = get_group_paths,
    output:
        results_all = os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_all.csv'),
        results_sig = os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_sig.csv'),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/gene_enrichment_analysis.yaml",
    log:
        "logs/rules/aggregate_{group}_{tool}_{db}.log"
    script:
        "../scripts/aggregate.py"
        
# visualize all results of the same group per database
rule visualize:
    input:
        results_all = os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_all.csv'),
    output:
        summary_plot = report(os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_summary.png'),
                             caption="../report/summary_plot.rst", 
                             category="{}_{}".format(config["project_name"], module_name),
                             subcategory="{group}",
                               labels={
                                  "name": "{tool}",
                                  "type": "summary plot",
                                  "misc": "{db}",
                              }),
        adjp_hm = os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_adjp_heatmap.pdf'),
        effect_hm = os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_effect_heatmap.pdf'),
    params:
        utils_path=workflow.source_path("../scripts/utils.R")
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/visualization.yaml",
    log:
        "logs/rules/visualize_{group}_{tool}_{db}.log"
    script:
        "../scripts/overview_plot.R"