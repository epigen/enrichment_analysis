
# summarize all results of the same group per database
rule summarize:
    input:
        enrichment_results = get_group_paths,
    output:
        adjpvalues = os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_adjp.csv'),
        oddsratios = os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_or.csv'),
    params:
        partition=config.get("partition"),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/gene_enrichment_analysis.yaml",
    log:
        "logs/rules/summarize_{group}_{tool}_{db}.log"
    script:
        "../scripts/summarize.py"
        
# summarize all results of the same group per database
rule visualize:
    input:
        adjpvalues = os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_adjp.csv'),
        oddsratios = os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_or.csv'),
    output:
        summary_plot = report(os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_summary.png'),
                             caption="../report/summary_plot.rst", 
                             category="{}_enrichment_analysis".format(config["project_name"]),
                             subcategory="{group}"),
        adjp_hm = report(os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_adjp_hm.pdf'),
                             caption="../report/summary_plot.rst", 
                             category="{}_enrichment_analysis".format(config["project_name"]),
                             subcategory="{group}"),
        or_hm = report(os.path.join(result_path,'{group}','{tool}','{db}','{group}_{db}_or_hm.pdf'),
                             caption="../report/summary_plot.rst", 
                             category="{}_enrichment_analysis".format(config["project_name"]),
                             subcategory="{group}"),
    params:
        partition=config.get("partition"),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/visualization.yaml",
    log:
        "logs/rules/visualize_{group}_{tool}_{db}.log"
    script:
        "../scripts/overview_plot.R"