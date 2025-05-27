
# load, convert and save local provided GMT & JSON databases to local resource folder
rule prepare_databases:
    input:
        get_db_path,
    output:
        db_file = os.path.join("resources",config["project_name"],module_name,"{database}.gmt"),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/gene_enrichment_analysis.yaml",
    log:
        os.path.join("logs","rules","prepare_databases_{database}.log"),
    script:
        "../scripts/prepare_databases_GSEApy.py"
