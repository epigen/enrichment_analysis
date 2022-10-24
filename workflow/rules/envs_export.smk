# one rule per used conda environment to document the exact versions and builds of the used software        
rule env_export:
    output:
        report(os.path.join(config["result_path"],'envs','enrichment_analysis','{env}.yaml'),
                      caption="../report/software.rst", 
                      category="Software", 
                      subcategory="{}_enrichment_analysis".format(config["project_name"])
                     ),
    conda:
        "../envs/{env}.yaml"
    resources:
        mem_mb=1000, #config.get("mem", "16000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","env_{env}.log"),
    params:
        partition=config.get("partition"),
    shell:
        """
        conda env export > {output}
        """
        
# add configuration files to report        
rule config_export:
    output:
        configs = report(os.path.join(config["result_path"],'configs','enrichment_analysis','{}_config.yaml'.format(config["project_name"])), 
                         caption="../report/configs.rst", 
                         category="Configuration", 
                         subcategory="{}_enrichment_analysis".format(config["project_name"])
                        )
    resources:
        mem_mb=1000, #config.get("mem", "16000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","config_export.log"),
    params:
        partition=config.get("partition"),
    run:
        with open(output["configs"], 'w') as outfile:
            yaml.dump(config, outfile)

# export used annotation file for documentation and reproducibility         
rule annot_export:
    input:
        config["annotation"],
    output:
        annot = report(os.path.join(config["result_path"],'configs','enrichment_analysis','{}_annot.csv'.format(config["project_name"])), 
                         caption="../report/configs.rst", 
                         category="Configuration", 
                         subcategory="{}_enrichment_analysis".format(config["project_name"])
                        )
    resources:
        mem_mb=1000, #config.get("mem", "16000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","annot_export.log"),
    params:
        partition=config.get("partition"),
    shell:
        """
        cp {input} {output}
        """
