
# load, convert and save local provided GMT & JSON databases to local resource folder
rule prepare_databases:
    input:
        get_db_path,
    output:
        db_file = os.path.join("resources", config["project_name"],"{database}.gmt"),
    params:
        partition = config.get("partition"),
    threads: config.get("threads", 1)
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/gene_enrichment_analysis.yaml",
    log:
        os.path.join("logs","rules","prepare_databases_{database}.log"),
    script:
        "../scripts/prepare_databases_GSEApy.py"

# # download enrichr databases to local json files using GSEApy
# rule load_enrichr_databases:
#     output:
#         result_file = os.path.join("resources", config["project_name"],"{db}.json"),
#     params:
#         database = lambda w: "{}".format(w.db),
#         partition=config.get("partition"),
#     threads: config.get("threads", 1)
#     resources:
#         mem_mb=config.get("mem", "16000"),
#     conda:
#         "../envs/gene_enrichment_analysis.yaml",
#     log:
#         os.path.join("logs","rules","get_Enrichr_databases_{db}.log"),
#     script:
#         "../scripts/get_Enrichr_databases_GSEApy.py"

# # copy local provided JSON databases to local resource folder
# rule copy_json_databases:
#     input:
#         get_json_db_path,
#     output:
#         db_file = os.path.join("resources", config["project_name"],"{db}.json"),
#     params:
#         partition=config.get("partition"),
#     threads: config.get("threads", 1)
#     resources:
#         mem_mb=1000, #config.get("mem", "16000"),
#     log:
#         os.path.join("logs","rules","copy_databases_{db}.log"),
#     shell:
#         """
#         cp {input} {output}
#         """
# # load, convert and save local provided GMT databases to local resource folder
# rule get_gmt_databases:
#     input:
#         get_gmt_db_path,
#     output:
#         result_file = os.path.join("resources", config["project_name"],"{db}.json"),
#     params:
#         database = lambda w: "{}".format(w.db),
#         partition = config.get("partition"),
#     threads: config.get("threads", 1)
#     resources:
#         mem_mb=config.get("mem", "16000"),
#     conda:
#         "../envs/gene_enrichment_analysis.yaml",
#     log:
#         os.path.join("logs","rules","get_gmt_databases_{db}.log"),
#     script:
#         "../scripts/get_gmt_databases_GSEApy.py"
        
# # download LOLA resources if needed
# rule load_lola_resources:
#     output:
#         lola_resources = directory(os.path.abspath(os.path.join("resources", config["project_name"], "LOLA"))),
#     params:
#         partition=config.get("partition"),
#     threads: config.get("threads", 1)
#     resources:
#         mem_mb=config.get("mem", "16000"),
#     log:
#         os.path.join("logs","rules","load_LOLA_resources.log"),
#     run:
#         print("start downloading and unpacking LOLA resources")
#         os.makedirs(output.lola_resources, exist_ok=True)

#         LOLA_path = output.lola_resources

#         URL_Core='http://big.databio.org/regiondb/LOLACoreCaches_180412.tgz'
#         URL_Ext='http://big.databio.org/regiondb/LOLAExtCaches_170206.tgz'

#         # download
#         getCore_str = 'wget --directory-prefix={} {}'.format(LOLA_path, URL_Core)
#         getExt_str = 'wget --directory-prefix={} {}'.format(LOLA_path, URL_Ext)

#         subprocess.run(getCore_str, shell=True)
#         subprocess.run(getExt_str, shell=True)

#         # unpack
#         unpackCore_str = 'tar zxvf {} -C {}'.format(os.path.join(LOLA_path, 'LOLACoreCaches_180412.tgz'), LOLA_path)
#         unpackExt_str = 'tar zxvf {} -C {}'.format(os.path.join(LOLA_path, 'LOLAExtCaches_170206.tgz'), LOLA_path)

#         subprocess.run(unpackCore_str, shell=True)
#         subprocess.run(unpackExt_str, shell=True)

#         # remove
#         removeCore_str = 'rm -f {}'.format(os.path.join(LOLA_path, 'LOLACoreCaches_180412.tgz'))
#         removeExt_str = 'rm -f {}'.format(os.path.join(LOLA_path, 'LOLAExtCaches_170206.tgz'))

#         subprocess.run(removeCore_str, shell=True)
#         subprocess.run(removeExt_str, shell=True)

#         print("finished downloading and unpacking LOLA resources")