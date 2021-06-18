# get functions for region dictionary
def get_region_path(wildcards):
    return regions[wildcards.region_set]['regions_bed']

def get_background_path(wildcards):
    return regions[wildcards.region_set]['background_bed']

# performs region enrichment analysis and generate plots
rule region_enrichment_analysis:
    input:
        regions=get_region_path,
        background=get_background_path,
    output:
#         result_folder=os.path.join(config['results_dir'],'{region_set}'),
        result_LOLA=report(directory(os.path.join(config['results_dir'],'{region_set}','LOLA')), patterns=["{name}.svg"], caption="../report/LOLA.rst", category="{region_set}", subcategory="LOLA"),
        result_GREAT=report(directory(os.path.join(config['results_dir'],'{region_set}','GREAT')), patterns=["{name}.svg"], caption="../report/GREAT.rst", category="{region_set}", subcategory="GREAT"),
        result_Enrichr=report(directory(os.path.join(config['results_dir'],'{region_set}','Enrichr')), patterns=["{name}.svg"], caption="../report/Enrichr.rst", category="{region_set}", subcategory="Enrichr"),
    params:
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/enrichment_analysis.yaml",
    log:
        "logs/rules/enrichment_analysis_{region_set}.log"
    script:
        "../scripts/region_enrichment_analysis.R"
