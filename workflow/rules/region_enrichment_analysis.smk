# get functions for region dictionary
def get_region_path(wildcards):
    return regions[wildcards.region_set]['regions_bed']

# for region enrichment
def get_background_path(wildcards):
    return regions[wildcards.region_set]['background_bed']

# for region to gene mapping of background region sets
def get_background_file(wildcards):
    return background_regions[wildcards.background_regions]['background_bed']

def get_background_genes(wildcards):
    return os.path.join(config['results_dir'],'background_genes', regions[wildcards.region_set]['background_name'],'GREAT_background_genes.txt')

# performs region enrichment analysis and generate plots w/ LOLA, GREAT & enrichR
rule region_enrichment_analysis:
    input:
        regions=get_region_path,
        background=get_background_path,
    output:
#         result_folder=os.path.join(config['results_dir'],'{region_set}'),
        result_LOLA=report(directory(os.path.join(config['results_dir'],'{region_set}','LOLA')), patterns=["{name}.svg"], caption="../report/LOLA.rst", category="{region_set}", subcategory="LOLA"),
        result_GREAT=report(directory(os.path.join(config['results_dir'],'{region_set}','GREAT')), patterns=["{name}.svg"], caption="../report/GREAT.rst", category="{region_set}", subcategory="GREAT"),
        result_Enrichr=report(directory(os.path.join(config['results_dir'],'{region_set}','Enrichr')), patterns=["{name}.svg"], caption="../report/Enrichr.rst", category="{region_set}", subcategory="Enrichr"),
        GREAT_genes=os.path.join(config['results_dir'],'{region_set}','GREAT','GREAT_genes.txt'),
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
        
        
        
# get gene mapping for all background region sets
rule get_background_genes:
    input:
        background=get_background_file,
    output:
        result_folder=os.path.join(config['results_dir'],'background_genes','{background_regions}'),
        result_file=os.path.join(config['results_dir'],'background_genes','{background_regions}', 'GREAT_background_genes.txt')
    params:
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/enrichment_analysis.yaml",
    log:
        "logs/rules/get_background_genes_{background_regions}.log"
    script:
        "../scripts/get_background_genes_GREAT.R"
        

# performs region enrichment analysis and generate plots w/ GSEApy
rule GSEApy_enrichment_analysis:
    input:
        query_genes=os.path.join(config['results_dir'],'{region_set}','GREAT','GREAT_genes.txt'),
        background_genes=get_background_genes,
        enrichr_databases = os.path.join("resources", "enrichr_databases.pkl"),
    output:
        result_GSEApy=report(directory(os.path.join(config['results_dir'],'{region_set}','GSEApy')), patterns=["{name}.svg"], caption="../report/GSEApy.rst", category="{region_set}", subcategory="GSEApy"),
    params:
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/enrichment_analysis.yaml",
    log:
        "logs/rules/GSEApy_enrichment_analysis_{region_set}.log"
    script:
        "../scripts/Enrichr_enrichment_GSEApy.py"
