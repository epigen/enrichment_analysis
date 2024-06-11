### utility GET/INPUT functions

# get user provided local database path
def get_db_path(wildcards):
    return database_dict[wildcards.database]

# get user provided LOLA database path
def get_lola_db_path(wildcards):
    return lola_db_dict[wildcards.database]

# get user provided pycisTarget database path
def get_pycistarget_db_path(wildcards):
    return pycistarget_db_dict[wildcards.database]

# get user provided RcisTarget database path
def get_rcistarget_db_path(wildcards):
    return rcistarget_db_dict[wildcards.database]

### for genomic region enrichment
# region set
def get_region_path(wildcards):
    if wildcards.region_set in regions_dict.keys():
        return regions_dict[wildcards.region_set]['features_path']
    elif wildcards.region_set in background_regions_dict.keys():
        return background_regions_dict[wildcards.region_set]['background_path']
    else:
        print("Region set not found")

# background region set
def get_background_region_path(wildcards):
    if wildcards.region_set in regions_dict.keys():
        return regions_dict[wildcards.region_set]['background_path']
    elif wildcards.region_set in background_regions_dict.keys():
        return background_regions_dict[wildcards.region_set]['background_path']
    else:
        print("Background region set not found")

### for ORA GSEA
# gene set
def get_gene_path(wildcards):
    if wildcards.gene_set in genes_dict.keys():
        return os.path.join(genes_dict[wildcards.gene_set]['features_path'])
    elif wildcards.gene_set in regions_dict.keys():
        return os.path.join(result_path, wildcards.gene_set,'GREAT','genes.txt')
    else:
        print("Gene set not found")
    
# background gene set
def get_background_gene_path(wildcards):
    if wildcards.gene_set in genes_dict.keys():
        return os.path.join(genes_dict[wildcards.gene_set]['background_path'])
    elif wildcards.gene_set in regions_dict.keys():
        return os.path.join(result_path, regions_dict[wildcards.gene_set]['background_name'],'GREAT','genes.txt')
    else:
        print("Background gene set not found")

### for preranked GSEA
def get_rnk_path(wildcards):
    return os.path.join(rnk_dict[wildcards.gene_set]['features_path'])

### for group summary & visualization
def get_group_paths(wildcards):
    feature_sets = list(annot.index[annot["group"]==wildcards.group])
    
    # for tool GREAT, LOLA or pycisTarget only consider region sets
    if wildcards.tool=="GREAT" or wildcards.tool=="LOLA" or wildcards.tool=="pycisTarget":
        feature_sets = [feature_set for feature_set in feature_sets if feature_set in regions_dict.keys()]
    if wildcards.tool=="ORA_GSEApy"or wildcards.tool=="RcisTarget":
        feature_sets = [feature_set for feature_set in feature_sets if feature_set in genes_dict.keys()] + [feature_set for feature_set in feature_sets if feature_set in regions_dict.keys()]
    if wildcards.tool=="preranked_GSEApy":
        feature_sets = [feature_set for feature_set in feature_sets if feature_set in rnk_dict.keys()]
    
    return expand(os.path.join(result_path,'{feature_set}','{tool}','{db}','{feature_set}_{db}.csv'),
                  feature_set=feature_sets,
                 tool=wildcards.tool,
                 db=wildcards.db,)
