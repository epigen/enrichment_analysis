
# load libraries
library("LOLA")
library("GenomicRanges")
library("data.table")

# configs

# input
query_regions <- snakemake@input[["regions"]]
background_regions <- snakemake@input[["background"]]
database_path <- snakemake@input[["database"]]

# output
result_path <- snakemake@output[["result"]]

# parameters
database_name <- snakemake@wildcards[["database"]]
genome <- snakemake@config[["genome"]] #"hg38" "mm10"
region_set <- snakemake@wildcards[["region_set"]]

### load data

# load query region sets
regionSet_query <- readBed(query_regions)

# load background/universe region sets (e.g., consensus region set)
regionSet_background <- readBed(background_regions)

# requires resources downloaded from: https://databio.org/regiondb
# requires simpleCache package installed
database <- loadRegionDB(file.path(database_path))

###### LOLA

# run LOLA
res <- runLOLA(regionSet_query, regionSet_background, database, cores=1)

# make description more descriptive
if (database_name=='LOLACore'){
    res$description <- paste(res$description, res$cellType, res$antibody, sep='.')
}else{
    res$description <- paste(res$description, res$filename, sep='.')
}
    
# format description that values are unique
res$description <- make.names(res$description, unique=TRUE)

# determine raw p-value
res$pValue <- ('^'(10,-1*res[['pValueLog']]))

# save results
fwrite(as.data.frame(res), file=file.path(result_path), row.names=FALSE)