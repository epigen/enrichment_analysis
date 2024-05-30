
# load libraries
library("LOLA")
library("GenomicRanges")
library("rGREAT")
library("data.table")

# configs
regions_file <- snakemake@input[["regions"]]
background_file <- snakemake@input[["background"]]
gene_file <- snakemake@output[["genes"]]
result_paths <- snakemake@output[["results"]]

genome <- snakemake@config[["genome"]] #"hg38"
databases <- snakemake@config[["great_dbs"]]
region_set <- snakemake@params[["region_set"]]

dir_result <- dirname(gene_file)

# make result directory, if not exist
dir.create(file.path(dir_result), showWarnings = FALSE, recursive = TRUE)

# load query and background/universe region sets (eg consensus region set)
regionSet_query <- unname(readBed(regions_file))

# create empty results if query region set is too large for GREAT and quit i.e. >500,000 https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655402/File+Size
if (length(regionSet_query)>500000){
    print("query regions set is too large (>500,000)")
    for (path in result_paths){
        file.create(path)
    }
    file.create(file.path(gene_file))
    quit(save = "no", status = 0)
}

# if query and background the same, do not use a background (used to map background regions to genes for downstream gene-based analyses)
if(regions_file!=background_file){
    regionSet_background <- unname(readBed(background_file))
    # check if background region set is not too large https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655402/File+Size
    if (length(regionSet_background)>1000000){
        print("background regions set is too large (>1,000,000) and set to NULL i.e. whole genome")
        regionSet_background <- NULL
    }
}else{
    regionSet_background <- NULL
}


###### GREAT
job = submitGreatJob(regionSet_query, species = genome, bg=regionSet_background)

# save job description
capture.output(job, file=file.path(dir_result, 'job_description.txt'), append=TRUE)

# plot, get and save gene-region association
pdf(file=file.path(dir_result, paste0('region_gene_associations.pdf')), width=20, height=15)
res = plotRegionGeneAssociationGraphs(job, request_interval=600)
# write.csv(res, file.path(dir_result, 'region_gene_associations.csv'), row.names=TRUE)
fwrite(as.data.frame(res), file=file.path(dir_result, 'region_gene_associations.csv'), row.names=TRUE)
dev.off()

# save unique associated genes by using mcols(), which returns a DataFrame object containing the metadata columns.
genes <- unique(mcols(res)$gene)
write(genes, file.path(gene_file))
# write.table(genes, file.path(gene_file), sep='\t', row.names=FALSE, col.names=FALSE, quote = FALSE)

# get & save enrichment results
if(regions_file!=background_file){
    for (db in databases){
        tb <- getEnrichmentTables(job,  ontology = db, category=NULL, download_by = 'tsv')
        tb$Desc <- paste(tb$Desc, tb$ID)
    #     write.table(tb[[db]], file.path(dir_result, db, paste0(region_set,'_',db,'.tsv')), sep='\t', row.names=FALSE)
#         write.csv(tb[[db]], file.path(dir_result, gsub(" ", "_", db), paste0(region_set,'_',gsub(" ", "_", db),'.csv')), row.names=FALSE)
        fwrite(as.data.frame(tb[[db]]), file=file.path(dir_result, gsub(" ", "_", db), paste0(region_set,'_',gsub(" ", "_", db),'.csv')), row.names=FALSE)
    }
}else{
    for (path in result_paths){
        file.create(path)
    }
}
