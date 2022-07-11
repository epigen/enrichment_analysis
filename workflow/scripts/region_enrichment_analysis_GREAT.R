
# load libraries
library(LOLA)
library(GenomicRanges)
library(rGREAT)

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
regionSet_query <- readBed(regions_file)

# if query and background the same, do not use a background (used to map background regions to genes for downstream gene-based analyses)
if(regions_file!=background_file){
    regionSet_background <- readBed(background_file)
}else{
    regionSet_background <- NULL
}

###### GREAT
job = submitGreatJob(regionSet_query, species = genome, bg=regionSet_background)

# save job description
capture.output(job, file=file.path(dir_result, 'job_description.txt'), append=TRUE)

# plot & get gene-region association
pdf(file=file.path(dir_result, paste0('region_gene_assciations.pdf')), width=20, height=15)
res = plotRegionGeneAssociationGraphs(job, request_interval=600)
dev.off()

# save unique associated genes
genes <- unique(mcols(res)$gene)
write(genes, file.path(gene_file))
# write.table(genes, file.path(gene_file), sep='\t', row.names=FALSE, col.names=FALSE, quote = FALSE)

# get & save enrichment results
if(regions_file!=background_file){
    for (db in databases){
        tb <- getEnrichmentTables(job,  ontology = db, category=NULL, download_by = 'tsv')
        tb$Desc <- paste(tb$Desc, tb$ID)
    #     write.table(tb[[db]], file.path(dir_result, db, paste0(region_set,'_',db,'.tsv')), sep='\t', row.names=FALSE)
        write.csv(tb[[db]], file.path(dir_result, gsub(" ", "_", db), paste0(region_set,'_',gsub(" ", "_", db),'.csv')), row.names=FALSE)
    }
}else{
    for (path in result_paths){
        file.create(path)
    }
}
