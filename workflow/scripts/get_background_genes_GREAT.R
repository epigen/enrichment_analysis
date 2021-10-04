

# Load the enrichment packages
library(LOLA)
library(GenomicRanges)
library(rGREAT)

background_regions <- snakemake@input[["background"]]
dir_result <- snakemake@output[["result_folder"]]
result_file <- snakemake@output[["result_file"]]
genome <- snakemake@config[["genome"]] #"hg38"


# make result directory, if not exist
dir.create(file.path(dir_result), showWarnings = FALSE, recursive = TRUE)

# load background/universe region sets (eg consensus region set)
regionSet_background = readBed(background_regions)

###### get & export background gene set for later use with Enrichr databases
job_bg = submitGreatJob(regionSet_background, species = genome)

# save job description
capture.output(job_bg, file=file.path(dir_result, 'GREAT_job_description_backgroundRegions.txt'), append=TRUE)

# plot & get gene-region association
pdf(file=file.path(dir_result, paste0('GREAT_backgroundRegion_gene_assciations.pdf')))
res_bg = plotRegionGeneAssociationGraphs(job_bg)
dev.off()

# save unique associated genes
genes <- unique(mcols(res_bg)$gene)
write.table(genes, file.path(result_file), sep='\t', row.names=FALSE, col.names=FALSE, quote = FALSE)