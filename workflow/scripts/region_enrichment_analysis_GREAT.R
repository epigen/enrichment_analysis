
# load libraries
library("GenomicRanges")
library("rGREAT")
library("data.table")
library("rtracklayer")

# configs

#input
regions_file <- snakemake@input[["regions"]]
background_file <- snakemake@input[["background"]]
database_path <- snakemake@input[["database"]]

# output
result_path <- snakemake@output[["result"]]

# parameters
genome <- snakemake@config[["genome"]]
great_params <- snakemake@config[["great_parameters"]]
cores_n <- snakemake@threads

# set genome
if (genome=="hg19" | genome=="hg38"){
    orgdb <- "org.Hs.eg.db"
}else if(genome=="mm9" | genome=="mm10"){
    orgdb <- "org.Mm.eg.db"
}


# Handle if file is empty
if (file.info(regions_file)$size < 1){
    # create empty file to keep track of output
    f <- file(result_path, open = "a")
    close(f)

    # quit R session
    quit(status = 0, save = "no")
}

# load query and background/universe region sets (e.g., consensus region set)
regionSet_query <- import(regions_file, format = "BED")
regionSet_background <- import(background_file, format = "BED")

# load database
database = read_gmt(file.path(database_path), from = "SYMBOL", to = "ENTREZ", orgdb = orgdb)

###### GREAT

# run GREAT
res <- great(gr = regionSet_query,
      gene_sets = database,
      tss_source = genome, 
      biomart_dataset = NULL,
      min_gene_set_size = great_params[["min_gene_set_size"]], #default: 5 
      mode = great_params[["mode"]],
      basal_upstream = great_params[["basal_upstream"]],
      basal_downstream = great_params[["basal_downstream"]],
      extension = great_params[["extension"]],
      extended_tss = NULL,
      background = regionSet_background, #default: NULL
      exclude = "gap",
      cores = cores_n, #default: 1
      verbose = TRUE #default: great_opt$verbose
     )

# get & save result table
tb <- getEnrichmentTable(res, min_region_hits = 0)
tb$description <- paste(tb$description, tb$id)
fwrite(as.data.frame(tb), file=file.path(result_path), row.names=FALSE)
