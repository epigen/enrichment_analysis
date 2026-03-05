
# load libraries
# library("LOLA")
library("GenomicRanges")
library("rGREAT")
library("data.table")
library("rtracklayer")

# configs

#input
regions_file <- snakemake@input[["regions"]]
database_path <- snakemake@input[["database"]]

# output
gene_path <- snakemake@output[["genes"]]
associations_table_path <- snakemake@output[["associations_table"]]
associations_plot_path <- snakemake@output[["associations_plot"]]

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

# stop early for empty region input
if (file.size(regions_file) == 0L){
    file.create(gene_path)
    file.create(associations_table_path)
    file.create(associations_plot_path)
    quit(save = "no", status = 0)
}

# load query region set
regionSet_query <- import(regions_file, format = "BED")

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
      background = NULL,
      exclude = "gap",
      cores = cores_n, #default: 1
      verbose = TRUE #default: great_opt$verbose
     )

# plot gene-region association
pdf(file=file.path(associations_plot_path), width=12, height=4)
plotRegionGeneAssociations(res)
dev.off()

# get and save gene-region association
associations <- getRegionGeneAssociations(res)
associations_df <- as.data.frame(associations)
# BED uses 0-based starts; export start in BED-style coordinates for direct comparability to input BED.
associations_df$start <- associations_df$start - 1
fwrite(associations_df, file=file.path(associations_table_path), row.names=TRUE)

# save unique associated genes by using mcols(), which returns a DataFrame object containing the metadata columns.
genes <- unique(unlist(mcols(associations)$annotated_genes))
write(genes, file.path(gene_path))
