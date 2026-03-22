
# load libraries
library("GenomicRanges")
library("rGREAT")
library("data.table")
library("rtracklayer")

# functions
annot_terms_with_features <- function(res, df) {
  regions <- vector("character", nrow(df))
  genes <- vector("character", nrow(df))

  needs_annotation <- rep(TRUE, nrow(df))
  if ("p_adjust" %in% names(df)) {
    needs_annotation <- needs_annotation & (is.na(df$p_adjust) | df$p_adjust <= snakemake@config[["adjp_th"]][["GREAT"]])
  }
  

  term_ids <- unique(df$id[needs_annotation])
  annotations <- setNames(vector("list", length(term_ids)), term_ids)

  for (term in term_ids) {
    gr <- getRegionGeneAssociations(res, term_id = term)
    annotations[[term]] <- list(
      regions = paste(
        paste0(seqnames(gr), ":", pmax(0L, start(gr) - 1L), "-", end(gr)),
        collapse = ","
      ),
      annotated_genes = paste(unique(unlist(gr$annotated_genes)), collapse = ",")
    )
  }

  for (i in which(needs_annotation)) {
    annotation <- annotations[[df$id[i]]]
    if (!is.null(annotation)) {
      regions[i] <- annotation$regions
      genes[i] <- annotation$annotated_genes
    }
  }

  df$regions <- regions
  df$annotated_genes <- genes

  return(df)
}
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

# stop early for empty query or background input
if (file.size(regions_file) == 0L || file.size(background_file) == 0L){
    file.create(result_path)
    quit(save = "no", status = 0)
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
# annotate (near) significant enrichments with features
tb <- annot_terms_with_features(res, tb)
fwrite(as.data.frame(tb), file=file.path(result_path), row.names=FALSE)
