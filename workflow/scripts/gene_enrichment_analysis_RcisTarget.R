
# load libraries
library("RcisTarget")
library("data.table")

process_genes <- function(input) {
  # Split the input string into individual gene names
  genes <- unlist(strsplit(input, "; "))
  
  # Remove bracketed text ending with a dot
  genes <- gsub("\\s*\\([^\\)]+\\)\\.", "", genes)
  
  # Remove duplicate gene names
  unique_genes <- unique(genes)
  
  # Join the unique gene names back into a single string
  output <- paste(unique_genes, collapse = "; ")
  
  return(output)
}

# configs

#input
genes_file <- snakemake@input[["genes"]]
background_file <- snakemake@input[["background_genes"]]
database_path <- snakemake@input[["database"]]
motif2tf_path <- snakemake@input[["motif2tf"]]

# output
result_path <- snakemake@output[["result"]]

# parameters
gene_set_name <- snakemake@wildcards[["gene_set"]]
rcistarget_params <- snakemake@config[["rcistarget_parameters"]]
cores_n <- snakemake@threads

print(rcistarget_params)

# load query and background gene sets
geneSets <- list()
geneSets[[gene_set_name]] <- readLines(genes_file)

background <- readLines(background_file)

# load database, filter for background and rer-rank
rankingsDb <- importRankings(database_path, columns = background)
motifRankings <- reRank(rankingsDb)
ranking_df <- getRanking(motifRankings)

# subset gene list for supported genes
geneSets[[gene_set_name]] <- intersect(colnames(ranking_df), geneSets[[gene_set_name]])

# load the motif to TF annotation
motifAnnot <- importAnnotations(motif2tf_path, motifsInRanking = ranking_df$features)

###### RcisTarget

# run RcisTarget
motifEnrichmentTable_wGenes <- cisTarget(geneSets = geneSets,
                                         motifRankings = motifRankings,
                                         motifAnnot = motifAnnot,
                                         motifAnnot_highConfCat = c(rcistarget_params[["motifAnnot_highConfCat"]]),
                                         motifAnnot_lowConfCat = c(rcistarget_params[["motifAnnot_lowConfCat"]]),
                                         highlightTFs = NULL,
                                         nesThreshold = rcistarget_params[["nesThreshold"]],
                                         aucMaxRank = rcistarget_params[["aucMaxRank_factor"]] * ncol(motifRankings),
                                         geneErnMethod = rcistarget_params[["geneErnMethod"]], 
                                         geneErnMaxRank = rcistarget_params[["geneErnMaxRank"]],
                                         nCores = cores_n,
                                         verbose = TRUE
                                        )

# format result table
motifEnrichmentTable_wGenes$description <- sapply(motifEnrichmentTable_wGenes$TF_highConf, process_genes)
motifEnrichmentTable_wGenes$description <- paste0(motifEnrichmentTable_wGenes$motif, " (",motifEnrichmentTable_wGenes$description,")")
# motifEnrichmentTable_wGenes$description <- paste0(motifEnrichmentTable_wGenes$motif, " (",motifEnrichmentTable_wGenes$TF_highConf,")")
motifEnrichmentTable_wGenes$name <- gene_set_name

# save result table
fwrite(as.data.frame(motifEnrichmentTable_wGenes), file=file.path(result_path), row.names=FALSE) #quote=FALSE
