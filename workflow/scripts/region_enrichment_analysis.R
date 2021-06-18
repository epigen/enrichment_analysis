

# Load the enrichment packages
library(LOLA)
library(GenomicRanges)
library(rGREAT)
library(enrichR)
library(ggplot2)
# library(DOSE)
library(svglite)

# plot functions

addline_format <- function(x,...){
    return(gsub('(.{1,40})(\\s|$)', '\\1\n', x))
}

do_enrichment_plot <- function(plot_data, title, x, y, size, colorBy, font.size, path, filename){
    ggplot(plot_data, aes_string(x=x, y=y, size=size, color=colorBy))  +
        geom_point() +
        scale_color_continuous(low="red", high="blue", name = colorBy, guide=guide_colorbar(reverse=TRUE)) +
        ggtitle(title) + 
#         theme_dose(font.size) +
        scale_size(range=c(3, 8)) +
        scale_y_discrete(label=addline_format) +
        theme(axis.text.y=element_text(vjust=0.6))   

    ggsave(
      filename,
      plot = last_plot(),
      device = "svg",
      path = path,
      scale = 1,
      dpi = 300,
      limitsize = FALSE,
        width = 200,
        height = 250,
        units = "mm"
    )
}

# configs
regionSet_name <- snakemake@wildcards[["region_set"]] #'majority_trend_opening'
query_regions <- snakemake@input[["regions"]] # file.path('results','atac','sorted','temporal_gradients','majority_trend_opening.bed')
background_regions <- snakemake@input[["background"]] #file.path('results','atac','all','consensus_regions.bed')
# dir_results <- snakemake@output[["result_folder"]]
dir_results_LOLA <- snakemake@output[["result_LOLA"]]
dir_results_GREAT <- snakemake@output[["result_GREAT"]]
dir_results_Enrichr <- snakemake@output[["result_Enrichr"]]
genome <- "hg38"

# make directories if not exist
# dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(dir_results_LOLA), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(dir_results_GREAT), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(dir_results_Enrichr), showWarnings = FALSE, recursive = TRUE)

# load query region sets
regionSet_query = readBed(query_regions)

# load background/universe region sets (eg consensus region set)
regionSet_background = readBed(background_regions)

# needs resources downloaded from: http://big.databio.org/regiondb/LOLACoreCaches_180412.tgz
# needs simpleCache package installed
LOLACore = loadRegionDB(file.path("resources/LOLA/nm/t1/resources/regions/LOLACore",genome))
# LOLAExt = loadRegionDB("resources/LOLAExt/hg38")
jaspar_motifs = loadRegionDB(file.path("resources/LOLA/scratch/ns5bc/resources/regions/LOLAExt",genome), collections = 'jaspar_motifs')
roadmap_epigenomics = loadRegionDB(file.path("resources/LOLA/scratch/ns5bc/resources/regions/LOLAExt",genome), collections = 'roadmap_epigenomics')

# run LOLA, save results and plot sign. results for every database LOLAExt=LOLAExt)#,
dblist <- list(LOLACore=LOLACore, jaspar_motifs=jaspar_motifs,roadmap_epigenomics=roadmap_epigenomics)

for (db in names(dblist)){
    # access like that: dblist[[db]]
    print(db)
    
    # run LOLA
    locResults = runLOLA(regionSet_query, regionSet_background, dblist[[db]], cores=1)

    # determine adjusted p-values via BH method -> not needed qvalue already provided
    # locResults[['adjPvalueLog']] <- -1*log10(p.adjust(('^'(10,-1*locResults[['pValueLog']])), method = 'BH'))

    # export all results
    dir.create(file.path(dir_results_LOLA, db), showWarnings = FALSE)
    writeCombinedEnrichment(locResults, outFolder= file.path(dir_results_LOLA, db), includeSplits=TRUE)
    


    # plot top 25 significant hits ordered by mean rank

    # select top 25 significant hits by meanRnk
    plot_data <- locResults[locResults$qValue<0.05,][order(meanRnk, decreasing=FALSE),][1:25,]
    
    if (db!='LOLACore'){
        plot_data$description <- plot_data$filename
    }

    # make column description consist of unique values
    plot_data$description <- make.names(plot_data$description, unique=TRUE)
    
    # fix and order correctly for plotting
    plot_data$description <- factor(plot_data$description, levels = plot_data$description[order(plot_data$meanRnk, decreasing=TRUE)])

    # plot
    do_enrichment_plot(plot_data=plot_data, 
                       title=paste0('LOLA results of ',regionSet_name,'\nin ',db), 
                       x='oddsRatio', 
                       y='description', 
                       size='support', 
                       colorBy='qValue', 
                       font.size=10, 
                       path=file.path(dir_results_LOLA, db), 
                       filename=paste0(regionSet_name,"_",db,"_top25.svg")
                      )

    if (db=='LOLACore'){
        #  plot top 100 significant hits aggregated by cellType
        plot_data <- locResults[locResults$qValue<0.05 & complete.cases(locResults$cellType),][order(meanRnk, decreasing=FALSE),][1:100,]

        # fix and order correctly for plotting
        plot_data$cellType <- factor(plot_data$cellType, levels = levels(factor(plot_data$cellType))[order(levels(factor(plot_data$cellType)), decreasing=TRUE)])

        do_enrichment_plot(plot_data=plot_data, 
                           title=paste0('LOLA results of ',regionSet_name,'\nin ',db), 
                           x='oddsRatio', 
                           y='cellType', 
                           size='support', 
                           colorBy='qValue', 
                           font.size=10, 
                           path=file.path(dir_results_LOLA, db), 
                           filename=paste0(regionSet_name,"_",db,"_top100_cellType.svg")
                          )
    }
}

job = submitGreatJob(regionSet_query, species = genome, bg=regionSet_background)

# save job description
capture.output(job, file=file.path(dir_results_GREAT, 'GREAT_job_description.txt'), append=TRUE)

tb = getEnrichmentTables(job)

# plot & get gene-region association
pdf(file=file.path(dir_results_GREAT, paste0('GREAT_region_gene_assciations.pdf')))
res = plotRegionGeneAssociationGraphs(job)
dev.off()

# save GREAT results
for (db in names(tb)){
    print(db)
    write.table(tb[[db]], file.path(dir_results_GREAT, paste0('GREAT_',db,'.tsv')), sep='\t', row.names=FALSE)
#     print(tb[[db]][tb[[db]][,'Hyper_Adjp_BH']<0.05, c('name','Hyper_Adjp_BH')])
}

# save unique associated genes
genes <- unique(mcols(res)$gene)
write.table(genes, file.path(dir_results_GREAT, 'GREAT_genes.txt'), sep='\t', row.names=FALSE, col.names=FALSE, quote = FALSE)

# plot singificiant hits per database
for (db in names(tb)){
#     print(db)
#     print(sum(tb[[db]]$Hyper_Adjp_BH<0.05))
    
    # filter for significant hits
    plot_data <- tb[[db]][tb[[db]]$Hyper_Adjp_BH<0.05,]
    # add description column
    plot_data$description <- paste(plot_data$ID, plot_data$name, sep=' ')
    # determine ranks
    plot_data$PValue_Rnk <- rank(plot_data$Hyper_Raw_PValue)
    plot_data$Fold_Rnk <- rank(-plot_data$Hyper_Fold_Enrichment)
    plot_data$Coverage_Rnk <- rank(-plot_data$Hyper_Region_Set_Coverage)
    # calculate mean rank
    plot_data$meanRnk <- rowMeans(plot_data[,c('PValue_Rnk', 'Fold_Rnk','Coverage_Rnk')])
    #prepare data for plotting
    plot_data <- plot_data[order(plot_data$meanRnk, decreasing=FALSE),][1:min(25,dim(plot_data)[1]),]
    plot_data$description <- factor(plot_data$description, levels = plot_data$description[order(plot_data$meanRnk, decreasing=TRUE)])
    #plot
    do_enrichment_plot(plot_data=plot_data, 
                   title=paste0('GREAT results of ',regionSet_name,'\nin ',db), 
                   x='Hyper_Fold_Enrichment', 
                   y='description', 
                   size='Hyper_Region_Set_Coverage', 
                   colorBy='Hyper_Adjp_BH', 
                   font.size=10, 
                   path=file.path(dir_results_GREAT), 
                   filename=paste0(regionSet_name,"_",db,"_top25.svg")
                  )

}

setEnrichrSite("Enrichr")

dbs <- listEnrichrDbs()
dbs <- dbs[['libraryName']]

enriched <- enrichr(genes, dbs)

# save Enrichr results
for (db in names(enriched)){
    print(db)
    write.table(enriched[[db]], file.path(dir_results_Enrichr, paste0('Enrichr_',db,'.tsv')), sep='\t', row.names=FALSE)
}


# utility for overlap determination
evaltext <- function(x){
    return(eval(parse(text=x)))
}

# plot singificiant hits per database
for (db in names(enriched)){
#     print(db)
    hits <- sum(enriched[[db]]$Adjusted.P.value<0.05)
#     print(hits)
    if (hits==0){
        next
    }
    
    # filter for significant hits
    plot_data <- enriched[[db]][enriched[[db]]$Adjusted.P.value<0.05,]
    # add description column
    plot_data$description <- plot_data$Term
    # determine ranks
    plot_data$PValue_Rnk <- rank(plot_data$P.value)
    plot_data$Fold_Rnk <- rank(-plot_data$Odds.Ratio)
    plot_data$Overlap_num <- as.numeric(lapply(plot_data$Overlap, evaltext))
    plot_data$Coverage_Rnk <- rank(-plot_data$Overlap_num)
    # calculate mean rank
    plot_data$meanRnk <- rowMeans(plot_data[,c('PValue_Rnk', 'Fold_Rnk','Coverage_Rnk')])
    #prepare data for plotting
    plot_data <- plot_data[order(plot_data$meanRnk, decreasing=FALSE),][1:min(25,dim(plot_data)[1]),]
    plot_data$description <- factor(plot_data$description, levels = plot_data$description[order(plot_data$meanRnk, decreasing=TRUE)])
#     #plot
    do_enrichment_plot(plot_data=plot_data, 
                   title=paste0('Enrichr results of ',regionSet_name,'\nin ',db), 
                   x='Odds.Ratio', 
                   y='description', 
                   size='Overlap_num', 
                   colorBy='Adjusted.P.value', 
                   font.size=10, 
                   path=file.path(dir_results_Enrichr), 
                   filename=paste0(regionSet_name,"_",db,"_top25.svg")
                  )
}


