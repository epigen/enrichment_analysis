

# Load the enrichment packages
library(LOLA)
library(GenomicRanges)

# configs
query_regions <- snakemake@input[["regions"]]
background_regions <- snakemake@input[["background"]]
genome <- snakemake@config[["genome"]] #"hg38" "mm10"
region_set <- snakemake@params[["region_set"]]
dir_result <- snakemake@params[["result_path"]]

### load data

# load query region sets
regionSet_query <- readBed(query_regions)

# load background/universe region sets (eg consensus region set)
regionSet_background <- readBed(background_regions)

# needs resources downloaded from: http://big.databio.org/regiondb/LOLACoreCaches_180412.tgz
# needs simpleCache package installed
LOLACore <- loadRegionDB(file.path("resources",snakemake@config[["project_name"]],"LOLA/nm/t1/resources/regions/LOLACore",genome))

if (genome=="hg19" | genome=="hg38"){
    jaspar_motifs <- loadRegionDB(file.path("resources",snakemake@config[["project_name"]],"LOLA/scratch/ns5bc/resources/regions/LOLAExt",genome), collections = 'jaspar_motifs')
    roadmap_epigenomics <- loadRegionDB(file.path("resources",snakemake@config[["project_name"]],"LOLA/scratch/ns5bc/resources/regions/LOLAExt",genome), collections = 'roadmap_epigenomics')
    dblist <- list(LOLACore=LOLACore, jaspar_motifs=jaspar_motifs,roadmap_epigenomics=roadmap_epigenomics)
} else{
    dblist <- list(LOLACore=LOLACore)
}

###### LOLA

# run LOLA and save results for every database
for (db in names(dblist)){
    
    # run LOLA
    locResults <- runLOLA(regionSet_query, regionSet_background, dblist[[db]], cores=1)

    # determine adjusted p-values via BH method -> not needed qvalue already provided
    # locResults[['adjPvalueLog']] <- -1*log10(p.adjust(('^'(10,-1*locResults[['pValueLog']])), method = 'BH'))

    # make description more descriptive
    if (db=='LOLACore'){
        locResults$description <- paste(locResults$description, locResults$cellType, locResults$antibody,sep='.')
    }else{
        locResults$description <- locResults$filename
    }
    
    # format description that values are unique
    locResults$description <- make.names(locResults$description, unique=TRUE)
    
    # determine raw p-value
    locResults$pValue <- ('^'(10,-1*locResults[['pValueLog']]))
    
    # export all results
    dir.create(file.path(dir_result, region_set, "LOLA", db), showWarnings = FALSE)
#     writeCombinedEnrichment(locResults, outFolder= file.path(dir_result, db), includeSplits=TRUE)
    
    write.csv(locResults, file.path(dir_result, region_set, "LOLA", db, paste0(region_set,'_',db,'.csv')), row.names=FALSE)
}    


#     if (db=='LOLACore'){
#         #  plot top 100 significant hits aggregated by cellType
#         plot_data <- locResults[locResults$qValue<0.05 & complete.cases(locResults$cellType),][order(meanRnk, decreasing=FALSE),][1:min(100,dim(plot_data)[1]),]

#         # fix and order correctly for plotting
#         plot_data$cellType <- factor(plot_data$cellType, levels = levels(factor(plot_data$cellType))[order(levels(factor(plot_data$cellType)), decreasing=TRUE)])

#         do_enrichment_plot(plot_data=plot_data, 
#                            title=paste0('LOLA results of ',regionSet_name,'\nin ',db), 
#                            x='oddsRatio', 
#                            y='cellType', 
#                            size='support', 
#                            colorBy='qValue', 
#                            font.size=10, 
#                            path=file.path(dir_results_LOLA, db), 
#                            filename=paste0(regionSet_name,"_",db,"_top100_cellType.svg")
#                           )
#     }
