
# load libraries
library("ggplot2")
library("svglite")
library("data.table")

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# configs

# input
enrichment_result_path <- snakemake@input[["enrichment_result"]]

# output
enrichment_plot_path <- snakemake@output[["enrichment_plot"]]

# parameters
tool <- snakemake@wildcards[["tool"]]
database <- snakemake@wildcards[["db"]]
feature_set <- snakemake@wildcards[["feature_set"]] 
plot_cols <- snakemake@config[["column_names"]][[tool]]

top_n <- plot_cols[["top_n"]]
pval_col <- plot_cols[["p_value"]]
adjp_col <- plot_cols[["adj_pvalue"]]
effect_col <- plot_cols[["effect_size"]]
overlap_col <- plot_cols[["overlap"]]
term_col <- plot_cols[["term"]]

# load enrichment result
if (file.size(enrichment_result_path) != 0L){
    enrichment_result <- data.frame(fread(file.path(enrichment_result_path), header=TRUE))
}else{
    file.create(enrichment_plot_path)
    quit(save = "no", status = 0)
}

# evaluate overlap numerically if necessary
if(class(enrichment_result[[overlap_col]])=="character"){
    enrichment_result[[overlap_col]] <- as.numeric(lapply(enrichment_result[[overlap_col]], evaltext))
}

# calculate comparable effect size either NES or odds-ratio/fold based
if (tool!="preranked_GSEApy" & tool!="pycisTarget" & tool!="RcisTarget"){
    # calculate log2(effect-size) and put in new column
    effect_col_new <- paste0("log2_",effect_col)
    enrichment_result[[effect_col_new]] <- log2(enrichment_result[[effect_col]])
    effect_col <- effect_col_new
}

# determine ranks
enrichment_result$PValue_Rnk <- if (tool!="pycisTarget" & tool!="RcisTarget") rank(enrichment_result[[pval_col]]) else rank(-enrichment_result[[pval_col]])
enrichment_result$Fold_Rnk <- rank(-abs(enrichment_result[[effect_col]]))
enrichment_result$Coverage_Rnk <- rank(-enrichment_result[[overlap_col]])
# calculate and sort by mean rank
enrichment_result$meanRnk <- rowMeans(enrichment_result[,c('PValue_Rnk', 'Fold_Rnk','Coverage_Rnk')])
enrichment_result <- enrichment_result[order(enrichment_result$meanRnk, decreasing=FALSE),]

# format term column that order is kept and values are unique
enrichment_result[[term_col]] <- make.names(enrichment_result[[term_col]], unique=TRUE)
enrichment_result[[term_col]] <- factor(enrichment_result[[term_col]], levels = enrichment_result[[term_col]])

# plot top_n terms by mean_rnk
do_enrichment_plot(plot_data=enrichment_result[1:top_n,], 
               title=paste0(tool, ' results of \n',feature_set,' in ',database), 
               x=effect_col, 
               y=term_col, 
               size=overlap_col, 
               colorBy=adjp_col, 
               font.size=10, 
               path=file.path(dirname(enrichment_plot_path)), 
               filename=paste0(feature_set,"_",database),
                   top_n = top_n
              )
