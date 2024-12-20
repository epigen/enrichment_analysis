
# load libraries
library("ggplot2")
library("svglite")
library("reshape2")
library("pheatmap")
library("data.table")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source(snakemake@params[["utils_path"]])

# configs

# input
results_all_path <- snakemake@input[["results_all"]]

# output
plot_path <- snakemake@output[["summary_plot"]]
adjp_hm_path <- snakemake@output[["adjp_hm"]]
effect_hm_path <- snakemake@output[["effect_hm"]]

# parameters
tool <- snakemake@wildcards[["tool"]]
database <- snakemake@wildcards[["db"]]
group <- snakemake@wildcards[["group"]]

term_col <- snakemake@config[["column_names"]][[tool]][["term"]]
adjp_col <- snakemake@config[["column_names"]][[tool]][["adj_pvalue"]]
effect_col <- snakemake@config[["column_names"]][[tool]][["effect_size"]]

top_n <- snakemake@config[["top_terms_n"]]
adjp_cap <- snakemake@config[["adjp_cap"]]
effect_cap <- if (tool=="preranked_GSEApy") snakemake@config[["nes_cap"]] else snakemake@config[["or_cap"]]
adjp_th <- as.numeric(snakemake@config[["adjp_th"]][[tool]])
cluster_flag <- as.logical(as.numeric(snakemake@config[["cluster_summary"]]))

# stop early if results are empty
if(file.size(results_all_path) == 0L){
    file.create(plot_path)
    file.create(adjp_hm_path)
    file.create(effect_hm_path)
    quit(save = "no", status = 0)
}

# load aggregated result dataframe
results_all <- data.frame(fread(file.path(results_all_path), header=TRUE))

# stop early if results consist of only one query
if(length(unique(results_all$name))==1){
    file.create(plot_path)
    file.create(adjp_hm_path)
    file.create(effect_hm_path)
    quit(save = "no", status = 0)
}

# remove empty terms
results_all <- results_all[results_all[[term_col]]!="",]

# determine top_n most significant terms (not necessarily statistically significant!)
top_terms <- c()
for (query in unique(results_all$name)){
    tmp_result <- results_all[results_all$name==query,]
    top_n_tmp <- min(top_n, nrow(tmp_result))
    
    if(tool=="pycisTarget" | tool=="RcisTarget"){
        tmp_terms <- tmp_result[order(-tmp_result[[adjp_col]]), term_col][1:top_n_tmp]
    }else{
        tmp_terms <- tmp_result[order(tmp_result[[adjp_col]]), term_col][1:top_n_tmp]
    }
    
    top_terms <- unique(c(top_terms,tmp_terms))
}

# make adjusted p-vale and effect-size (odds-ratio or normalized enrichment scores) dataframes
adjp_df <- dcast(results_all, as.formula(paste(term_col, "~ name")), value.var = adjp_col)
rownames(adjp_df) <- adjp_df[[term_col]]
adjp_df[[term_col]] <- NULL
adjp_df <- adjp_df[,unique(results_all$name)]

effect_df <- dcast(results_all, as.formula(paste(term_col, "~ name")), value.var = effect_col)
rownames(effect_df) <- effect_df[[term_col]]
effect_df[[term_col]] <- NULL
effect_df <- effect_df[,unique(results_all$name)]

# filter by top terms
adjp_df <- adjp_df[top_terms,]
effect_df <- effect_df[top_terms,]
                                    
# fill NA for effect_df with 1 or 0 (i.e., neutral enrichment) and for adjp_df with 1 (i.e., no significance)
effect_df[is.na(effect_df)] <- if (tool=="preranked_GSEApy" | tool=="pycisTarget" | tool=="RcisTarget") 0 else 1
adjp_df[is.na(adjp_df)] <- if (tool=="pycisTarget" | tool=="RcisTarget") 0 else 1
                                     
# make stat. sign. annotation for effect-size plot later
adjp_annot <- adjp_df
if(tool=="pycisTarget" | tool=="RcisTarget"){
    adjp_annot[adjp_df >= adjp_th] <- "*"
    adjp_annot[adjp_df < adjp_th] <- ""
}else{
    adjp_annot[adjp_df <= adjp_th] <- "*"
    adjp_annot[adjp_df > adjp_th] <- ""
}

# log2 transform odds ratios
if (tool!="preranked_GSEApy" & tool!="pycisTarget" & tool!="RcisTarget"){
    effect_df <- log2(effect_df)
}

# transform and cap adjp and effect size
if(tool!="pycisTarget" & tool!="RcisTarget"){
    # cap effect_df for plotting depending on tool  abs(log2(or)) < or_cap OR abs(NES) < nes_cap
    effect_df[effect_df > effect_cap] <- effect_cap
    effect_df[effect_df < -effect_cap] <- -effect_cap

    # log10 transform adjp & cap -log10(adjpvalue) < adjp_cap
    adjp_df <- -log10(adjp_df)
    adjp_df[adjp_df > adjp_cap] <- adjp_cap
}
                                     
# plot hierarchically clustered heatmap for adjp and effect
width_hm <- 0.2 * dim(adjp_df)[2] + 5
height_hm <- 0.2 * dim(adjp_df)[1] + 3

pheatmap(adjp_df,
         display_numbers=adjp_annot,
         main= if (tool=="pycisTarget" | tool=="RcisTarget") adjp_col else"-log10(adj. p-values)",
         treeheight_row = 10,
         treeheight_col = 10,
         fontsize = 6,
         fontsize_number = 10,
         cluster_cols = cluster_flag,
         silent=TRUE,
         width=width_hm,
         height=height_hm,
         angle_col=45, 
         cellwidth=10,
         cellheight=10,
         filename=adjp_hm_path,
         breaks= if (max(adjp_df)==0) seq(0, 1, length.out=200) else seq(0, max(adjp_df), length.out=200),
         color=colorRampPalette(c("white", "red"))(200)
        )

pheatmap(effect_df,
         display_numbers=adjp_annot,
         main = if (tool=="preranked_GSEApy" | tool=="pycisTarget" | tool=="RcisTarget") effect_col else paste0("log2(",effect_col,")"),
         treeheight_row = 10,
         treeheight_col = 10,
         fontsize = 6,
         fontsize_number = 10,
         cluster_cols = cluster_flag,
         silent=TRUE,
         width=width_hm,
         height=height_hm,
         angle_col=45, 
         cellwidth=10,
         cellheight=10,
         filename=effect_hm_path,
         breaks=seq(-max(abs(effect_df)), max(abs(effect_df)), length.out=200),
         color=colorRampPalette(c("blue", "white", "red"))(200)
        )

# perform hierarchical clustering on the effect-sizes (NES or log2 odds ratios) of the terms and reorder dataframe
hc_rows <- hclust(dist(effect_df))
hc_row_names <- rownames(effect_df)[hc_rows$order]
if (cluster_flag){
    hc_cols <- hclust(dist(t(effect_df)))
    hc_col_names <- colnames(effect_df)[hc_cols$order]
    effect_df <- effect_df[hc_rows$order, hc_cols$order]
} else{
    effect_df <- effect_df[hc_rows$order,]
}

# add a column for the terms
effect_df$terms <- rownames(effect_df)
# melt data frame for plotting
plot_df <- melt(data=effect_df,
                id.vars="terms",
                measure.vars=colnames(adjp_df),
                variable.name = "feature_set",
               value.name = "effect")

# add adjusted p-values to plot dataframe
plot_df$adjp <- apply(plot_df, 1, function(x) adjp_df[x['terms'], x['feature_set']])

# set effect-size and adjusted p-value conditional to NA (odds ratios == 0 and NES == 0) for plotting
plot_df$effect[plot_df$effect==0] <- NA
plot_df$adjp[plot_df$adjp==0] <- NA

# ensure that the order of terms and feature sets is kept
plot_df$terms <- factor(plot_df$terms,levels=hc_row_names)
if (cluster_flag){
    plot_df$feature_set <- factor(plot_df$feature_set, levels=hc_col_names)
}

# stat. significance star df
if(tool=="pycisTarget" | tool=="RcisTarget"){
    adjp_star_df <- plot_df[(!is.na(plot_df$adjp)) & (plot_df$adjp >= adjp_th),]
}else{
    adjp_star_df <- plot_df[(!is.na(plot_df$adjp)) & (plot_df$adjp >= -log10(adjp_th)),]
}

# plot
enr_plot <- ggplot(plot_df, aes(x=feature_set, y=terms, fill=effect, size=adjp))+ 
        geom_point(shape=21, stroke=0.25) +
        geom_point(data = adjp_star_df, aes(x=feature_set, y=terms), shape=8, size=0.5, color = "black", alpha = 0.5) + # stars for statistical significance
        scale_fill_gradient2(midpoint=0, low="royalblue4", mid="white", high="firebrick2", space ="Lab", name = if (tool=="preranked_GSEApy" | tool=="pycisTarget" | tool=="RcisTarget") effect_col else paste0("log2(",effect_col,")")) +
        scale_y_discrete(label=addline_format) + 
        scale_size_continuous(range = c(1,5), name = if (tool=="pycisTarget" | tool=="RcisTarget") adjp_col else paste("-log10(",adjp_col,")")) +
        ggtitle(paste(tool, database, group, sep='\n')) +
        clean_theme() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
              axis.title.x=element_blank(),
              axis.title.y=element_blank()) + guides(size = guide_legend(reverse=TRUE))

width <- 0.15 * dim(adjp_df)[2] + 3
height <- 0.2 * dim(adjp_df)[1] + 2
# options(repr.plot.width=width, repr.plot.height=height)
# enr_plot
                      
# save plot
ggsave_new(filename = paste0(group,'_',database,'_summary'),
               results_path = dirname(plot_path),
               plot = enr_plot,
               width = width,
               height = height
              )
