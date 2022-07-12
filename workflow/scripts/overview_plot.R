
# load libraries
library(ggplot2)
library(svglite)
library(reshape2)
library(pheatmap)

# source utility functions
source("workflow/scripts/utils.R")

# configs
adjp_path <- snakemake@input[["adjpvalues"]] #"/research/home/sreichl/projects/genomic_region_enrichment/test/enrichment_analysis/testgroup/GSEApy/GO_Biological_Process_2021/testgroup_GO_Biological_Process_2021_adjp.csv"
or_path <- snakemake@input[["oddsratios"]] #"/research/home/sreichl/projects/genomic_region_enrichment/test/enrichment_analysis/testgroup/GSEApy/GO_Biological_Process_2021/testgroup_GO_Biological_Process_2021_or.csv"

plot_path <- snakemake@output[["summary_plot"]] #"/research/home/sreichl/projects/genomic_region_enrichment/test/enrichment_analysis/testgroup/GSEApy/GO_Biological_Process_2021/testgroup_GO_Biological_Process_2021_summary.png"
adjp_hm_path <- snakemake@output[["adjp_hm"]] #"/research/home/sreichl/projects/genomic_region_enrichment/test/enrichment_analysis/testgroup/GSEApy/GO_Biological_Process_2021/testgroup_GO_Biological_Process_2021_adjp_hm.pdf"
or_hm_path <- snakemake@output[["or_hm"]] #"/research/home/sreichl/projects/genomic_region_enrichment/test/enrichment_analysis/testgroup/GSEApy/GO_Biological_Process_2021/testgroup_GO_Biological_Process_2021_or_hm.pdf"

tool <- snakemake@wildcards[["tool"]] #"GSEApy"
database <- snakemake@wildcards[["db"]] #"GO_Biological_Process_2021"
group <- snakemake@wildcards[["group"]] #"testgroup"

top_n <- snakemake@config[["top_terms_n"]] #5
adjp_cap <- snakemake@config[["adjp_cap"]] #4

# load summaries
adjp_df <- read.csv(adjp_path, row.names = 1, header= TRUE)
or_df <- read.csv(or_path, row.names = 1, header= TRUE)

# stop early if results are empty
if(dim(adjp_df)[1]<2){
    file.create(plot_path)
    quit()
}

# determine top_n most significant terms (not necessarily statistically significant!)
top_terms <- unique(as.vector(sapply(adjp_df, function(x) rownames(adjp_df)[order(x)][1:top_n])))

# filter by top terms
adjp_df <- adjp_df[top_terms,]
or_df <- or_df[top_terms,]
                                    
# fill or NA or <1 with 1 (ie neutral or negative enrichment)
or_df[is.na(or_df)] <- 1
or_df[or_df<1] <- 1
adjp_df[is.na(adjp_df)] <- 1

# log2 transform odds ratios
or_df <- log2(or_df)

# log10 transform adjp: calculate & cap -log10(adjpvalue) < adjp_cap
adjp_df <- -log10(adjp_df)
adjp_df[adjp_df>adjp_cap] <- adjp_cap
                                     
# plot hierarchically clustered heatmap for adjp and or
width_hm <- 0.5 * dim(adjp_df)[2] + 4
height_hm <- 0.2 * dim(adjp_df)[1] + 1

pheatmap(adjp_df,
         main="-log10(adj. p-values)",
         treeheight_row = 10,
         treeheight_col = 10,
         fontsize = 6,
         silent=TRUE,
        width=width_hm,
        height=height_hm,
       angle_col=45, 
         cellwidth=10,
         cellheight=10,
         filename=adjp_hm_path)
                                     
pheatmap(or_df,
         main="log2(odds ratios)",
         treeheight_row = 10,
         treeheight_col = 10,
         fontsize = 6,
         silent=TRUE,
        width=width_hm,
        height=height_hm,
       angle_col=45, 
         cellwidth=10,
         cellheight=10,
         filename=or_hm_path)


# perform hierarchical clustering on the log2 odds ratios of the terms and reorder DF
hc_rows <- hclust(dist(or_df))
hc_row_names <- rownames(or_df)[hc_rows$order]
hc_cols <- hclust(dist(t(or_df)))
hc_col_names <- colnames(or_df)[hc_cols$order]
or_df <- or_df[hc_rows$order, hc_cols$order]

# add a column for the terms
or_df$terms <- rownames(or_df)
# melt data frame for plotting
plot_df <- melt(data=or_df,
                id.vars="terms",
                measure.vars=colnames(adjp_df),
                variable.name = "feature_set",
               value.name = "or")

# add adjusted p-values to plot dataframe
plot_df$adjp <- apply(plot_df, 1, function(x) adjp_df[x['terms'], x['feature_set']])

# set odds ratios <=0 to NA for plotting
plot_df$or[plot_df$or<=0] <- NA

# ensure that the order of terms and feature sets is kept
plot_df$terms <- factor(plot_df$terms,levels=rev(unique(plot_df$terms)))
plot_df$feature_set <- factor(plot_df$feature_set, levels=rev(unique(plot_df$feature_set)))

# plot
enr_plot <- ggplot(plot_df, aes(x=feature_set, y=terms, fill=adjp, size=or))+ 
geom_point(shape=21, stroke=0.25) +
scale_fill_gradient(low="grey", high="red", breaks = c(1, 2, 3, 4), limits = c(0, 4), name="-log10(adjp)") +
scale_y_discrete(label=addline_format) + 
scale_size_continuous(range = c(0.5,5), , name="log2(or)") +
ggtitle(paste(tool, database, group, sep='\n'))+
clean_theme()+
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
      axis.title.x=element_blank(),
      axis.title.y=element_blank())

width <- 0.5 * dim(adjp_df)[2] + 3
height <- 0.5 * length(top_terms) + 2
# options(repr.plot.width=width, repr.plot.height=height)
# enr_plot
                      
# save plot
ggsave_new(filename=paste0(group,'_',database,'_summary'),
               results_path=dirname(plot_path),
               plot=enr_plot,
               width=width,
               height=height
              )
