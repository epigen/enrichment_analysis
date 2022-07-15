### Utility Functions

# numeric overlap determination
evaltext <- function(x){
    return(eval(parse(text=x)))
}

# extended ggsave
ggsave_new <- function(filename, results_path, plot, width=5, height=5, units="in"){
    # make result directory if not exist
    if (!dir.exists(results_path)){
        dir.create(results_path, recursive = TRUE)
    }
    
    for (format in c('png','svg')){
        ggsave(
          paste0(filename,'.',format),
          plot = plot,
          device = format,
          path = file.path(results_path),
          scale = 1,
          dpi = 300,
            width = width,
            height = height,
          limitsize = FALSE,
            units=units
        )
    }
}

# plot functions

# add new lines in text of plots
addline_format <- function(x,...){
    return(gsub('(.{1,40})(\\s|\\.|$)', '\\1\n', x))
}

# draw enrichment plot
do_enrichment_plot <- function(plot_data, title, x, y, size, colorBy, font.size, path, filename, top_n){
    enr_p <- ggplot(plot_data, aes_string(x=x, y=y, size=size, color=colorBy))  +
        geom_point() +
        scale_color_continuous(low="red", high="blue", name = colorBy, guide=guide_colorbar(reverse=TRUE)) +
        ggtitle(title) + 
#         theme_dose(font.size) +
        scale_size(range=c(3, 8)) +
        scale_y_discrete(label=addline_format, limits=rev) +
        theme(axis.text.y=element_text(vjust=0.6))
    
    ggsave_new(filename = filename, 
           results_path=path, 
           plot=enr_p, 
           width=200, 
           height=10*top_n,
              units = "mm")
}

# default plotting theme
clean_theme <- function(){
    
    # settings
    font <- "Arial"
    size <- 6
    
    theme_bw(
        base_size=size,
        base_family = font
    ) %+replace% 
    
    theme(
      #grid elements
#       panel.grid.major = element_blank(),    #strip major gridlines
#       panel.grid.minor = element_blank(),    #strip minor gridlines
#       axis.ticks = element_blank(),          #strip axis ticks
      
      #text elements
      plot.title = element_text(             #title
                   family = font,            #set font family
                   size = size,                #set font size
                   face = 'bold',            #bold typeface
                   hjust = 0.5,                #center align
                   vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
                   family = font,            #font family
                   size = size),               #font size
      
      plot.caption = element_text(           #caption
                   family = font,            #font family
                   size = size,                 #font size
                   hjust = 0.5),               #center align
      
      axis.title = element_text(             #axis titles
                   family = font,            #font family
                   size = size),               #font size
      
      axis.text = element_text(              #axis text
                   family = font,            #axis famuly
                   size = size),                #font size
      
#       axis.text.x = element_text(            #margin for axis text
#                     margin=margin(5, b = 10))
    )
}