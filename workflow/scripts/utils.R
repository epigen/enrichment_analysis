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
    
    for (format in c('png')){
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
    return(gsub('(.{1,40})(\\s|\\.|$|_)', '\\1\n', x))
}

# default plotting theme
clean_theme <- function(){
    
    # settings
    font <- "Arial"
    size <- 4
    
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