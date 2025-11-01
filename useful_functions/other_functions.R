### my_ggMarginal function
### inspired by https://stackoverflow.com/questions/8545035/scatterplot-with-marginal-histograms-in-ggplot2
my_ggMarginal <- function(
        df, 
        x_var=NULL, y_var=NULL,
        color_var=NULL,
        my_xlim=NULL, my_ylim=NULL, 
        my_color_scheme=NULL, 
        my_title=NULL, my_subtitle=NULL,
        combine_plots=TRUE, # if FALSE we return a list of plots not the patchwork
        ## the ... will go in geom_point
        ...) {
    
    
    ## checks:
    if(is.null(x_var) | is.null(y_var)) {
        stop("\n\nYou must supply x_var and y_var\n\n")
    }
    if(!x_var %in% colnames(df)) {
        stop("\n\nYour x variable name ", x_var, " is not in the data frame you supplied\n\n")
    }
    if(!y_var %in% colnames(df)) {
        stop("\n\nYour y variable name ", y_var, " is not in the data frame you supplied\n\n")
    }
    if(!color_var %in% colnames(df)) {
        stop("\n\nYour colors variable name ", color_var, " is not in the data frame you supplied\n\n")
    }
    
    ## set up the plot aes without colors for groups
    if (is.null(color_var)) {
        p1 <- df %>%
            ggplot(aes(x=.data[[x_var]], y=.data[[y_var]])) 
        ## make the density plots
        dens_x <- df %>%  
            ggplot(aes(x = .data[[x_var]])) 
        dens_y <- df %>%  
            ggplot(aes(y = .data[[y_var]])) 
    } else {
        ## set up the plot aes with colors for groups
        p1 <- df %>%
            ggplot(aes(x=.data[[x_var]], y=.data[[y_var]], 
                       color = .data[[color_var]])) 
        ## make the density plots
        dens_x <- df %>%  
            ggplot(aes(x = .data[[x_var]], 
                       color = .data[[color_var]], fill=.data[[color_var]])) 
        dens_y <- df %>%  
            ggplot(aes(y = .data[[y_var]], 
                       color = .data[[color_var]], fill=.data[[color_var]])) 
    }
    
    ### add the geoms
    p1 <- p1 +
        geom_point(...) +
        theme_classic() 
    dens_x <- dens_x + 
        geom_density(show.legend = FALSE, alpha=0.3) +
        theme_void()
    dens_y <- dens_y + 
        geom_density(show.legend = FALSE, alpha=0.3) +
        theme_void() 
    
    ## customize if needed
    ## I pull the two coord_cartesian x and y together to avoid a warning 
    if(!is.null(my_xlim) & !is.null(my_ylim)) { 
        p1 <- p1 + coord_cartesian(xlim=my_xlim, ylim=my_ylim) 
    }
    
    if(!is.null(my_xlim) ) { 
        if(is.null(my_ylim)) {
            p1 <- p1 + coord_cartesian(xlim=my_xlim) 
        }
        dens_x <- dens_x + coord_cartesian(xlim=my_xlim) 
    }
    if(!is.null(my_ylim)) { 
        if(is.null(my_xlim)) {
            p1 <- p1 + coord_cartesian(ylim=my_ylim) 
        }
        dens_y <- dens_y + coord_cartesian(ylim=my_ylim) 
    }
    if(!is.null(my_color_scheme)) {
        p1 <- p1 + scale_color_manual(values=my_color_scheme)
        dens_x <- dens_x + scale_color_manual(values=my_color_scheme)
        dens_y <- dens_y + scale_color_manual(values=my_color_scheme)
    }
    ## title and subtitle go above the dens_x plot
    if(!is.null(my_title)) {
        dens_x <- dens_x + labs(title=my_title)
    }
    if(!is.null(my_subtitle)) {
        dens_x <- dens_x + labs(subtitle=my_subtitle)
    }
    
    ### maybe we don't combine
    if(!combine_plots) {
        return(list(main_plot=p1, 
                    dens_x=dens_x, 
                    dens_y=dens_y))
    }
    
    ### combine using patchwork
    all_plots <- (dens_x + plot_spacer() + p1 + dens_y) + 
        plot_layout(ncol = 2, nrow = 2, 
                    widths = c(6, 1), heights = c(1, 6), guides="collect") 
    return(all_plots)
}



#### blank_plot and blank_ggplot - can be useful as placeholders when using patchwork to combine plots

blank_plot <- function() {
    plot(1:10,1:10,"n",bty="n",xaxt="n",yaxt="n", xlab="", ylab="")
}

blank_ggplot <- function() {
    tibble(x=1:10, y=1:10) %>% 
        ggplot() +
        theme_void()
}



### lists the biggest num (default 10) items in memory
# all: can list hidden objects, too
# see also memory.profile()
lsmem <- function (num=10, all=FALSE, units="Kb") {
    z <- sapply(ls(.GlobalEnv), function(x) {
        z1 <- object.size(get(x, envir = .GlobalEnv))
        z1
    })
    z <- as.matrix(rev(sort(z))[1:num])
    didIformat <- "no"
    if (units=="Gb") {z <- round(z/ 1024^3,1) ; didIformat <- "yes"}
    if (units=="Mb") {z <- round(z/ 1024^2,1) ; didIformat <- "yes" }
    if (units=="Kb") {z <- round(z/ 1024,1)   ; didIformat <- "yes" }
    if (units=="bytes") {didIformat <- "yes" }
    if (didIformat=="no") {
        cat("\n\nWARNING - did not recognize the units you asked for. Options are Gb, Mb, Kb\n\n")
        units <- "bytes"
    }
    colnames(z)[1] <- units
    z
}