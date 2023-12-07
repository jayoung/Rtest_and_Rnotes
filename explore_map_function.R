library(tidyverse)
library(patchwork)

temp <- iris %>% 
    group_by(Species) %>% 
    map(plots=ggplot(data=.) +
          aes(x=Petal.Width, y=Petal.Length) + 
           geom_point() + 
           ggtitle(unique(.$Species)))


temp$plots[[1]]


myPlotFunc <- function(x) {
    ggplot(data=x) +
        aes(x=Petal.Width, y=Petal.Length) + 
        geom_point() + 
        ggtitle(unique(x$Species))
}

temp <- iris %>% 
    group_split(Species) %>% 
    map(function(x) {
        ggplot(data=x) +
            aes(x=Petal.Width, y=Petal.Length) + 
            geom_point() + 
            ggtitle(unique(x$Species))
    })
