library(tidyverse)
library(patchwork)

myPlotFunc <- function(dat) {
    ggplot(dat,
           aes(x = Petal.Width, y = Petal.Length)) +
        geom_point() +
        ggtitle(unique(dat$Species))
}

## all data
iris %>% myPlotFunc()

## group_split then map to plot each (purrr::map)
temp <- iris %>%
    group_split(Species) %>%
    map(myPlotFunc)
temp[[1]]


## or use an anonymous function with the ~
temp <- iris %>%
    group_split(Species) %>%
    map(~ggplot(.,
               aes(x = Petal.Width, y = Petal.Length)) +
            geom_point() +
            ggtitle(unique(.$Species)))
temp[[1]] + temp[[2]] + temp[[3]]



