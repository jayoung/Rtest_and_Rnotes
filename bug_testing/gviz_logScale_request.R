# April 15, 2024
# https://github.com/ivanek/Gviz/issues/91

library(tidyverse)
library(patchwork)

dat <- data.frame(x=c(1,2,3,4),
                  y=c(1,10,100,1000) )

p1 <- dat %>% 
    ggplot(aes(x=x, y=y)) +
    geom_col() +
    scale_y_log10() +
    labs(title="I like this y-axis")

p2 <- dat %>% 
    ggplot(aes(x=x, y=log10(y))) +
    geom_col()  +
    labs(title="better than this one")

p1 + p2


ggsave(filename="bug_testing/gviz_logScale_request.png", height=5,width=5)
