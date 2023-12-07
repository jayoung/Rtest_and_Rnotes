library(tidyverse)

m <- mtcars %>% as_tibble()

# show counts in each category
m %>% 
    count(cyl)

# barplots of counts in each category
m %>% 
    ggplot(aes(x=cyl)) +
    geom_bar()


###########  for Phoebe's gtf/gff question

library(rtracklayer)
?GFFFile