## https://github.com/YuLab-SMU/ggbreak/issues/82

library(ggplot2)
library(patchwork)
library(ggbreak)

set.seed(2019-01-19)

d <- data.frame(
    x = 1:20,
    y = c(rnorm(5) + 4, rnorm(5) + 20, rnorm(5) + 5, rnorm(5) + 22)
)

p <- ggplot(d, aes(x, y)) + 
    geom_col() +
    labs(title="p - no axis break")
x <- p +
    scale_y_break(c(7, 17)) +
    labs(title="x - has axis break")


## axis breaks are lost with x + p
x + p


## axis breaks are retained with p + x
p + x



sessionInfo()

# library(reprex)