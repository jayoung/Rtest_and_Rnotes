# from https://www.njtierney.com/post/2023/11/10/how-to-get-good-with-r/


## define some functions
mean_center <- function(variable){
    variable - mean(variable, na.rm = TRUE)
}
scale_sd <- function(variable){
    variable / sd(variable, na.rm = TRUE)
}
center_scale <- function(variable){
    centered <- mean_center(variable)
    centered_scaled <- scale_sd(centered)
    centered_scaled
}


## example data

dat <- tibble::tibble(
    weight = rnorm(10, mean = 80, sd = 6),
    height = rnorm(10, mean = 165, sd = 10),
)

## use functions
dat$weight_0 <- center_scale(dat$weight)

dat$height_0 <- center_scale(dat$height)


## add browser() to a function:
# need to remember to remove it later
mean_center <- function(variable){
    browser()
    variable - mean(variable, na.rm = TRUE)
}
mean_center(1:10)

## or use debug(myFunc) and undebug(myFunc)
mean_center <- function(variable){
    variable - mean(variable, na.rm = TRUE)
}
debug(mean_center)
mean_center(1:10)
undebug(mean_center)


## or debugonce() does debug(myFunc) and then turns off debugging


### example data 2, using tidyverse
library(tidyverse)
diamonds %>%
    mutate(
        price_per_carat = price / carat
    ) %>% 
    group_by(
        cut
    ) %>% 
    summarise(
        price_mean = mean(price_per_carat),
        price_sd = sd(price_per_carat),
        mean_color = mean(color)
    )

## need to have >1 step for debugging to be effective
myFunc <- function(myTbl) {
    output <- myTbl
    output <- output %>%
        mutate(
            price_per_carat = price / carat
        ) 
    output <- output %>% 
        group_by(
            cut
        )
    output <- output %>% 
        summarise(
            price_mean = mean(price_per_carat),
            price_sd = sd(price_per_carat),
            mean_color = mean(color)
        )
    return(output)
}
debug(myFunc)
myFunc(diamonds)
undebug(myFunc)




## need to have >1 step for debugging to be effective
myFunc2 <- function(myTbl) {
    browser()
    output <- myTbl
    output <- output %>%
        mutate(
            price_per_carat = price / carat
        ) 
    output <- output %>% 
        group_by(
            cut
        )
    output <- output %>% 
        summarise(
            price_mean = mean(price_per_carat),
            price_sd = sd(price_per_carat),
            mean_color = mean(color)
        )
    return(output)
}
myFunc2(diamonds)

