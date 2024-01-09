########### basic tidyverse

library(tidyverse)

m <- mtcars %>% as_tibble()

## show counts in each category
m %>% 
    count(cyl)

## barplots of counts in each category
m %>% 
    ggplot(aes(x=cyl)) +
    geom_bar()


## tribble - a way to manually create small tibbles
tribble(
    ~colA, ~colB,
    "a",   1,
    "b",   2,
    "c",   3
)
# same as:
tibble(colA=c("A","B","C"), 
       colB=1:3)


## the purrr::map functions are a bit like lapply / apply
1:10 %>%
    map(rnorm, n = 10)

# do something to each column
mtcars %>% map_dbl(sum)

# split into list, then map.  
# map_dbl simplify output, in this case to a dbl (or numeric)
# same for map_lgl(), map_int() and map_chr()
mtcars %>%
    split(.$cyl) %>%
    map(~ lm(mpg ~ wt, data = .x)) %>%
    map(summary) %>%
    map_dbl("r.squared")

# map_dfr tries to return a data.frame (binding by rows)  (map_dfc is similar but binds by columns)
mtcars %>%
    split(.$cyl) %>%
    map(~ lm(mpg ~ wt, data = .x)) %>%
    map_dfr(~ as.data.frame(t(as.matrix(coef(.)))))


########### for Phoebe's gtf/gff question

library(rtracklayer)
?GFFFile


########### "bang-bang" or !!
# this performs "name injection"
# explained more here https://dplyr.tidyverse.org/articles/programming.html#name-injection

## in dplyr::rename example:
#  !! interprets the variable and then we have := (instead of =)

newVarName <- "sepalLen_new"
iris %>% 
    dplyr::rename(!!newVarName := Sepal.Length) %>% 
    head()


########### The 'embracing' operator (`{{ }}`) 
## useful when we want to pass in variable(s) that we want to be interpreted before it's used
## related to !! and !!!


library(dplyr)

### this does not work
get_var0 <- function(data, column, value) {
    data %>% filter(column == value)
}
get_var0(mtcars, cyl, 6)
#> Error: Problem with `filter()` input `..1`.
#> x object 'cyl' not found
#> i Input `..1` is `column == value`.
get_var0(mtcars, "cyl", 6)
# returns empty tbl, should be 7

### this works
get_var1 <- function(data, column, value) {
    data %>% filter({{ column }} == value)
}
get_var1(mtcars, cyl, 6)


# more discussion here: 
# https://rlang.r-lib.org/reference/embrace-operator.html
# https://rlang.r-lib.org/reference/topic-data-mask.html
# https://adv-r.hadley.nz/quasiquotation.html

## similar, but not identical to !! ("bang-bang")
# https://www.r-bloggers.com/2019/07/bang-bang-how-to-program-with-dplyr/

## similar, but not identical to !!! ("bang-bang-bang")
# https://www.reddit.com/r/Rlanguage/comments/g5m5bh/what_does_the_bang_bang_bang_do/?rdt=45180

?rlang::`!!`
?rlang::`!!!`

### !!! example: (from https://www.reddit.com/r/Rlanguage/comments/g5m5bh/what_does_the_bang_bang_bang_do/?rdt=45180)

# Let's say we want to select the first 3 columns of a data frame. So you can do something like this:
test_df <- tibble(a = 1, b = 1, c = 1, d = 1)
test_df %>%
  select(1, 2, 3)

# Easy enough. Now let's say we have a list of values that we want to use to replicate the code above. Now if you pass this to the select function it fails:
our_list <- list(1, 2, 3)
test_df %>%
    select(our_list)
# That's because that code essentially translates to this, which doesn't work.
test_df %>%
  select(list(1, 2, 3))

# What we need to do is "unpack" the list using !!!.
test_df %>%
    select(!!!our_list)
# which translates to this:
test_df %>%
    select(1, 2, 3)
