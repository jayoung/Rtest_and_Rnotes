miscellaneous_testCode
================
Janet Young

2025-08-22

# basic tidyverse

``` r
m <- mtcars %>% as_tibble()

## show counts in each category
m %>% 
    count(cyl)
```

    ## # A tibble: 3 × 2
    ##     cyl     n
    ##   <dbl> <int>
    ## 1     4    11
    ## 2     6     7
    ## 3     8    14

## separate_longer trick (for list-like columns), and str_replace_all trick (multiple find-replaces)

From [here](https://rfortherestofus.com/2024/04/seperate-fcts)

``` r
details <- readr::read_csv("https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2022/2022-01-25/details.csv",
                           show_col_types = FALSE)
board_games <- details |>
    select(id, name = primary, boardgamecategory)

board_games %>% 
    head(3)
```

    ## # A tibble: 3 × 3
    ##      id name        boardgamecategory                                  
    ##   <dbl> <chr>       <chr>                                              
    ## 1 30549 Pandemic    ['Medical']                                        
    ## 2   822 Carcassonne ['City Building', 'Medieval', 'Territory Building']
    ## 3    13 Catan       ['Economic', 'Negotiation']

``` r
board_games |>
    separate_longer_delim(cols = boardgamecategory, delim = ", ") |>
    mutate(
        boardgamecategory = str_replace_all(
            boardgamecategory,
            c(
                # pattern_to_replace = replacement
                # We have to wrap the names in backticks because all of the things we want to replace are special characters and R doesn’t like them in vector names if they’re not in backticks
                `[` = "",
                `]` = "",
                `"` = "",
                `'` = ""
            ) |> 
                coll() ## coll() tells str_replace_all() that we explicitly do not want to use regular expressions.
        ))
```

    ## # A tibble: 56,915 × 3
    ##       id name        boardgamecategory 
    ##    <dbl> <chr>       <chr>             
    ##  1 30549 Pandemic    Medical           
    ##  2   822 Carcassonne City Building     
    ##  3   822 Carcassonne Medieval          
    ##  4   822 Carcassonne Territory Building
    ##  5    13 Catan       Economic          
    ##  6    13 Catan       Negotiation       
    ##  7 68448 7 Wonders   Ancient           
    ##  8 68448 7 Wonders   Card Game         
    ##  9 68448 7 Wonders   City Building     
    ## 10 68448 7 Wonders   Civilization      
    ## # ℹ 56,905 more rows

## barplots of counts in each category

``` r
m %>% 
    ggplot(aes(x=cyl)) +
    geom_bar()
```

![](miscellaneous_testCode_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## tribble - a way to manually create small tibbles, row-wise

``` r
tribble(
    ~colA, ~colB,
    "a",   1,
    "b",   2,
    "c",   3
)
```

    ## # A tibble: 3 × 2
    ##   colA   colB
    ##   <chr> <dbl>
    ## 1 a         1
    ## 2 b         2
    ## 3 c         3

``` r
# same as:
tibble(colA=c("A","B","C"), 
       colB=1:3)
```

    ## # A tibble: 3 × 2
    ##   colA   colB
    ##   <chr> <int>
    ## 1 A         1
    ## 2 B         2
    ## 3 C         3

## the purrr::map functions are a bit like lapply / apply

``` r
1:10 %>%
    map(rnorm, n = 10)
```

    ## [[1]]
    ##  [1]  1.3603895 -1.4114722 -1.0326192  1.2545541  1.3357514  0.0215660
    ##  [7] -0.4046363 -0.8406554  0.3354220 -1.3492406
    ## 
    ## [[2]]
    ##  [1] 2.164291 3.102322 1.562906 2.301881 3.624128 1.099051 4.085332 2.602163
    ##  [9] 1.281087 1.131071
    ## 
    ## [[3]]
    ##  [1] 2.374676 2.851518 3.192542 3.756049 2.779175 4.756329 3.603554 2.468738
    ##  [9] 3.760330 4.457496
    ## 
    ## [[4]]
    ##  [1] 3.967028 5.176950 4.572501 3.101835 4.800949 3.863271 4.918282 3.395852
    ##  [9] 4.658805 4.439076
    ## 
    ## [[5]]
    ##  [1] 4.424805 4.316911 4.696888 5.511118 5.626091 5.101665 5.191551 5.336532
    ##  [9] 5.883760 5.616675
    ## 
    ## [[6]]
    ##  [1] 6.503958 6.170519 5.857271 5.479236 5.982719 4.050278 4.480401 6.318529
    ##  [9] 6.058498 6.383746
    ## 
    ## [[7]]
    ##  [1] 7.871752 5.106919 9.122487 7.108571 9.814153 7.334294 6.854749 6.888949
    ##  [9] 7.705922 7.304566
    ## 
    ## [[8]]
    ##  [1] 7.315753 9.581842 8.256553 8.612141 7.281892 8.788097 7.879163 7.977126
    ##  [9] 8.482422 9.303036
    ## 
    ## [[9]]
    ##  [1]  8.671588  8.636743  8.768620 10.070894  8.759060  7.809230  8.248732
    ##  [8]  7.975862  8.868676  8.062077
    ## 
    ## [[10]]
    ##  [1]  9.621413  9.115014 10.345091 10.842908 11.192355 10.015732  9.781999
    ##  [8] 11.130238 11.751150  9.672034

``` r
# do something to each column
mtcars %>% map_dbl(sum)
```

    ##      mpg      cyl     disp       hp     drat       wt     qsec       vs 
    ##  642.900  198.000 7383.100 4694.000  115.090  102.952  571.160   14.000 
    ##       am     gear     carb 
    ##   13.000  118.000   90.000

``` r
# split into list, then map.  
# map_dbl simplify output, in this case to a dbl (or numeric)
# same for map_lgl(), map_int() and map_chr()
mtcars %>%
    split(.$cyl) %>%
    map(~ lm(mpg ~ wt, data = .x)) %>%
    map(summary) %>%
    map_dbl("r.squared")
```

    ##         4         6         8 
    ## 0.5086326 0.4645102 0.4229655

``` r
# map_dfr tries to return a data.frame (binding by rows)  (map_dfc is similar but binds by columns)
mtcars %>%
    split(.$cyl) %>%
    map(~ lm(mpg ~ wt, data = .x)) %>%
    map_dfr(~ as.data.frame(t(as.matrix(coef(.)))))
```

    ##   (Intercept)        wt
    ## 1    39.57120 -5.647025
    ## 2    28.40884 -2.780106
    ## 3    23.86803 -2.192438

`setNames()` is a way to set names on a list - can be used in a pipe as
an alternative to adding a separate `names(x) <- y` after a set of piped
commands

# variable names - using variables

## bang-bang

“bang-bang” , a.ka. !! performs “name injection”

it’s explained more
[here](https://dplyr.tidyverse.org/articles/programming.html#name-injection)

To see more info: ?rlang::`!!`

here I use it in dplyr::rename - we use `!!` to interpret the variable
and then we have to use `:=` (instead of =)

``` r
newVarName <- "sepalLen_new"
iris %>% 
    dplyr::rename(!!newVarName := Sepal.Length) %>% 
    head()
```

    ##   sepalLen_new Sepal.Width Petal.Length Petal.Width Species
    ## 1          5.1         3.5          1.4         0.2  setosa
    ## 2          4.9         3.0          1.4         0.2  setosa
    ## 3          4.7         3.2          1.3         0.2  setosa
    ## 4          4.6         3.1          1.5         0.2  setosa
    ## 5          5.0         3.6          1.4         0.2  setosa
    ## 6          5.4         3.9          1.7         0.4  setosa

There’s also something called
[“bang-bang-bang”](https://www.reddit.com/r/Rlanguage/comments/g5m5bh/what_does_the_bang_bang_bang_do/?rdt=45180)

## embracing operator

The ‘embracing’ operator (`{{ }}`) is useful when we want to pass in
variable(s) that we want to be interpreted before being used. It’s
related to !! and !!!

more discussion
[here](https://rlang.r-lib.org/reference/embrace-operator.html),
[here](https://rlang.r-lib.org/reference/topic-data-mask.html) and
[here](https://adv-r.hadley.nz/quasiquotation.html)

The embracing operator is [similar
to](https://www.r-bloggers.com/2019/07/bang-bang-how-to-program-with-dplyr/),
but not identical to !! (“bang-bang”). It’s also [similar
to](https://www.reddit.com/r/Rlanguage/comments/g5m5bh/what_does_the_bang_bang_bang_do/?rdt=45180),
but not identical to !!! (“bang-bang-bang”)

Example:

``` r
### this wouldn't work
## first define a function
get_var0 <- function(data, column, value) {
    data %>% filter(column == value)
}
## then use it - gives an error this way
# get_var0(mtcars, cyl, 6)
#> Error: Problem with `filter()` input `..1`.
#> x object 'cyl' not found
#> i Input `..1` is `column == value`.

## gives empty tbl output this way (wrong, should be 7):
# get_var0(mtcars, "cyl", 6)

### this DOES work
get_var1 <- function(data, column, value) {
    data %>% filter({{ column }} == value)
}
get_var1(mtcars, cyl, 6)
```

    ##                 mpg cyl  disp  hp drat    wt  qsec vs am gear carb
    ## Mazda RX4      21.0   6 160.0 110 3.90 2.620 16.46  0  1    4    4
    ## Mazda RX4 Wag  21.0   6 160.0 110 3.90 2.875 17.02  0  1    4    4
    ## Hornet 4 Drive 21.4   6 258.0 110 3.08 3.215 19.44  1  0    3    1
    ## Valiant        18.1   6 225.0 105 2.76 3.460 20.22  1  0    3    1
    ## Merc 280       19.2   6 167.6 123 3.92 3.440 18.30  1  0    4    4
    ## Merc 280C      17.8   6 167.6 123 3.92 3.440 18.90  1  0    4    4
    ## Ferrari Dino   19.7   6 145.0 175 3.62 2.770 15.50  0  1    5    6

### !!! example:

from
[here](https://www.reddit.com/r/Rlanguage/comments/g5m5bh/what_does_the_bang_bang_bang_do/?rdt=45180)

To see more info: ?rlang::`!!!`

``` r
## Let's say we want to select the first 3 columns of a data frame. So you can do something like this:
test_df <- tibble(a = 1, b = 1, c = 1, d = 1)
test_df %>%
    select(1, 2, 3)
```

    ## # A tibble: 1 × 3
    ##       a     b     c
    ##   <dbl> <dbl> <dbl>
    ## 1     1     1     1

``` r
## Easy enough. Now let's say we have a list of values that we want to use to replicate the code above. Now if you pass this to the select function it fails:
our_list <- list(1, 2, 3)
# test_df %>%
#     select(our_list)
## That's because that code essentially translates to this, which doesn't work.
# test_df %>%
#     select(list(1, 2, 3))

# What we need to do is "unpack" the list using !!!.
test_df %>%
    select(!!!our_list)
```

    ## # A tibble: 1 × 3
    ##       a     b     c
    ##   <dbl> <dbl> <dbl>
    ## 1     1     1     1

``` r
# which translates to this:
# test_df %>%
#     select(1, 2, 3)
```

# adding metadata / attributes

I had a use case where I wanted to associate a few bits of info with a
tibble of data. The actual example is that my tibble contained
genome-wide t-test results, and I wanted to record eaxctly what method
I’d used to do the t-tests)

Note that `attr()` and `atrributes()` are different functions

``` r
# make example tibble:
a <- 1:5
b <- tibble(a, a * 2)

### add attributes under the methods tag (there can be:
attr(b, "methods") <- "here's some text about the method"
attr(b, "details") <- "here's some more detail about the method"

## access the metadata:
attr(b, "methods")
```

    ## [1] "here's some text about the method"

``` r
## show all attributes:
attributes(b)
```

    ## $class
    ## [1] "tbl_df"     "tbl"        "data.frame"
    ## 
    ## $row.names
    ## [1] 1 2 3 4 5
    ## 
    ## $names
    ## [1] "a"     "a * 2"
    ## 
    ## $methods
    ## [1] "here's some text about the method"
    ## 
    ## $details
    ## [1] "here's some more detail about the method"

<https://github.com/sfirke/packagemetrics?tab=readme-ov-file>

``` r
devtools::install_github("sfirke/packagemetrics")
```

    ## Using GitHub PAT from the git credential store.

    ## Skipping install of 'packagemetrics' from a github remote, the SHA1 (431cd4f7) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
library("packagemetrics")
dplyr_and_dt <- package_list_metrics(c("dplyr", "data.table"))
# data frame
```

``` r
metrics_table(dplyr_and_dt)
```

    ## Warning in gradient(as.numeric(x), ...): NAs introduced by coercion

<table class="table table-condensed">

<thead>

<tr>

<th style="text-align:right;">

package
</th>

<th style="text-align:right;">

published
</th>

<th style="text-align:right;">

dl_last_month
</th>

<th style="text-align:right;">

stars
</th>

<th style="text-align:right;">

tidyverse_happy
</th>

<th style="text-align:right;">

has_tests
</th>

<th style="text-align:right;">

vignette
</th>

<th style="text-align:right;">

last_commit
</th>

<th style="text-align:right;">

last_issue_closed
</th>

<th style="text-align:right;">

contributors
</th>

<th style="text-align:right;">

depends_count
</th>

<th style="text-align:right;">

reverse_count
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

<span style="font-weight: bold">dplyr </span>
</td>

<td style="text-align:right;">

2023-11-17
</td>

<td style="text-align:right;">

<span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: #56A33E; width: 100.00%">1346628</span>
</td>

<td style="text-align:right;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">4.9k</span>
</td>

<td style="text-align:right;">

<span style="color: purple"> <i class="glyphicon glyphicon-glass"></i>
</span>
</td>

<td style="text-align:right;">

<span style="color: green"> <i class="glyphicon glyphicon-ok"></i>
</span>
</td>

<td style="text-align:right;">

<span style="color: green"> <i class="glyphicon glyphicon-ok"></i>
</span>
</td>

<td style="text-align:right;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">21.5</span>
</td>

<td style="text-align:right;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #f06b13"></span>
</td>

<td style="text-align:right;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #1cc2e3">271</span>
</td>

<td style="text-align:right;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #1cc2e3">1</span>
</td>

<td style="text-align:right;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #1cc2e3">4845</span>
</td>

</tr>

<tr>

<td style="text-align:right;">

<span style="font-weight: bold">data.table</span>
</td>

<td style="text-align:right;">

2025-07-10
</td>

<td style="text-align:right;">

<span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: #56A33E; width: 52.20%">702982</span>
</td>

<td style="text-align:right;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">3.8k</span>
</td>

<td style="text-align:right;">

<span style="color: white"> <i class="glyphicon glyphicon-glass"></i>
</span>
</td>

<td style="text-align:right;">

<span style="color: red"> <i class="glyphicon glyphicon-remove"></i>
</span>
</td>

<td style="text-align:right;">

<span style="color: green"> <i class="glyphicon glyphicon-ok"></i>
</span>
</td>

<td style="text-align:right;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">
</span>
</td>

<td style="text-align:right;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #f06b13"></span>
</td>

<td style="text-align:right;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">167</span>
</td>

<td style="text-align:right;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #1cc2e3">1</span>
</td>

<td style="text-align:right;">

<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">1783</span>
</td>

</tr>

</tbody>

</table>
