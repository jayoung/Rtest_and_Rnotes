miscellaneous_testCode
================
Janet Young

2024-10-24

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

## barplots of counts in each category

``` r
m %>% 
    ggplot(aes(x=cyl)) +
    geom_bar()
```

![](miscellaneous_testCode_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

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
    ##  [1]  1.3675113  1.7277559  1.2197515  1.2385978  1.7158086  0.6170572
    ##  [7]  0.9591685  0.1974981  1.7599008 -2.5328191
    ## 
    ## [[2]]
    ##  [1] 1.6723304 3.4909170 2.4674111 2.4541835 1.0656303 0.9991537 3.3682911
    ##  [8] 0.7938045 1.0490553 1.2877373
    ## 
    ## [[3]]
    ##  [1] 3.855963 2.766923 4.198429 4.844597 3.248347 3.130018 3.014173 3.124036
    ##  [9] 1.044734 2.555878
    ## 
    ## [[4]]
    ##  [1] 2.570264 3.104092 4.960435 4.795026 2.634247 3.431036 2.744862 4.519494
    ##  [9] 4.679337 3.363718
    ## 
    ## [[5]]
    ##  [1] 4.293316 5.284948 4.428178 4.883492 4.543018 3.989738 5.256036 4.641714
    ##  [9] 3.326366 3.894331
    ## 
    ## [[6]]
    ##  [1] 4.074849 4.557792 5.870403 4.814594 8.762370 5.244945 6.106749 5.136882
    ##  [9] 4.921977 6.352757
    ## 
    ## [[7]]
    ##  [1] 7.591354 6.723504 5.875923 6.782765 9.621716 6.756472 6.780931 7.403761
    ##  [9] 6.644453 7.355384
    ## 
    ## [[8]]
    ##  [1] 6.734067 7.775191 8.191314 7.501943 8.275420 6.549329 6.987847 9.774830
    ##  [9] 7.576214 8.282886
    ## 
    ## [[9]]
    ##  [1]  6.776016  8.652995  7.862085  8.563511  9.680577 10.817991  8.074243
    ##  [8]  8.732927  9.233615  8.966438
    ## 
    ## [[10]]
    ##  [1] 11.025289 11.539361 10.809502 10.904747 11.521591  9.431730  9.061578
    ##  [8]  8.814568 10.811812 11.039821

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

# some tidbit for Phoebe’s gtf/gff question

``` r
library(rtracklayer)
```

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     intersect, setdiff, union

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, table,
    ##     tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     second, second<-

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:lubridate':
    ## 
    ##     %within%

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## Loading required package: GenomeInfoDb

``` r
?GFFFile
```

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
# devtools::install_github("sfirke/packagemetrics")
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
<span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: #56A33E; width: 100.00%">1759397</span>
</td>
<td style="text-align:right;">
<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">4.8k</span>
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
<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">11.4</span>
</td>
<td style="text-align:right;">
<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">0.1</span>
</td>
<td style="text-align:right;">
<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #1cc2e3">271</span>
</td>
<td style="text-align:right;">
<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #1cc2e3">1</span>
</td>
<td style="text-align:right;">
<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #1cc2e3">4294</span>
</td>
</tr>
<tr>
<td style="text-align:right;">
<span style="font-weight: bold">data.table</span>
</td>
<td style="text-align:right;">
2024-10-10
</td>
<td style="text-align:right;">
<span style="display: inline-block; direction: rtl; unicode-bidi: plaintext; border-radius: 4px; padding-right: 2px; background-color: #56A33E; width: 52.23%">918913</span>
</td>
<td style="text-align:right;">
<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">3.6k</span>
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
<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">0.1</span>
</td>
<td style="text-align:right;">
<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">151</span>
</td>
<td style="text-align:right;">
<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #1cc2e3">1</span>
</td>
<td style="text-align:right;">
<span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">1684</span>
</td>
</tr>
</tbody>
</table>
