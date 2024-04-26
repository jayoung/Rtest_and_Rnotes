mutate_across
================
Janet Young

2024-04-26

``` r
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.0     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

Get example data

``` r
penguins <- palmerpenguins::penguins
```

Use mutate across (example from the help)

``` r
## make a function and use it in mutate across:
divide_by_ten <- function(x) x / 10

## do mutate(across())
penguins |>
    mutate(
        across(
            .cols = ends_with("mm"),   ## everything() chooses all columns
            .fns = divide_by_ten,
            ## use .names arg to figure out what output colnames will be. Empty names arg means values in current columns get replaced
            .names = '{str_remove(.col, "_mm")}_cm'
        )
    ) |>
    select(-ends_with("mm")) # remove old columns here
```

    ## # A tibble: 344 × 8
    ##    species island    body_mass_g sex     year bill_length_cm bill_depth_cm
    ##    <fct>   <fct>           <int> <fct>  <int>          <dbl>         <dbl>
    ##  1 Adelie  Torgersen        3750 male    2007           3.91          1.87
    ##  2 Adelie  Torgersen        3800 female  2007           3.95          1.74
    ##  3 Adelie  Torgersen        3250 female  2007           4.03          1.8 
    ##  4 Adelie  Torgersen          NA <NA>    2007          NA            NA   
    ##  5 Adelie  Torgersen        3450 female  2007           3.67          1.93
    ##  6 Adelie  Torgersen        3650 male    2007           3.93          2.06
    ##  7 Adelie  Torgersen        3625 female  2007           3.89          1.78
    ##  8 Adelie  Torgersen        4675 male    2007           3.92          1.96
    ##  9 Adelie  Torgersen        3475 <NA>    2007           3.41          1.81
    ## 10 Adelie  Torgersen        4250 <NA>    2007           4.2           2.02
    ## # ℹ 334 more rows
    ## # ℹ 1 more variable: flipper_length_cm <dbl>

More practise:

``` r
myFunc <- function(x) {x^2}
penguins %>% 
    mutate(across( ends_with("_mm"),
                   myFunc,
                   # .names = '{str_remove(.col, "_mm")}_cm'
                   .names = '{.col}_sq'
    )) %>% 
    select(union (ends_with("_mm"), ends_with("_sq") )) 
```

    ## # A tibble: 344 × 6
    ##    bill_length_mm bill_depth_mm flipper_length_mm bill_length_mm_sq
    ##             <dbl>         <dbl>             <int>             <dbl>
    ##  1           39.1          18.7               181             1529.
    ##  2           39.5          17.4               186             1560.
    ##  3           40.3          18                 195             1624.
    ##  4           NA            NA                  NA               NA 
    ##  5           36.7          19.3               193             1347.
    ##  6           39.3          20.6               190             1544.
    ##  7           38.9          17.8               181             1513.
    ##  8           39.2          19.6               195             1537.
    ##  9           34.1          18.1               193             1163.
    ## 10           42            20.2               190             1764 
    ## # ℹ 334 more rows
    ## # ℹ 2 more variables: bill_depth_mm_sq <dbl>, flipper_length_mm_sq <dbl>

# A bit more on using `select()`

How to select columns based on more than one condition - use
`union`/`intersect`, or `|`/`&`

``` r
penguins |> 
    select(union (ends_with("_mm"), ends_with("_g") )) 
```

    ## # A tibble: 344 × 4
    ##    bill_length_mm bill_depth_mm flipper_length_mm body_mass_g
    ##             <dbl>         <dbl>             <int>       <int>
    ##  1           39.1          18.7               181        3750
    ##  2           39.5          17.4               186        3800
    ##  3           40.3          18                 195        3250
    ##  4           NA            NA                  NA          NA
    ##  5           36.7          19.3               193        3450
    ##  6           39.3          20.6               190        3650
    ##  7           38.9          17.8               181        3625
    ##  8           39.2          19.6               195        4675
    ##  9           34.1          18.1               193        3475
    ## 10           42            20.2               190        4250
    ## # ℹ 334 more rows

``` r
## same as this:
# penguins |> 
# select( ends_with("_mm") | ends_with("_sq")) 
```

``` r
penguins |>
    select(intersect (ends_with("_mm"), matches("_length") ))
```

    ## # A tibble: 344 × 2
    ##    bill_length_mm flipper_length_mm
    ##             <dbl>             <int>
    ##  1           39.1               181
    ##  2           39.5               186
    ##  3           40.3               195
    ##  4           NA                  NA
    ##  5           36.7               193
    ##  6           39.3               190
    ##  7           38.9               181
    ##  8           39.2               195
    ##  9           34.1               193
    ## 10           42                 190
    ## # ℹ 334 more rows

``` r
## same as this:
# penguins |> 
#     select( ends_with("_mm") & matches("_length") ) 
```

How to select columns based on class of data within them, using
`where()`:

Notice that you don’t actually call is.numeric in your code (i.e. there
are no paretheses after is.numeric). That’s because you partner up the
function is.numeric with the where function. This way, the where
function just knows the name of its partner (in this case is.numeric )
and can let its partner do the checking whenever where feels like it.

``` r
penguins |>
    select(
        where(is.numeric)
    )
```

    ## # A tibble: 344 × 5
    ##    bill_length_mm bill_depth_mm flipper_length_mm body_mass_g  year
    ##             <dbl>         <dbl>             <int>       <int> <int>
    ##  1           39.1          18.7               181        3750  2007
    ##  2           39.5          17.4               186        3800  2007
    ##  3           40.3          18                 195        3250  2007
    ##  4           NA            NA                  NA          NA  2007
    ##  5           36.7          19.3               193        3450  2007
    ##  6           39.3          20.6               190        3650  2007
    ##  7           38.9          17.8               181        3625  2007
    ##  8           39.2          19.6               195        4675  2007
    ##  9           34.1          18.1               193        3475  2007
    ## 10           42            20.2               190        4250  2007
    ## # ℹ 334 more rows
