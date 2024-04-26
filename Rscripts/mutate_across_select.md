mutate_across
================
Janet Young

2024-04-26

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
    select(-ends_with("mm")) %>% # remove old columns here
    head(3)
```

    ## # A tibble: 3 × 8
    ##   species island    body_mass_g sex     year bill_length_cm bill_depth_cm
    ##   <fct>   <fct>           <int> <fct>  <int>          <dbl>         <dbl>
    ## 1 Adelie  Torgersen        3750 male    2007           3.91          1.87
    ## 2 Adelie  Torgersen        3800 female  2007           3.95          1.74
    ## 3 Adelie  Torgersen        3250 female  2007           4.03          1.8 
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
    select(union (ends_with("_mm"), ends_with("_sq") )) %>% 
    head(3)
```

    ## # A tibble: 3 × 6
    ##   bill_length_mm bill_depth_mm flipper_length_mm bill_length_mm_sq
    ##            <dbl>         <dbl>             <int>             <dbl>
    ## 1           39.1          18.7               181             1529.
    ## 2           39.5          17.4               186             1560.
    ## 3           40.3          18                 195             1624.
    ## # ℹ 2 more variables: bill_depth_mm_sq <dbl>, flipper_length_mm_sq <dbl>

# A bit more on using `select()`

How to select columns based on more than one condition - use
`union`/`intersect`, or `|`/`&`

``` r
penguins |> 
    select(union (ends_with("_mm"), ends_with("_g") ))  %>% 
    head(3)
```

    ## # A tibble: 3 × 4
    ##   bill_length_mm bill_depth_mm flipper_length_mm body_mass_g
    ##            <dbl>         <dbl>             <int>       <int>
    ## 1           39.1          18.7               181        3750
    ## 2           39.5          17.4               186        3800
    ## 3           40.3          18                 195        3250

``` r
## same as this:
# penguins |> 
# select( ends_with("_mm") | ends_with("_sq")) 
```

``` r
penguins |>
    select(intersect (ends_with("_mm"), matches("_length") )) %>% 
    head(3)
```

    ## # A tibble: 3 × 2
    ##   bill_length_mm flipper_length_mm
    ##            <dbl>             <int>
    ## 1           39.1               181
    ## 2           39.5               186
    ## 3           40.3               195

``` r
## same as this:
# penguins |> 
#     select( ends_with("_mm") & matches("_length") )  %>% 
    # head(3)
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
    ) %>% 
    head(3)
```

    ## # A tibble: 3 × 5
    ##   bill_length_mm bill_depth_mm flipper_length_mm body_mass_g  year
    ##            <dbl>         <dbl>             <int>       <int> <int>
    ## 1           39.1          18.7               181        3750  2007
    ## 2           39.5          17.4               186        3800  2007
    ## 3           40.3          18                 195        3250  2007
