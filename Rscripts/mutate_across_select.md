mutate_across
================
Janet Young

2024-04-26

Get example data

``` r
penguins <- palmerpenguins::penguins
```

Use mutate across with a simple built-in function (example from the
help)

``` r
penguins |>
    mutate(
        across(ends_with("mm"), log10,
               .names = '{.col}_log10'
        )
    ) |>
    select(matches("depth") & (ends_with("mm") | ends_with("log10"))) |> # remove old columns here
    head(3)
```

    ## # A tibble: 3 × 2
    ##   bill_depth_mm bill_depth_mm_log10
    ##           <dbl>               <dbl>
    ## 1          18.7                1.27
    ## 2          17.4                1.24
    ## 3          18                  1.26

Use mutate across with a user-defined function (example from the help)

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
    select(matches("depth") & (ends_with("mm") | ends_with("cm"))) |> # remove old columns here
    head(3)
```

    ## # A tibble: 3 × 2
    ##   bill_depth_mm bill_depth_cm
    ##           <dbl>         <dbl>
    ## 1          18.7          1.87
    ## 2          17.4          1.74
    ## 3          18            1.8

Same thing but without defining the function

``` r
penguins |>
    mutate(
        across(
            .cols = ends_with("mm"),   ## everything() chooses all columns
            .fns = function(x) x / 10,
            ## use .names arg to figure out what output colnames will be. Empty names arg means values in current columns get replaced
            .names = '{str_remove(.col, "_mm")}_cm'
        )
    ) |>
    select(matches("depth") & (ends_with("mm") | ends_with("cm"))) |> # remove old columns here
    head(3)
```

    ## # A tibble: 3 × 2
    ##   bill_depth_mm bill_depth_cm
    ##           <dbl>         <dbl>
    ## 1          18.7          1.87
    ## 2          17.4          1.74
    ## 3          18            1.8

More practise:

``` r
myFunc <- function(x) {x^2}
penguins |> 
    mutate(across( ends_with("_mm"),
                   myFunc,
                   .names = '{.col}_sq'
    )) |> 
    select(union (ends_with("_mm"), ends_with("_sq") )) |> 
    head(3)
```

    ## # A tibble: 3 × 6
    ##   bill_length_mm bill_depth_mm flipper_length_mm bill_length_mm_sq
    ##            <dbl>         <dbl>             <int>             <dbl>
    ## 1           39.1          18.7               181             1529.
    ## 2           39.5          17.4               186             1560.
    ## 3           40.3          18                 195             1624.
    ## # ℹ 2 more variables: bill_depth_mm_sq <dbl>, flipper_length_mm_sq <dbl>

We can also use the ~ style way to specify a function on the fly, in
which case we use `.x` as a standin for the data:

``` r
penguins |> 
    mutate(across( ends_with("_mm"),
                   ~ (.x ^ 2),
                   .names = '{.col}_sq'
    )) |> 
    select(intersect(matches("length"), union(ends_with("_mm"), ends_with("_sq") ))) |> 
    head(3)
```

    ## # A tibble: 3 × 4
    ##   bill_length_mm flipper_length_mm bill_length_mm_sq flipper_length_mm_sq
    ##            <dbl>             <int>             <dbl>                <dbl>
    ## 1           39.1               181             1529.                32761
    ## 2           39.5               186             1560.                34596
    ## 3           40.3               195             1624.                38025

# A bit more on using `select()`

How to select columns based on more than one condition - use
`union`/`intersect`, or `|`/`&`

``` r
## two alternative ways to write the same thing:
# penguins |>
#     select( ends_with("_mm") | ends_with("_sq")) |> 
#     head(3)
penguins |> 
    select(union (ends_with("_mm"), ends_with("_g") )) |> 
    head(3)
```

    ## # A tibble: 3 × 4
    ##   bill_length_mm bill_depth_mm flipper_length_mm body_mass_g
    ##            <dbl>         <dbl>             <int>       <int>
    ## 1           39.1          18.7               181        3750
    ## 2           39.5          17.4               186        3800
    ## 3           40.3          18                 195        3250

``` r
## two alternative ways to write the same thing:
# penguins |>
#     select( ends_with("_mm") & matches("_length") )  |>
#     head(3)
penguins |>
    select(intersect (ends_with("_mm"), matches("_length") )) |> 
    head(3)
```

    ## # A tibble: 3 × 2
    ##   bill_length_mm flipper_length_mm
    ##            <dbl>             <int>
    ## 1           39.1               181
    ## 2           39.5               186
    ## 3           40.3               195

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
    ) |> 
    head(3)
```

    ## # A tibble: 3 × 5
    ##   bill_length_mm bill_depth_mm flipper_length_mm body_mass_g  year
    ##            <dbl>         <dbl>             <int>       <int> <int>
    ## 1           39.1          18.7               181        3750  2007
    ## 2           39.5          17.4               186        3800  2007
    ## 3           40.3          18                 195        3250  2007

# summarize(across())

Use this instead to get column means, medians etc

``` r
penguins |> 
    summarise(across(where(is.numeric), 
                     ~ mean(.x, na.rm=TRUE)))
```

    ## # A tibble: 1 × 5
    ##   bill_length_mm bill_depth_mm flipper_length_mm body_mass_g  year
    ##            <dbl>         <dbl>             <dbl>       <dbl> <dbl>
    ## 1           43.9          17.2              201.       4202. 2008.
