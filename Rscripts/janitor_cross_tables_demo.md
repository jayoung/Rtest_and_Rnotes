janitor package tabulation demo
================
Janet Young

2025-12-10

# janitor tabulation

<https://cran.r-project.org/web/packages/janitor/vignettes/tabyls.html>

<https://hutchdatascience.org/data_snacks/r_snacks/janitor.html>

``` r
### read example breakfast cereal data
cereals_raw <- readr::read_csv(here("Rscripts/janitor_cross_tables_demo_data/cereal.csv"),
                               show_col_types = FALSE)


### clean data
manu_labels <- c("American Home"="A",
                 "General Mills"="G",
                 "Kelloggs"="K",
                 "Nabisco" = "N",
                 "Post" = "P",
                 "Quaker Oats" = "Q", 
                 "Ralston Purina" = "R")
cereals <- cereals_raw |> 
    dplyr::rename(manufacturer=mfr) |>
    janitor::clean_names() |>
    mutate(shelf = factor(shelf, ordered=TRUE)) |>
    mutate(type=str_replace(type, "H","hot")) |>
    mutate(type=str_replace(type, "C","cold")) |>
    mutate(across(c("manufacturer", "type"), as.factor)) |>
    mutate(manufacturer = forcats::fct_recode(manufacturer, !!!manu_labels))

## show first 5 rows
cereals |> 
    slice_head(n=5)|>
    knitr::kable() |> 
    kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

name
</th>

<th style="text-align:left;">

manufacturer
</th>

<th style="text-align:left;">

type
</th>

<th style="text-align:right;">

calories
</th>

<th style="text-align:right;">

protein
</th>

<th style="text-align:right;">

fat
</th>

<th style="text-align:right;">

sodium
</th>

<th style="text-align:right;">

fiber
</th>

<th style="text-align:right;">

carbo
</th>

<th style="text-align:right;">

sugars
</th>

<th style="text-align:right;">

potass
</th>

<th style="text-align:right;">

vitamins
</th>

<th style="text-align:left;">

shelf
</th>

<th style="text-align:right;">

weight
</th>

<th style="text-align:right;">

cups
</th>

<th style="text-align:right;">

rating
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

100% Bran
</td>

<td style="text-align:left;">

Nabisco
</td>

<td style="text-align:left;">

cold
</td>

<td style="text-align:right;">

70
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

130
</td>

<td style="text-align:right;">

10
</td>

<td style="text-align:right;">

5
</td>

<td style="text-align:right;">

6
</td>

<td style="text-align:right;">

280
</td>

<td style="text-align:right;">

25
</td>

<td style="text-align:left;">

3
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

0.33
</td>

<td style="text-align:right;">

68.40297
</td>

</tr>

<tr>

<td style="text-align:left;">

100% Natural Bran
</td>

<td style="text-align:left;">

Quaker Oats
</td>

<td style="text-align:left;">

cold
</td>

<td style="text-align:right;">

120
</td>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

5
</td>

<td style="text-align:right;">

15
</td>

<td style="text-align:right;">

2
</td>

<td style="text-align:right;">

8
</td>

<td style="text-align:right;">

8
</td>

<td style="text-align:right;">

135
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:left;">

3
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1.00
</td>

<td style="text-align:right;">

33.98368
</td>

</tr>

<tr>

<td style="text-align:left;">

All-Bran
</td>

<td style="text-align:left;">

Kelloggs
</td>

<td style="text-align:left;">

cold
</td>

<td style="text-align:right;">

70
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

260
</td>

<td style="text-align:right;">

9
</td>

<td style="text-align:right;">

7
</td>

<td style="text-align:right;">

5
</td>

<td style="text-align:right;">

320
</td>

<td style="text-align:right;">

25
</td>

<td style="text-align:left;">

3
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

0.33
</td>

<td style="text-align:right;">

59.42551
</td>

</tr>

<tr>

<td style="text-align:left;">

All-Bran with Extra Fiber
</td>

<td style="text-align:left;">

Kelloggs
</td>

<td style="text-align:left;">

cold
</td>

<td style="text-align:right;">

50
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

140
</td>

<td style="text-align:right;">

14
</td>

<td style="text-align:right;">

8
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

330
</td>

<td style="text-align:right;">

25
</td>

<td style="text-align:left;">

3
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

0.50
</td>

<td style="text-align:right;">

93.70491
</td>

</tr>

<tr>

<td style="text-align:left;">

Almond Delight
</td>

<td style="text-align:left;">

Ralston Purina
</td>

<td style="text-align:left;">

cold
</td>

<td style="text-align:right;">

110
</td>

<td style="text-align:right;">

2
</td>

<td style="text-align:right;">

2
</td>

<td style="text-align:right;">

200
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

14
</td>

<td style="text-align:right;">

8
</td>

<td style="text-align:right;">

-1
</td>

<td style="text-align:right;">

25
</td>

<td style="text-align:left;">

3
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

0.75
</td>

<td style="text-align:right;">

34.38484
</td>

</tr>

</tbody>

</table>

default tabyl gives a column called “percent” that’s actually a fraction
(out of 1) rather than a percent. Argh.

``` r
cereals |> 
    tabyl(shelf) |> 
    dplyr::rename(fraction=percent) |>
    knitr::kable() |> 
    kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

shelf
</th>

<th style="text-align:right;">

n
</th>

<th style="text-align:right;">

fraction
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

1
</td>

<td style="text-align:right;">

20
</td>

<td style="text-align:right;">

0.2597403
</td>

</tr>

<tr>

<td style="text-align:left;">

2
</td>

<td style="text-align:right;">

21
</td>

<td style="text-align:right;">

0.2727273
</td>

</tr>

<tr>

<td style="text-align:left;">

3
</td>

<td style="text-align:right;">

36
</td>

<td style="text-align:right;">

0.4675325
</td>

</tr>

</tbody>

</table>

But it’s easy to reformat that as an actual %

``` r
cereals |> 
    tabyl(shelf) |> 
    adorn_pct_formatting() |>
    knitr::kable() |> 
    kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

shelf
</th>

<th style="text-align:right;">

n
</th>

<th style="text-align:left;">

percent
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

1
</td>

<td style="text-align:right;">

20
</td>

<td style="text-align:left;">

26.0%
</td>

</tr>

<tr>

<td style="text-align:left;">

2
</td>

<td style="text-align:right;">

21
</td>

<td style="text-align:left;">

27.3%
</td>

</tr>

<tr>

<td style="text-align:left;">

3
</td>

<td style="text-align:right;">

36
</td>

<td style="text-align:left;">

46.8%
</td>

</tr>

</tbody>

</table>

``` r
cereals 
```

    ## # A tibble: 77 × 16
    ##    name      manufacturer type  calories protein   fat sodium fiber carbo sugars
    ##    <chr>     <fct>        <fct>    <dbl>   <dbl> <dbl>  <dbl> <dbl> <dbl>  <dbl>
    ##  1 100% Bran Nabisco      cold        70       4     1    130  10     5        6
    ##  2 100% Nat… Quaker Oats  cold       120       3     5     15   2     8        8
    ##  3 All-Bran  Kelloggs     cold        70       4     1    260   9     7        5
    ##  4 All-Bran… Kelloggs     cold        50       4     0    140  14     8        0
    ##  5 Almond D… Ralston Pur… cold       110       2     2    200   1    14        8
    ##  6 Apple Ci… General Mil… cold       110       2     2    180   1.5  10.5     10
    ##  7 Apple Ja… Kelloggs     cold       110       2     0    125   1    11       14
    ##  8 Basic 4   General Mil… cold       130       3     2    210   2    18        8
    ##  9 Bran Chex Ralston Pur… cold        90       2     1    200   4    15        6
    ## 10 Bran Fla… Post         cold        90       3     0    210   5    13        5
    ## # ℹ 67 more rows
    ## # ℹ 6 more variables: potass <dbl>, vitamins <dbl>, shelf <ord>, weight <dbl>,
    ## #   cups <dbl>, rating <dbl>

A more complex example - a two way table

``` r
cereals |>
    janitor::tabyl(manufacturer, type) |>
    janitor::adorn_percentages(denominator = "row") |>
    janitor::adorn_pct_formatting() |>
    janitor::adorn_ns() |>
    knitr::kable() |> 
    kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

manufacturer
</th>

<th style="text-align:left;">

cold
</th>

<th style="text-align:left;">

hot
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

American Home
</td>

<td style="text-align:left;">

0.0% (0)
</td>

<td style="text-align:left;">

100.0% (1)
</td>

</tr>

<tr>

<td style="text-align:left;">

General Mills
</td>

<td style="text-align:left;">

100.0% (22)
</td>

<td style="text-align:left;">

0.0% (0)
</td>

</tr>

<tr>

<td style="text-align:left;">

Kelloggs
</td>

<td style="text-align:left;">

100.0% (23)
</td>

<td style="text-align:left;">

0.0% (0)
</td>

</tr>

<tr>

<td style="text-align:left;">

Nabisco
</td>

<td style="text-align:left;">

83.3% (5)
</td>

<td style="text-align:left;">

16.7% (1)
</td>

</tr>

<tr>

<td style="text-align:left;">

Post
</td>

<td style="text-align:left;">

100.0% (9)
</td>

<td style="text-align:left;">

0.0% (0)
</td>

</tr>

<tr>

<td style="text-align:left;">

Quaker Oats
</td>

<td style="text-align:left;">

87.5% (7)
</td>

<td style="text-align:left;">

12.5% (1)
</td>

</tr>

<tr>

<td style="text-align:left;">

Ralston Purina
</td>

<td style="text-align:left;">

100.0% (8)
</td>

<td style="text-align:left;">

0.0% (0)
</td>

</tr>

</tbody>

</table>

# Fisnihed

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Tahoe 26.1
    ## 
    ## Matrix products: default
    ## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] janitor_2.2.1    here_1.0.2       kableExtra_1.4.0 lubridate_1.9.4 
    ##  [5] forcats_1.0.0    stringr_1.5.2    dplyr_1.1.4      purrr_1.1.0     
    ##  [9] readr_2.1.5      tidyr_1.3.1      tibble_3.3.0     ggplot2_3.5.2   
    ## [13] tidyverse_2.0.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.6         generics_0.1.4     xml2_1.4.0         stringi_1.8.7     
    ##  [5] hms_1.1.3          digest_0.6.37      magrittr_2.0.4     evaluate_1.0.5    
    ##  [9] grid_4.5.2         timechange_0.3.0   RColorBrewer_1.1-3 fastmap_1.2.0     
    ## [13] rprojroot_2.1.1    viridisLite_0.4.2  scales_1.4.0       textshaping_1.0.3 
    ## [17] cli_3.6.5          crayon_1.5.3       rlang_1.1.6        bit64_4.6.0-1     
    ## [21] withr_3.0.2        yaml_2.3.10        parallel_4.5.2     tools_4.5.2       
    ## [25] tzdb_0.5.0         vctrs_0.6.5        R6_2.6.1           lifecycle_1.0.4   
    ## [29] snakecase_0.11.1   bit_4.6.0          vroom_1.6.6        pkgconfig_2.0.3   
    ## [33] pillar_1.11.1      gtable_0.3.6       glue_1.8.0         systemfonts_1.3.1 
    ## [37] xfun_0.53          tidyselect_1.2.1   rstudioapi_0.17.1  knitr_1.50        
    ## [41] farver_2.1.2       htmltools_0.5.8.1  rmarkdown_2.29     svglite_2.2.1     
    ## [45] compiler_4.5.2
