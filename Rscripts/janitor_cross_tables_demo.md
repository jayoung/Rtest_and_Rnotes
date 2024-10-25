janitor package tabulation demo
================
Janet Young

2024-10-16

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

head(cereals)
```

    ## # A tibble: 6 × 16
    ##   name       manufacturer type  calories protein   fat sodium fiber carbo sugars
    ##   <chr>      <fct>        <fct>    <dbl>   <dbl> <dbl>  <dbl> <dbl> <dbl>  <dbl>
    ## 1 100% Bran  Nabisco      cold        70       4     1    130  10     5        6
    ## 2 100% Natu… Quaker Oats  cold       120       3     5     15   2     8        8
    ## 3 All-Bran   Kelloggs     cold        70       4     1    260   9     7        5
    ## 4 All-Bran … Kelloggs     cold        50       4     0    140  14     8        0
    ## 5 Almond De… Ralston Pur… cold       110       2     2    200   1    14        8
    ## 6 Apple Cin… General Mil… cold       110       2     2    180   1.5  10.5     10
    ## # ℹ 6 more variables: potass <dbl>, vitamins <dbl>, shelf <ord>, weight <dbl>,
    ## #   cups <dbl>, rating <dbl>

default tabyl gives a column called “percent” that’s actually a fraction
(out of 1) rather than a percent. Argh.

``` r
cereals %>% 
    tabyl(shelf) %>% 
    dplyr::rename(fraction=percent)
```

    ##  shelf  n  fraction
    ##      1 20 0.2597403
    ##      2 21 0.2727273
    ##      3 36 0.4675325

But it’s easy to reformat that as an actual %

``` r
cereals %>% 
    tabyl(shelf) %>% 
      adorn_pct_formatting()
```

    ##  shelf  n percent
    ##      1 20   26.0%
    ##      2 21   27.3%
    ##      3 36   46.8%

more complex example

``` r
cereals |>
    janitor::tabyl(manufacturer, type) |>
    janitor::adorn_percentages(denominator = "row") |>
    janitor::adorn_pct_formatting() |>
    janitor::adorn_ns() %>%
    knitr::kable() %>% 
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
