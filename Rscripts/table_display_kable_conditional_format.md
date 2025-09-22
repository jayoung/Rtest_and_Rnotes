table_display_flextable_demo
================
Janet Young

2025-09-22

# Goal

Get conditional formatting using kable to show up in a github_document
rendering.

Problem - it shows up nicely in the Rstudio interface and in the knit
window that pops up, but not when we render it on the github site.

## mtcars example

Set up very simple data table:

``` r
my_dat <- mtcars %>% 
    slice_head(n=6) %>% 
    select(am, carb, gear, mpg, drat)
```

Show without conditional formatting:

``` r
my_dat %>% 
    kable(caption="Example kable table without conditional formatting") %>% 
    kable_styling(full_width = FALSE)
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Example kable table without conditional formatting
</caption>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

am
</th>

<th style="text-align:right;">

carb
</th>

<th style="text-align:right;">

gear
</th>

<th style="text-align:right;">

mpg
</th>

<th style="text-align:right;">

drat
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Mazda RX4
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

21.0
</td>

<td style="text-align:right;">

3.90
</td>

</tr>

<tr>

<td style="text-align:left;">

Mazda RX4 Wag
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

21.0
</td>

<td style="text-align:right;">

3.90
</td>

</tr>

<tr>

<td style="text-align:left;">

Datsun 710
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

22.8
</td>

<td style="text-align:right;">

3.85
</td>

</tr>

<tr>

<td style="text-align:left;">

Hornet 4 Drive
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

21.4
</td>

<td style="text-align:right;">

3.08
</td>

</tr>

<tr>

<td style="text-align:left;">

Hornet Sportabout
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

2
</td>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

18.7
</td>

<td style="text-align:right;">

3.15
</td>

</tr>

<tr>

<td style="text-align:left;">

Valiant
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

18.1
</td>

<td style="text-align:right;">

2.76
</td>

</tr>

</tbody>

</table>

Show without conditional formatting adding html options:

``` r
my_dat %>% 
    kable("html",
        # format = "latex", 
        # booktabs = TRUE,
          caption="Example kable table without conditional formatting add html option", escape = FALSE) %>% 
    kable_styling(full_width = FALSE)
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Example kable table without conditional formatting add html option
</caption>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

am
</th>

<th style="text-align:right;">

carb
</th>

<th style="text-align:right;">

gear
</th>

<th style="text-align:right;">

mpg
</th>

<th style="text-align:right;">

drat
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Mazda RX4
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

21.0
</td>

<td style="text-align:right;">

3.90
</td>

</tr>

<tr>

<td style="text-align:left;">

Mazda RX4 Wag
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

21.0
</td>

<td style="text-align:right;">

3.90
</td>

</tr>

<tr>

<td style="text-align:left;">

Datsun 710
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

22.8
</td>

<td style="text-align:right;">

3.85
</td>

</tr>

<tr>

<td style="text-align:left;">

Hornet 4 Drive
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

21.4
</td>

<td style="text-align:right;">

3.08
</td>

</tr>

<tr>

<td style="text-align:left;">

Hornet Sportabout
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

2
</td>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

18.7
</td>

<td style="text-align:right;">

3.15
</td>

</tr>

<tr>

<td style="text-align:left;">

Valiant
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

18.1
</td>

<td style="text-align:right;">

2.76
</td>

</tr>

</tbody>

</table>

With conditional formatting:

``` r
my_dat %>% 
    mutate(drat = cell_spec(drat, 
                            format = "html", 
                            background=ifelse(drat > 3, "red", "green"))) %>% 
    kable("html",
        # format = "latex", 
        # booktabs = TRUE,
          caption="Example kable table with conditional formatting", escape = FALSE) %>% 
    kable_styling(full_width = FALSE)
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Example kable table with conditional formatting
</caption>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

am
</th>

<th style="text-align:right;">

carb
</th>

<th style="text-align:right;">

gear
</th>

<th style="text-align:right;">

mpg
</th>

<th style="text-align:left;">

drat
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Mazda RX4
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

21.0
</td>

<td style="text-align:left;">

<span style="     border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: red !important;">3.9</span>
</td>

</tr>

<tr>

<td style="text-align:left;">

Mazda RX4 Wag
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

21.0
</td>

<td style="text-align:left;">

<span style="     border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: red !important;">3.9</span>
</td>

</tr>

<tr>

<td style="text-align:left;">

Datsun 710
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

22.8
</td>

<td style="text-align:left;">

<span style="     border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: red !important;">3.85</span>
</td>

</tr>

<tr>

<td style="text-align:left;">

Hornet 4 Drive
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

21.4
</td>

<td style="text-align:left;">

<span style="     border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: red !important;">3.08</span>
</td>

</tr>

<tr>

<td style="text-align:left;">

Hornet Sportabout
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

2
</td>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

18.7
</td>

<td style="text-align:left;">

<span style="     border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: red !important;">3.15</span>
</td>

</tr>

<tr>

<td style="text-align:left;">

Valiant
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

18.1
</td>

<td style="text-align:left;">

<span style="     border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: green !important;">2.76</span>
</td>

</tr>

</tbody>

</table>

## kableExtra cookbook example

Example from here
<https://sharlagelfand.github.io/kableExtra-cookbook/how-to.html#>

``` r
data_noun <- comma_survey %>%
    filter(!is.na(care_data)) %>%
    select(education, care_data)


data_noun_percent <- data_noun %>%
    mutate(education = fct_explicit_na(education, na_level = "Education unknown")) %>%
    group_by(education, care_data) %>%
    summarise(n = n()) %>%
    mutate(
        prop = n / sum(n),
        percent = percent(prop)
    ) %>%
    ungroup() %>%
    select(-n)
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `education = fct_explicit_na(education, na_level = "Education
    ##   unknown")`.
    ## Caused by warning:
    ## ! `fct_explicit_na()` was deprecated in forcats 1.0.0.
    ## ℹ Please use `fct_na_value_to_level()` instead.

    ## `summarise()` has grouped output by 'education'. You can override using the
    ## `.groups` argument.

``` r
# data_noun_percent
data_noun_percent_wide <- data_noun_percent %>%
    select(-prop) %>%
    pivot_wider(
        names_from = "care_data",
        values_from = "percent"
    )

data_noun_percent_wide
```

    ## # A tibble: 6 × 5
    ##   education                        `Not at all` `Not much` Some   `A lot`
    ##   <ord>                            <chr>        <chr>      <chr>  <chr>  
    ## 1 Less than high school degree     27.3%        45.5%      18.2%  9.1%   
    ## 2 High school degree               28%          34%        30%    8%     
    ## 3 Some college or Associate degree 20.0%        40.3%      29.8%  9.8%   
    ## 4 Bachelor degree                  18.90%       35.47%     35.76% 9.88%  
    ## 5 Graduate degree                  13.0%        37.3%      31.9%  17.8%  
    ## 6 Education unknown                18.5%        30.8%      32.3%  18.5%

``` r
data_noun_percent_wide %>%
    kable(
        col.names = c("Education", "Not at all", "Not much", "Some", "A lot"),
        align = c("lrrrr")
    ) %>%
    add_header_above(header = c(" " = 1, "How much, if at all, do you care about the debate over the use of the word 'data' as a singular or plural noun?" = 4))
```

<table>

<thead>

<tr>

<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">

</th>

<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

How much, if at all, do you care about the debate over the use of the
word ‘data’ as a singular or plural noun?

</div>

</th>

</tr>

<tr>

<th style="text-align:left;">

Education
</th>

<th style="text-align:right;">

Not at all
</th>

<th style="text-align:right;">

Not much
</th>

<th style="text-align:right;">

Some
</th>

<th style="text-align:right;">

A lot
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Less than high school degree
</td>

<td style="text-align:right;">

27.3%
</td>

<td style="text-align:right;">

45.5%
</td>

<td style="text-align:right;">

18.2%
</td>

<td style="text-align:right;">

9.1%
</td>

</tr>

<tr>

<td style="text-align:left;">

High school degree
</td>

<td style="text-align:right;">

28%
</td>

<td style="text-align:right;">

34%
</td>

<td style="text-align:right;">

30%
</td>

<td style="text-align:right;">

8%
</td>

</tr>

<tr>

<td style="text-align:left;">

Some college or Associate degree
</td>

<td style="text-align:right;">

20.0%
</td>

<td style="text-align:right;">

40.3%
</td>

<td style="text-align:right;">

29.8%
</td>

<td style="text-align:right;">

9.8%
</td>

</tr>

<tr>

<td style="text-align:left;">

Bachelor degree
</td>

<td style="text-align:right;">

18.90%
</td>

<td style="text-align:right;">

35.47%
</td>

<td style="text-align:right;">

35.76%
</td>

<td style="text-align:right;">

9.88%
</td>

</tr>

<tr>

<td style="text-align:left;">

Graduate degree
</td>

<td style="text-align:right;">

13.0%
</td>

<td style="text-align:right;">

37.3%
</td>

<td style="text-align:right;">

31.9%
</td>

<td style="text-align:right;">

17.8%
</td>

</tr>

<tr>

<td style="text-align:left;">

Education unknown
</td>

<td style="text-align:right;">

18.5%
</td>

<td style="text-align:right;">

30.8%
</td>

<td style="text-align:right;">

32.3%
</td>

<td style="text-align:right;">

18.5%
</td>

</tr>

</tbody>

</table>

``` r
k <- data_noun_percent_wide %>%
    kable(
        col.names = c(
            "Education",
            names(data_noun_percent_wide)[-1]
        ),
        align = c("l", rep("r", ncol(data_noun_percent_wide) - 1))
    )

k
```

| Education                        | Not at all | Not much |   Some | A lot |
|:---------------------------------|-----------:|---------:|-------:|------:|
| Less than high school degree     |      27.3% |    45.5% |  18.2% |  9.1% |
| High school degree               |        28% |      34% |    30% |    8% |
| Some college or Associate degree |      20.0% |    40.3% |  29.8% |  9.8% |
| Bachelor degree                  |     18.90% |   35.47% | 35.76% | 9.88% |
| Graduate degree                  |      13.0% |    37.3% |  31.9% | 17.8% |
| Education unknown                |      18.5% |    30.8% |  32.3% | 18.5% |

``` r
## this makes the highest proportion in each row bold

data_noun_bold_highest_prop <- data_noun_percent %>%
    group_by(education) %>%
    mutate(highest_prop = max(prop)) %>%
    mutate(percent = cell_spec(percent, format = "html", bold = prop == highest_prop))

# data_noun_bold_highest_prop

data_noun_bold_highest_prop_wide <- data_noun_bold_highest_prop %>%
    select(-prop, -highest_prop) %>%
    pivot_wider(
        names_from = "care_data",
        values_from = "percent"
    )

k_header <- c(1, ncol(data_noun_bold_highest_prop_wide) - 1)
names(k_header) <- c(" ", "How much, if at all, do you care about the debate over the use of the word 'data' as a singular or plural noun?")

data_noun_bold_highest_prop_wide %>%
    kable("html",
          col.names = c("Education", names(data_noun_bold_highest_prop_wide)[-1]),
          align = c("l", rep("r", ncol(data_noun_bold_highest_prop_wide) - 1)),
          escape = FALSE
    ) %>%
    add_header_above(header = k_header)
```

<table>

<thead>

<tr>

<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">

</th>

<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

How much, if at all, do you care about the debate over the use of the
word ‘data’ as a singular or plural noun?

</div>

</th>

</tr>

<tr>

<th style="text-align:left;">

Education
</th>

<th style="text-align:right;">

Not at all
</th>

<th style="text-align:right;">

Not much
</th>

<th style="text-align:right;">

Some
</th>

<th style="text-align:right;">

A lot
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Less than high school degree
</td>

<td style="text-align:right;">

<span style="     ">27.3%</span>
</td>

<td style="text-align:right;">

<span style=" font-weight: bold;    ">45.5%</span>
</td>

<td style="text-align:right;">

<span style="     ">18.2%</span>
</td>

<td style="text-align:right;">

<span style="     ">9.1%</span>
</td>

</tr>

<tr>

<td style="text-align:left;">

High school degree
</td>

<td style="text-align:right;">

<span style="     ">28%</span>
</td>

<td style="text-align:right;">

<span style=" font-weight: bold;    ">34%</span>
</td>

<td style="text-align:right;">

<span style="     ">30%</span>
</td>

<td style="text-align:right;">

<span style="     ">8%</span>
</td>

</tr>

<tr>

<td style="text-align:left;">

Some college or Associate degree
</td>

<td style="text-align:right;">

<span style="     ">20.0%</span>
</td>

<td style="text-align:right;">

<span style=" font-weight: bold;    ">40.3%</span>
</td>

<td style="text-align:right;">

<span style="     ">29.8%</span>
</td>

<td style="text-align:right;">

<span style="     ">9.8%</span>
</td>

</tr>

<tr>

<td style="text-align:left;">

Bachelor degree
</td>

<td style="text-align:right;">

<span style="     ">18.90%</span>
</td>

<td style="text-align:right;">

<span style="     ">35.47%</span>
</td>

<td style="text-align:right;">

<span style=" font-weight: bold;    ">35.76%</span>
</td>

<td style="text-align:right;">

<span style="     ">9.88%</span>
</td>

</tr>

<tr>

<td style="text-align:left;">

Graduate degree
</td>

<td style="text-align:right;">

<span style="     ">13.0%</span>
</td>

<td style="text-align:right;">

<span style=" font-weight: bold;    ">37.3%</span>
</td>

<td style="text-align:right;">

<span style="     ">31.9%</span>
</td>

<td style="text-align:right;">

<span style="     ">17.8%</span>
</td>

</tr>

<tr>

<td style="text-align:left;">

Education unknown
</td>

<td style="text-align:right;">

<span style="     ">18.5%</span>
</td>

<td style="text-align:right;">

<span style="     ">30.8%</span>
</td>

<td style="text-align:right;">

<span style=" font-weight: bold;    ">32.3%</span>
</td>

<td style="text-align:right;">

<span style="     ">18.5%</span>
</td>

</tr>

</tbody>

</table>

# Finished

``` r
sessionInfo()
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sequoia 15.6.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
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
    ##  [1] scales_1.4.0          fivethirtyeight_0.6.2 kableExtra_1.4.0     
    ##  [4] here_1.0.2            lubridate_1.9.4       forcats_1.0.0        
    ##  [7] stringr_1.5.2         dplyr_1.1.4           purrr_1.1.0          
    ## [10] readr_2.1.5           tidyr_1.3.1           tibble_3.3.0         
    ## [13] ggplot2_4.0.0         tidyverse_2.0.0      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6       compiler_4.5.1     tidyselect_1.2.1   xml2_1.4.0        
    ##  [5] textshaping_1.0.3  systemfonts_1.2.3  yaml_2.3.10        fastmap_1.2.0     
    ##  [9] R6_2.6.1           generics_0.1.4     knitr_1.50         rprojroot_2.1.1   
    ## [13] svglite_2.2.1      pillar_1.11.1      RColorBrewer_1.1-3 tzdb_0.5.0        
    ## [17] rlang_1.1.6        utf8_1.2.6         stringi_1.8.7      xfun_0.53         
    ## [21] S7_0.2.0           viridisLite_0.4.2  timechange_0.3.0   cli_3.6.5         
    ## [25] withr_3.0.2        magrittr_2.0.4     digest_0.6.37      grid_4.5.1        
    ## [29] rstudioapi_0.17.1  hms_1.1.3          lifecycle_1.0.4    vctrs_0.6.5       
    ## [33] evaluate_1.0.5     glue_1.8.0         farver_2.1.2       rmarkdown_2.29    
    ## [37] tools_4.5.1        pkgconfig_2.0.3    htmltools_0.5.8.1
