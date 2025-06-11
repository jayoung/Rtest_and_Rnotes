Rmarkdown tips and tricks
================
Janet Young

2025-06-11

# Bulleted lists

Rules for bulleted lists:

- there should be an EMPTY LINE before the list starts
- top-level list items start with an `*`, `-`, or `+` and don’t need an
  indent. I think I’ll use `*` to keep things uniform.
  - next-level list items are indented. Can also start with `*`, `-`, or
    `+` but I think I’ll use `+` to keep things uniform.
    - third-level indent
- you can add more items

Rules for numbered lists:

1.  there should be an EMPTY LINE before the list starts
2.  top-level list items start with an asterisk and don’t need an indent
    - next-level list items start with a letter
    - second next-level item here
    - it doesn’t seem easy to have a/b/c tags for second level. Maybe
      [this link](https://pandoc.org/MANUAL.html#ordered-lists) explains
      how, but it’s unclear.
3.  you can add more items

Without the blank line before the list, you might be able to make it
work somehow, but you need lots of extra spaces at the end of pretty
much every line

Note - Rstudio’s command-I reindent option messes up bulleted lists (at
least in early 2025). It’s a bug, reported
[here](https://github.com/rstudio/rstudio/issues/13211)

Here’s a bulleted list WITHOUT the empty line and WITHOUT extra spaces -
it looks really bad: \* there should be an EMPTY LINE before the list
starts \* top-level list items start with an `*`, `-`, or `+` and don’t
need an indent. I think I’ll use `*` to keep things uniform. +
next-level list items are indented. Can also start with `*`, `-`, or `+`
but I think I’ll use `+` to keep things uniform. - third-level indent \*
you can add more items

Here’s a bulleted list WITHOUT the empty line but WITH extra spaces. It
still looks bad!  
\* there should be an EMPTY LINE before the list starts.  
\* top-level list items start with an `*`, `-`, or `+` and don’t need an
indent. I think I’ll use `*` to keep things uniform.  
+ next-level list items are indented. Can also start with `*`, `-`, or
`+` but I think I’ll use `+` to keep things uniform.  
- third-level indent.  
\* you can add more items.

# Controlling chunk output

If we are showing code in an Rmd document, this is a good header to use
for chunks that load libraries, to suppress any output but still show
the code:

```` ```{r, echo=TRUE, results='hide', warning=FALSE, message=FALSE``` ````

And here’s a chunk where I use that header:

``` r
library(tidyverse)
```

# Finished

``` r
sessionInfo()
```

    ## R version 4.5.0 (2025-04-11)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sequoia 15.5
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
    ##  [1] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
    ##  [5] purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1   
    ##  [9] ggplot2_3.5.2   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6       compiler_4.5.0     tidyselect_1.2.1   scales_1.4.0      
    ##  [5] yaml_2.3.10        fastmap_1.2.0      R6_2.6.1           generics_0.1.4    
    ##  [9] knitr_1.50         pillar_1.10.2      RColorBrewer_1.1-3 tzdb_0.5.0        
    ## [13] rlang_1.1.6        stringi_1.8.7      xfun_0.52          timechange_0.3.0  
    ## [17] cli_3.6.5          withr_3.0.2        magrittr_2.0.3     digest_0.6.37     
    ## [21] grid_4.5.0         rstudioapi_0.17.1  hms_1.1.3          lifecycle_1.0.4   
    ## [25] vctrs_0.6.5        evaluate_1.0.3     glue_1.8.0         farver_2.1.2      
    ## [29] rmarkdown_2.29     tools_4.5.0        pkgconfig_2.0.3    htmltools_0.5.8.1
