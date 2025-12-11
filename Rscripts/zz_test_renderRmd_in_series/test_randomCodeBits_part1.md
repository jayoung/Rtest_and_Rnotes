test\_randomCodeBits\_part1
================
Janet Young

2025-11-07

This is the first in a set of two linked Rmd scripts, to help me test my
`render_Rmd_series.perl` script.

The render\_Rmd\_series.perl script is set up so that it will not run
downstream scripts if anything in the chain fails.

Make sure to check the resulting

Here’s the command I use to run the series of scripts:

    render_Rmd_series.perl test_randomCodeBits_part1.Rmd test_randomCodeBits_part2.Rmd

To remove all output:

    rm -r iris_plus_groups* z* *Rerr.txt *Rout.txt *_files *.md slurm-*out

This script (part1) creates a fake dataset, which we will make a plot of
in part2.

``` r
iris_plus_groups <- iris |> 
    as_tibble() |> 
    mutate(group = sample(1:2, size=nrow(iris), replace=TRUE)) |> 
    mutate(group=paste0("group_", group)) 
```

Here’s some code that will produce an error, if we were to change `eval`
to `TRUE.` The render\_Rmd\_series.perl script is set up so that it will
not run downstream scripts if anything in the chain fails.

``` r
xxxx
```

# Finished

``` r
save(iris_plus_groups, 
     file=here("Rscripts/zz_test_renderRmd_in_series/iris_plus_groups.Rdata"))
```

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 18.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: FlexiBLAS OPENBLAS;  LAPACK version 3.11.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] here_1.0.1      lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1  
    ##  [5] dplyr_1.1.4     purrr_1.0.4     readr_2.1.5     tidyr_1.3.1    
    ##  [9] tibble_3.3.0    ggplot2_3.5.1   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6      compiler_4.4.1    tidyselect_1.2.1  scales_1.3.0     
    ##  [5] yaml_2.3.8        fastmap_1.2.0     R6_2.6.1          generics_0.1.4   
    ##  [9] knitr_1.50        rprojroot_2.0.4   munsell_0.5.1     pillar_1.11.0    
    ## [13] tzdb_0.4.0        rlang_1.1.5       stringi_1.8.7     xfun_0.53        
    ## [17] timechange_0.3.0  cli_3.6.3         withr_3.0.2       magrittr_2.0.3   
    ## [21] digest_0.6.36     grid_4.4.1        hms_1.1.3         lifecycle_1.0.4  
    ## [25] vctrs_0.6.5       evaluate_0.24.0   glue_1.8.0        colorspace_2.1-0 
    ## [29] rmarkdown_2.27    tools_4.4.1       pkgconfig_2.0.3   htmltools_0.5.8.1
