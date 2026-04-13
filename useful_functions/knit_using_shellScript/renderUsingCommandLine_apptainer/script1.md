script1
================
Janet Young

2026-04-13

``` r
mtcars_efficient <- mtcars[which(mtcars$mpg >= 30),]
```

This script filters the `mtcars` dataset (contains 32 rows) to keep only
cars where mpg is \>= 30. There are 4 rows after this filtering.

Save `mtcars_efficient` for use in script2.

``` r
# it is saved in the same location as the script
save(mtcars_efficient, file="mtcars_efficient.Rdata")
```

Show R and package version information

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Etc/UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.5.2    fastmap_1.2.0     cli_3.6.3         ragg_1.5.0       
    ##  [5] htmltools_0.5.8.1 tools_4.5.2       otel_0.2.0        yaml_2.3.10      
    ##  [9] rmarkdown_2.30    knitr_1.51        digest_0.6.37     xfun_0.57        
    ## [13] textshaping_1.0.4 lifecycle_1.0.5   systemfonts_1.3.1 rlang_1.1.5      
    ## [17] evaluate_1.0.5
